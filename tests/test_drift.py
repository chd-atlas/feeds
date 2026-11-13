"""The weekly drift check: comparison, the Orphadata header probe, and failure handling.

No network and no hvantk: the hvantk runner, the byte fetcher and ``gh`` are injected.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

_PATH = Path(__file__).resolve().parent.parent / "scripts" / "drift.py"
_spec = importlib.util.spec_from_file_location("feeds_drift_script", _PATH)
assert _spec is not None and _spec.loader is not None
drift = importlib.util.module_from_spec(_spec)
sys.modules["feeds_drift_script"] = drift
_spec.loader.exec_module(drift)

_JDBOR = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<JDBOR date="2026-06-23 07:57:18" version="1.3.42 / 4.1.8 [2025-03-03]" copyright="x">\n'
)


def _hvantk_observed(**over: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "probe_version": 2,
        "source_version": None,
        "headers": {"f.csv": ["a", "b"]},
        "checksums": {"f.csv": "aaaa"},
        "extras": {"content_length": "10"},
        "informational": {"last_modified": "Mon"},
        "fetched_at": "2026-11-13T00:00:00+00:00",
    }
    return base | over


def test_a_changed_fingerprint_is_reported_by_key() -> None:
    base = _hvantk_observed()
    assert drift.compare(base, _hvantk_observed()) == []
    # fetched_at and informational change every run and are never compared
    noise = _hvantk_observed(informational={"last_modified": "Tue"}, fetched_at="later")
    assert drift.compare(base, noise) == []
    changed = _hvantk_observed(
        checksums={"f.csv": "bbbb"},
        extras={"content_length": "11"},
        headers={"f.csv": ["a", "c"]},
    )
    assert drift.compare(base, changed) == [
        "checksums/f.csv",
        "extras/content_length",
        "headers/f.csv",
    ]
    # a key present on one side only is a change
    assert drift.compare(base, {k: v for k, v in base.items() if k != "extras"}) == [
        "extras/content_length"
    ]
    assert drift.compare(base, _hvantk_observed(source_version="4563 experiments")) == [
        "source_version"
    ]


def test_an_orphadata_header_is_read_from_the_first_bytes() -> None:
    seen: list[tuple[str, int]] = []

    def fetch(url: str, nbytes: int) -> bytes:
        seen.append((url, nbytes))
        # a truncated document: the root element is never closed
        return (_JDBOR + "x" * 10_000).encode()[:nbytes]

    got = drift.orphadata_probe("https://example.org/p4.xml", fetch)
    assert got == {"version": "1.3.42 / 4.1.8 [2025-03-03]", "date": "2026-06-23 07:57:18"}
    assert seen == [("https://example.org/p4.xml", 4096)]
    for body in (b"<html>moved</html>", b'<JDBOR date="d"/>', b'<JDBOR date=" " version="1"/>'):
        with pytest.raises(drift.ProbeFailed):
            drift.orphadata_probe("u", lambda _u, _n, b=body: b)


def _write_baselines(root: Path, fingerprints: dict[str, dict[str, Any]]) -> None:
    for source_id, fp in fingerprints.items():
        folder = root / source_id
        folder.mkdir(parents=True)
        (folder / "drift.json").write_text(
            json.dumps(
                {"captured_on": "2026-11-01", "package": f"{source_id}/v1", "fingerprint": fp}
            )
        )


def _fixture(tmp_path: Path) -> tuple[Path, dict[str, dict[str, Any]]]:
    fps: dict[str, dict[str, Any]] = {
        "hgnc": drift.fingerprint_of_hvantk(_hvantk_observed()),
        "orphanet_product4": {"version": "1", "date": "d"},
    }
    _write_baselines(tmp_path, fps)
    return tmp_path, fps


def _specs() -> dict[str, Any]:
    return {
        "hgnc": drift.HvantkSource("hgnc:lookup", "hgnc/v1"),
        "orphanet_product4": drift.OrphadataSource(
            "https://example.org/p4.xml", "orphanet_product4/v1"
        ),
    }


def _hvantk_ok(observed: dict[str, Any]) -> Callable[[str, int], tuple[int, str]]:
    def run(dataset: str, timeout: int) -> tuple[int, str]:
        status = "clean"
        return 0, json.dumps([{"dataset_name": dataset, "status": status, "observed": observed}])

    return run


def _fetch_for(version: str) -> Callable[[str, int], bytes]:
    return lambda _u, n: f'<JDBOR date="d" version="{version}">'.encode()[:n]


def test_a_probe_failure_is_reported_and_never_read_as_clean(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    root, _ = _fixture(tmp_path)
    clean = _hvantk_ok(_hvantk_observed())
    args = ["--baselines", str(root)]

    def go(hvantk: Any, fetch: Any, extra: list[str] | None = None) -> int:
        return int(
            drift.main(
                [*args, *(extra or [])], hvantk=hvantk, fetch=fetch, gh=_no_gh, sources=_specs()
            )
        )

    assert go(clean, _fetch_for("1")) == 0
    assert go(clean, _fetch_for("2")) == 1  # drifted, not failed
    capsys.readouterr()

    def probe_failed(dataset: str, timeout: int) -> tuple[int, str]:
        entry = {"dataset_name": dataset, "status": "probe_failed", "probe_error": "timed out"}
        return 2, json.dumps([entry])

    # one failed probe beside a drifted one: exit 2, and the report names the failure
    assert go(probe_failed, _fetch_for("2")) == 2
    out = capsys.readouterr().out
    assert "hgnc: PROBE FAILED" in out and "timed out" in out
    assert "orphanet_product4: DRIFTED" in out

    def garbage(dataset: str, timeout: int) -> tuple[int, str]:
        return 0, "not json"

    def crashed(dataset: str, timeout: int) -> tuple[int, str]:
        return 3, ""

    def no_observed(dataset: str, timeout: int) -> tuple[int, str]:
        return 0, json.dumps([{"dataset_name": dataset, "status": "clean"}])

    def boom(_u: str, _n: int) -> bytes:
        raise OSError("network down")

    for bad in (garbage, crashed, no_observed):
        assert go(bad, _fetch_for("1")) == 2
    assert go(clean, boom) == 2
    # a source with no baseline cannot be called clean either
    (root / "hgnc" / "drift.json").unlink()
    assert go(clean, _fetch_for("1")) == 2
    assert "no baseline" in capsys.readouterr().out


def _no_gh(argv: list[str]) -> str:
    raise AssertionError(f"gh was called: {argv}")


def test_issues_are_opened_only_with_the_flag(tmp_path: Path) -> None:
    root, _ = _fixture(tmp_path)
    calls: list[list[str]] = []

    def gh(argv: list[str]) -> str:
        calls.append(argv)
        if argv[:2] == ["issue", "list"]:
            return json.dumps([{"number": 7, "title": "drift: hgnc"}])
        return ""

    changed = _hvantk_observed(checksums={"f.csv": "bbbb"})

    def run(flags: list[str]) -> int:
        return int(
            drift.main(
                ["--baselines", str(root), *flags],
                hvantk=_hvantk_ok(changed),
                fetch=_fetch_for("1"),
                gh=gh,
                sources=_specs(),
            )
        )

    assert run([]) == 1
    assert calls == []
    assert run(["--open-issues"]) == 1
    # hgnc drifted and already has an open issue: updated, not duplicated; product 4 is clean
    verbs = [c[:2] for c in calls]
    assert verbs == [["issue", "list"], ["issue", "edit"]]
    edit = calls[1]
    body = edit[edit.index("--body") + 1]
    assert edit[2] == "7" and "hgnc/v1" in body and "2026-11-01" in body
    assert "checksums/f.csv" in body and "nothing changed since" in body

    # no open issue: one is created, titled drift: <id>
    calls.clear()

    def gh_none(argv: list[str]) -> str:
        calls.append(argv)
        return "[]"

    drift.main(
        ["--baselines", str(root), "--open-issues"],
        hvantk=_hvantk_ok(changed),
        fetch=_fetch_for("1"),
        gh=gh_none,
        sources=_specs(),
    )
    create = calls[-1]
    assert create[:2] == ["issue", "create"]
    assert create[create.index("--title") + 1] == "drift: hgnc"


def test_capture_writes_baselines_only_when_every_probe_answered(tmp_path: Path) -> None:
    root, _ = _fixture(tmp_path)

    def boom(_u: str, _n: int) -> bytes:
        raise OSError("down")

    def cap(fetch: Any) -> int:
        return int(
            drift.main(
                ["--baselines", str(root), "--capture"],
                hvantk=_hvantk_ok(_hvantk_observed(checksums={"f.csv": "cccc"})),
                fetch=fetch,
                gh=_no_gh,
                sources=_specs(),
                today="2026-12-01",
            )
        )

    before = (root / "hgnc" / "drift.json").read_text()
    assert cap(boom) == 2
    assert (root / "hgnc" / "drift.json").read_text() == before
    assert cap(_fetch_for("9")) == 0
    written = json.loads((root / "hgnc" / "drift.json").read_text())
    assert written["captured_on"] == "2026-12-01"
    assert written["package"] == "hgnc/v1"  # the spec's pinned package
    assert written["fingerprint"]["checksums"] == {"f.csv": "cccc"}
    assert "fetched_at" not in written["fingerprint"]


def _stream_for(digests: dict[str, tuple[str, int]]) -> Callable[[str], tuple[str, int]]:
    return lambda url: digests[url.rsplit("/", 1)[1]]


def _files() -> Any:
    return drift.FilesSource(
        {"a.tsv": "https://example.org/x/a.tsv", "b.xml": "https://example.org/x/b.xml"},
        "files_src/v1",
    )


def test_the_expression_atlas_source_is_fingerprinted_by_its_own_input_files() -> None:
    from sources.expression_atlas_e_mtab_6814 import recipe

    spec = drift.SOURCES["expression_atlas_e_mtab_6814"]
    assert isinstance(spec, drift.FilesSource)
    assert spec.urls == {
        recipe.MATRIX_NAME: recipe.MATRIX_URL,
        recipe.CONFIGURATION_NAME: recipe.CONFIGURATION_URL,
    }
    baseline = json.loads(
        (_PATH.parent.parent / "sources" / recipe.SOURCE / "drift.json").read_text()
    )
    got = {n: v["sha256"] for n, v in baseline["fingerprint"].items()}
    assert got == {
        recipe.MATRIX_NAME: recipe.MATRIX_SHA256,
        recipe.CONFIGURATION_NAME: recipe.CONFIGURATION_SHA256,
    }
    digests = {"a.tsv": ("aa", 5), "b.xml": ("bb", 7)}
    fp = drift.probe(_files(), _unused_hvantk, _unused_fetch, _stream_for(digests))
    assert fp == {"a.tsv": {"sha256": "aa", "size": 5}, "b.xml": {"sha256": "bb", "size": 7}}
    assert drift.compare(fp, fp) == []
    digests["b.xml"] = ("cc", 7)
    changed = drift.probe(_files(), _unused_hvantk, _unused_fetch, _stream_for(digests))
    assert drift.compare(fp, changed) == ["b.xml/sha256"]

    def down(_url: str) -> tuple[str, int]:
        raise OSError("down")

    with pytest.raises(drift.ProbeFailed):
        drift.probe(_files(), _unused_hvantk, _unused_fetch, down)


def _unused_hvantk(dataset: str, timeout: int) -> tuple[int, str]:
    raise AssertionError("hvantk was called")


def _unused_fetch(url: str, nbytes: int) -> bytes:
    raise AssertionError("fetch was called")


def _run(
    root: Path,
    *,
    flags: list[str] | None = None,
    hvantk: Any = None,
    fetch: Any = None,
    gh: Any = _no_gh,
) -> int:
    return int(
        drift.main(
            ["--baselines", str(root), *(flags or [])],
            hvantk=hvantk or _hvantk_ok(_hvantk_observed()),
            fetch=fetch or _fetch_for("1"),
            gh=gh,
            sources=_specs(),
        )
    )


@pytest.mark.parametrize("fingerprint", ["text", ["a"], None, 3])
def test_a_baseline_without_a_mapping_fingerprint_fails_that_source_only(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], fingerprint: Any
) -> None:
    root, _ = _fixture(tmp_path)
    (root / "hgnc" / "drift.json").write_text(
        json.dumps({"captured_on": "d", "package": "p", "fingerprint": fingerprint})
    )
    assert _run(root) == 2
    out = capsys.readouterr().out
    assert "hgnc: PROBE FAILED" in out
    assert "orphanet_product4: clean" in out
    (root / "hgnc" / "drift.json").write_text("[1]")  # not even an object
    assert _run(root) == 2
    assert "orphanet_product4: clean" in capsys.readouterr().out


def test_a_gh_failure_is_reported_continues_and_makes_the_run_incomplete(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    root, _ = _fixture(tmp_path)
    calls: list[list[str]] = []

    def gh(argv: list[str]) -> str:
        calls.append(argv)
        raise OSError("gh not found")

    # both sources drift; the first gh failure must not stop the second from being tried
    code = _run(
        root,
        flags=["--open-issues"],
        hvantk=_hvantk_ok(_hvantk_observed(checksums={"f.csv": "bbbb"})),
        fetch=_fetch_for("2"),
        gh=gh,
    )
    out = capsys.readouterr().out
    assert code == 2
    assert len(calls) == 2
    assert "hgnc: ISSUE FAILED" in out and "orphanet_product4: ISSUE FAILED" in out
    assert "hgnc: DRIFTED" in out and "orphanet_product4: DRIFTED" in out

    def bad_json(argv: list[str]) -> str:
        return "not json"

    assert (
        _run(
            root,
            flags=["--open-issues"],
            hvantk=_hvantk_ok(_hvantk_observed(checksums={"f.csv": "bbbb"})),
            gh=bad_json,
        )
        == 2
    )


def test_an_empty_fingerprint_is_a_probe_failure_even_against_an_empty_baseline(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    root = tmp_path
    _write_baselines(root, {"hgnc": {}, "orphanet_product4": {"version": "1", "date": "d"}})
    only_noise = {
        "informational": {"last_modified": "Mon"},
        "fetched_at": "2026-11-13T00:00:00+00:00",
    }
    assert _run(root, hvantk=_hvantk_ok(only_noise)) == 2
    assert "hgnc: PROBE FAILED" in capsys.readouterr().out
    empty = drift.FilesSource({}, "e/v1")
    with pytest.raises(drift.ProbeFailed):
        drift.probe(empty, _unused_hvantk, _unused_fetch, lambda _u: ("", 0))


def test_hvantk_answering_for_another_dataset_is_a_probe_failure(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    root, _ = _fixture(tmp_path)

    def wrong(dataset: str, timeout: int) -> tuple[int, str]:
        entry = {"dataset_name": "other:thing", "status": "clean"}
        return 0, json.dumps([entry | {"observed": _hvantk_observed()}])

    def unnamed(dataset: str, timeout: int) -> tuple[int, str]:
        return 0, json.dumps([{"status": "clean", "observed": _hvantk_observed()}])

    for bad in (wrong, unnamed):
        assert _run(root, hvantk=bad) == 2
        assert "hgnc: PROBE FAILED" in capsys.readouterr().out


class _FakeResponse:
    def __init__(self, body: bytes, status: int = 200, length: int | None = None) -> None:
        self._body, self.status = body, status
        self.headers: dict[str, str] = {} if length is None else {"Content-Length": str(length)}

    def read(self, n: int = -1) -> bytes:
        chunk, self._body = self._body[:n], self._body[n:]
        return chunk

    def __enter__(self) -> _FakeResponse:
        return self

    def __exit__(self, *_a: object) -> None:
        return None


def test_a_truncated_or_partial_download_is_not_fingerprinted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import hashlib

    def serve(response: _FakeResponse) -> None:
        monkeypatch.setattr(drift.urllib.request, "urlopen", lambda *_a, **_k: response)

    serve(_FakeResponse(b"abc", 200, 3))
    assert drift.http_stream("u") == (hashlib.sha256(b"abc").hexdigest(), 3)
    serve(_FakeResponse(b"abc", 200, None))  # no Content-Length: nothing to check against
    assert drift.http_stream("u")[1] == 3
    for bad in (_FakeResponse(b"abc", 200, 5_000_000), _FakeResponse(b"abc", 206, 3)):
        serve(bad)
        with pytest.raises(OSError):
            drift.http_stream("u")
