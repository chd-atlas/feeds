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
