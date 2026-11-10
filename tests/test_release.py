from __future__ import annotations

import hashlib
import subprocess
from collections.abc import Sequence
from pathlib import Path

import pytest

from feeds.package import write_package
from feeds.release import ReleaseRefused, publish
from tests.conftest import package_kwargs

COMMIT = "b" * 40


class Fake:
    def __init__(
        self,
        *,
        branches: str | None = None,
        view_error: str | None = None,
        release_dir: Path | None = None,
    ) -> None:
        self.branches, self.view_error, self.release_dir, self.calls = (
            branches,
            view_error,
            release_dir,
            [],
        )  # type: ignore[var-annotated]

    def __call__(self, command: Sequence[str]) -> str:
        self.calls.append(list(command))
        if command[:2] == ["git", "branch"]:
            return self.branches if self.branches is not None else "  origin/dev\n"
        if command[:3] == ["gh", "release", "view"]:
            if self.view_error is not None:
                raise subprocess.CalledProcessError(1, command, stderr=self.view_error)
            if self.release_dir is None:
                raise subprocess.CalledProcessError(1, command, stderr="release not found")
            return ""
        if command[:3] == ["gh", "release", "download"]:
            assert self.release_dir is not None
            dest = Path(command[command.index("-D") + 1])
            for p in self.release_dir.iterdir():
                (dest / p.name).write_bytes(p.read_bytes())
        return ""


@pytest.fixture
def dist(tmp_path: Path) -> Path:
    return write_package(tmp_path / "dist", **package_kwargs())


def _digest(directory: Path) -> dict[str, str]:
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.iterdir()}


def test_release_refuses_an_unpushed_recipe_commit_and_a_changed_tag(
    dist: Path, tmp_path: Path
) -> None:
    fake = Fake(branches="  origin/dev-old\n  origin/other\n")
    with pytest.raises(ReleaseRefused, match="origin/dev"):
        publish(dist, dry_run=False, runner=fake)
    assert not any(c[:3] == ["gh", "release", "create"] for c in fake.calls)
    assert fake.calls[0] == ["git", "fetch", "origin", "dev"]

    published = tmp_path / "published"
    published.mkdir()
    for p in dist.iterdir():
        (published / p.name).write_bytes(p.read_bytes())
    assert publish(dist, dry_run=False, runner=Fake(release_dir=published)) == "unchanged"

    (published / "burden.tsv.gz").write_bytes(b"other")
    fake = Fake(release_dir=published)
    with pytest.raises(ReleaseRefused, match="burden.tsv.gz"):
        publish(dist, dry_run=False, runner=fake)
    assert not any(c[:3] == ["gh", "release", "create"] for c in fake.calls)

    (published / "burden.tsv.gz").write_bytes((dist / "burden.tsv.gz").read_bytes())
    (published / "extra.txt").write_text("x")
    with pytest.raises(ReleaseRefused, match="extra.txt"):
        publish(dist, dry_run=False, runner=Fake(release_dir=published))


def test_a_new_tag_is_created_at_the_recipe_commit_and_a_dry_run_only_prints(dist: Path) -> None:
    dry = Fake()
    line = publish(dist, runner=dry)
    assert line.startswith("dry-run: gh release create demo/2026-11-03 ")
    assert not any(c[:3] == ["gh", "release", "create"] for c in dry.calls)

    live = Fake()
    assert publish(dist, dry_run=False, runner=live) == "demo/2026-11-03"
    create = live.calls[-1]
    assert create[:4] == ["gh", "release", "create", "demo/2026-11-03"]
    assert create[create.index("--target") + 1] == COMMIT
    assert {Path(a).name for a in create if a.startswith(str(dist))} == _digest(dist).keys()


def test_release_refuses_a_dist_with_a_stray_or_missing_file(dist: Path) -> None:
    (dist / "stray.txt").write_text("x")
    with pytest.raises(ReleaseRefused, match="stray.txt"):
        publish(dist, runner=Fake())
    (dist / "stray.txt").unlink()
    (dist / "burden.tsv.gz").unlink()
    with pytest.raises(ReleaseRefused, match="burden.tsv.gz"):
        publish(dist, runner=Fake())


def test_only_a_not_found_view_means_no_release(dist: Path) -> None:
    fake = Fake(view_error="HTTP 401: Bad credentials")
    with pytest.raises(subprocess.CalledProcessError):
        publish(dist, dry_run=False, runner=fake)
    assert not any(c[:3] == ["gh", "release", "create"] for c in fake.calls)


def test_release_refuses_a_file_that_no_longer_matches_the_manifest_even_in_a_dry_run(
    dist: Path,
) -> None:
    (dist / "burden.tsv.gz").write_bytes(b"tampered")
    fake = Fake()
    with pytest.raises(ReleaseRefused, match="sha256"):
        publish(dist, runner=fake)
    assert fake.calls == []


def test_the_dry_run_command_is_shell_quoted(dist: Path, tmp_path: Path) -> None:
    spaced = tmp_path / "with space"
    dist.rename(spaced)
    line = publish(spaced, runner=Fake())
    assert f"'{spaced}/manifest.json'" in line
