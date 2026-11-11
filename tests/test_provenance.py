from __future__ import annotations

import json
import subprocess
from datetime import date
from pathlib import Path
from typing import Any

import pytest

from feeds.package import PackageRefused
from feeds.provenance import (
    BUILD_PATHS,
    hvantk_sidecar,
    recipe_commit,
    require_hvantk_release,
    stamp,
)

GIT = ["git", "-c", "user.name=t", "-c", "user.email=t@example.org"]


def _git(root: Path, *args: str) -> str:
    done = subprocess.run(
        [*GIT, "-C", str(root), *args], check=True, capture_output=True, text=True
    )
    return done.stdout.strip()


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    root = tmp_path / "checkout"
    (root / "src").mkdir(parents=True)
    (root / "sources").mkdir()
    (root / "src" / "a.py").write_text("x = 1\n")
    (root / "pyproject.toml").write_text("[project]\n")
    (root / ".gitignore").write_text("__pycache__/\n")
    (root / "README.md").write_text("readme\n")
    _git(root.parent, "init", "-q", str(root))
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "first")
    return root


def test_a_clean_checkout_gives_its_head(repo: Path) -> None:
    assert recipe_commit(repo) == _git(repo, "rev-parse", "HEAD")
    # README and ignored files are outside what can change a package's bytes.
    (repo / "README.md").write_text("edited\n")
    (repo / "src" / "__pycache__").mkdir()
    (repo / "src" / "__pycache__" / "a.pyc").write_bytes(b"x")
    assert recipe_commit(repo) == _git(repo, "rev-parse", "HEAD")


@pytest.mark.parametrize("path", ["src/a.py", "pyproject.toml", "src/new.py", "sources/new.py"])
def test_uncommitted_recipe_code_is_refused(repo: Path, path: str) -> None:
    (repo / path).write_text("changed\n")  # modified, or untracked for a new file
    with pytest.raises(PackageRefused, match="uncommitted changes.*clean checkout"):
        recipe_commit(repo)


def test_a_staged_but_uncommitted_change_is_refused(repo: Path) -> None:
    (repo / "src" / "a.py").write_text("y = 2\n")
    _git(repo, "add", "src/a.py")
    with pytest.raises(PackageRefused, match="uncommitted"):
        recipe_commit(repo)


def test_a_directory_that_is_not_a_repository_is_refused(tmp_path: Path) -> None:
    with pytest.raises(PackageRefused, match="cannot read the recipe commit"):
        recipe_commit(tmp_path)


def test_every_build_path_is_checked() -> None:
    assert set(BUILD_PATHS) == {"src", "sources", "schemas", "pyproject.toml", "uv.lock"}


def test_stamp_sets_the_commit_and_refuses_a_different_supplied_one() -> None:
    commit = "c" * 40
    stamped, got = stamp({"host": "h"}, lambda: commit)
    assert (stamped, got) == ({"host": "h", "recipe_commit": commit}, commit)
    assert stamp({"recipe_commit": commit}, lambda: commit)[0] == {"recipe_commit": commit}
    with pytest.raises(PackageRefused, match="omit it"):
        stamp({"recipe_commit": "d" * 40}, lambda: commit)
    with pytest.raises(PackageRefused, match="40-hex"):
        stamp({}, lambda: "HEAD")
    original: dict[str, Any] = {"host": "h"}
    stamp(original, lambda: commit)
    assert original == {"host": "h"}  # the caller's dict is not edited


@pytest.mark.parametrize("value", [None, "", "  ", 4])
def test_hvantk_release_is_required(value: object) -> None:
    with pytest.raises(PackageRefused, match="hvantk_release"):
        require_hvantk_release({"hvantk_release": value})
    assert require_hvantk_release({"hvantk_release": "96ceafb4"}) == "96ceafb4"


def test_the_hvantk_sidecar_gives_the_commit_and_the_date(tmp_path: Path) -> None:
    table = tmp_path / "t.tsv.gz"
    sidecar = tmp_path / "t.tsv.gz.provenance.json"
    commit = "4" * 40
    sidecar.write_text(
        json.dumps({"builder_commit": commit, "build_timestamp": "2026-11-10T09:30:00+00:00"})
    )
    assert hvantk_sidecar(table) == (commit, date(2026, 11, 10))
    sidecar.write_text(json.dumps({"builder_commit": "hvantk 0.6.0", "build_timestamp": "x"}))
    with pytest.raises(PackageRefused, match="'hvantk 0.6.0'"):
        hvantk_sidecar(table)
    sidecar.unlink()
    with pytest.raises(PackageRefused, match="missing or unreadable"):
        hvantk_sidecar(table)
