"""Publish a built package as the GitHub release ``<package>/<version>``."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
from collections.abc import Callable, Sequence
from pathlib import Path

Runner = Callable[[Sequence[str]], str]


class ReleaseRefused(Exception):
    """The release would be unsafe or ambiguous; nothing was published."""


def _run(command: Sequence[str]) -> str:
    done = subprocess.run(command, check=True, capture_output=True, text=True)
    return done.stdout


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def publish(
    dist: Path,
    *,
    repo: str = "chd-atlas/feeds",
    dry_run: bool = True,
    runner: Runner = _run,
) -> str:
    """Return ``"unchanged"``, ``"dry-run: <command>"`` or the tag that was created."""
    manifest = json.loads((dist / "manifest.json").read_text(encoding="utf-8"))
    commit: str = manifest["built_by"]["recipe_commit"]
    tag = f"{manifest['package']}/{manifest['version']}"
    expected = {"manifest.json", *(t["file"] for t in manifest["tables"])}
    present = {p.name for p in dist.iterdir()}
    if present != expected:
        raise ReleaseRefused(
            f"dist differs from the manifest: unexpected {sorted(present - expected)}, "
            f"missing {sorted(expected - present)}"
        )
    runner(["git", "fetch", "origin", "dev"])
    branches = [
        b.strip() for b in runner(["git", "branch", "-r", "--contains", commit]).splitlines()
    ]
    if "origin/dev" not in branches:
        raise ReleaseRefused(f"recipe commit {commit} is not on origin/dev; push it first")
    try:
        runner(["gh", "release", "view", tag, "-R", repo])
        exists = True
    except subprocess.CalledProcessError:
        exists = False
    if exists:
        with tempfile.TemporaryDirectory() as tmp:
            runner(["gh", "release", "download", tag, "-R", repo, "-D", tmp])
            theirs = {p.name: _sha(p) for p in Path(tmp).iterdir()}
        ours = {name: _sha(dist / name) for name in sorted(expected)}
        differing = sorted(n for n in ours.keys() | theirs.keys() if ours.get(n) != theirs.get(n))
        if differing:
            raise ReleaseRefused(f"release {tag} exists with different assets: {differing}")
        return "unchanged"
    assets = [str(dist / name) for name in sorted(expected)]
    command = [
        "gh", "release", "create", tag, "-R", repo, "--target", commit,
        "--title", tag, "--notes", f"Source package {tag}.", *assets,
    ]  # fmt: skip
    if dry_run:
        line = "dry-run: " + " ".join(command)
        print(line)
        return line
    runner(command)
    return tag
