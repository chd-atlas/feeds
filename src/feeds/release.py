"""Publish a built package as the GitHub release ``<package>/<version>``."""

from __future__ import annotations

import hashlib
import json
import re
import shlex
import subprocess
import tempfile
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

Runner = Callable[[Sequence[str]], str]
REPO_ROOT = Path(__file__).resolve().parents[2]


class ReleaseRefused(Exception):
    """The release would be unsafe or ambiguous; nothing was published."""


def _run(command: Sequence[str]) -> str:
    done = subprocess.run(command, check=True, capture_output=True, text=True)
    return done.stdout


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def release_notes(manifest: dict[str, Any], tag: str) -> str:
    """The release page's text: what the package is and whose terms its data carry."""
    source = manifest["source"]
    return (
        f"Source package {tag}: {source['name']}.\n\n"
        f"Licence of the data: {source['licence']}. Redistribution: "
        f"{source['redistribution']}.\n"
        "These are the source's terms, not this repository's: the repository's Apache-2.0 "
        "LICENSE covers the code only. The evidence for the terms is `licence_evidence` in "
        "`manifest.json`.\n"
    )


def _names(url: str, repo: str) -> bool:
    return re.search(rf"(?:^|[:/]){re.escape(repo)}(?:\.git)?/?$", url.strip()) is not None


def publish(
    dist: Path,
    *,
    repo: str = "chd-atlas/feeds",
    dry_run: bool = True,
    runner: Runner = _run,
    root: Path = REPO_ROOT,
) -> str:
    """Return ``"unchanged"``, ``"dry-run: <command>"`` or the tag that was created.

    Contacts the network even when ``dry_run``: ``git fetch origin dev`` in ``root``
    (the feeds checkout, not the working directory), ``gh release view`` and, when the
    tag exists, ``gh release download``. ``dry_run`` only withholds ``gh release create``.
    """
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
    for entry in manifest["tables"]:
        if _sha(dist / entry["file"]) != entry["sha256"]:
            raise ReleaseRefused(f"{entry['file']} does not match the manifest's sha256")
    git = ["git", "-C", str(root)]
    origin = runner([*git, "remote", "get-url", "origin"]).strip()
    if not _names(origin, repo):
        raise ReleaseRefused(f"origin of {root} is {origin!r}, not {repo}; refusing to publish")
    runner([*git, "fetch", "origin", "dev"])
    try:
        listing = runner([*git, "branch", "-r", "--contains", commit])
    except subprocess.CalledProcessError as exc:
        if "no such commit" not in f"{exc.stderr or ''}{exc.output or ''}":
            raise
        raise ReleaseRefused(
            f"recipe commit {commit} is unknown to {root}'s origin; push it first"
        ) from exc
    branches = [b.strip() for b in listing.splitlines()]
    if "origin/dev" not in branches:
        raise ReleaseRefused(f"recipe commit {commit} is not on origin/dev; push it first")
    try:
        runner(["gh", "release", "view", tag, "-R", repo])
        exists = True
    except subprocess.CalledProcessError as exc:
        if "release not found" not in f"{exc.stderr or ''}{exc.output or ''}":
            raise
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
        "--title", tag, "--notes", release_notes(manifest, tag), *assets,
    ]  # fmt: skip
    if dry_run:
        line = "dry-run: " + shlex.join(command)
        print(line)
        return line
    runner(command)
    return tag
