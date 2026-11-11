"""Provenance the build records for itself, so no operator types it.

``recipe_commit`` is read from the checkout that runs the build, and a build
from a checkout with uncommitted recipe code is refused: the commit a manifest
names must be the code that produced its bytes. ``hvantk_builder_commit`` reads
hvantk's own sidecar for the commit of the hvantk code that ran.
"""

from __future__ import annotations

import json
import re
import subprocess
from collections.abc import Callable
from datetime import date, datetime
from pathlib import Path
from typing import Any

from feeds.package import PackageRefused

REPO_ROOT = Path(__file__).resolve().parents[2]
# Everything whose bytes can change a package. README and tests are not listed.
BUILD_PATHS = ("src", "sources", "schemas", "pyproject.toml", "uv.lock")
COMMIT = re.compile(r"[0-9a-f]{40}")

CommitOf = Callable[[], str]


def _git(root: Path, *args: str) -> str:
    try:
        done = subprocess.run(
            ["git", "-C", str(root), *args], check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, "stderr", "") or str(exc)
        raise PackageRefused(
            f"cannot read the recipe commit from {root}: {detail.strip()}"
        ) from exc
    return done.stdout


def recipe_commit(root: Path) -> str:
    """``git rev-parse HEAD`` of ``root``, refused when the build paths are not clean.

    Untracked files under the build paths count as dirty; ignored files do not.
    """
    dirty = _git(root, "status", "--porcelain", "--", *BUILD_PATHS).strip()
    if dirty:
        shown = "; ".join(line.strip() for line in dirty.splitlines()[:5])
        raise PackageRefused(
            f"{root} has uncommitted changes under {', '.join(BUILD_PATHS)} ({shown}); "
            "build from a clean checkout of the commit you will publish"
        )
    commit = _git(root, "rev-parse", "HEAD").strip()
    if not COMMIT.fullmatch(commit):
        raise PackageRefused(f"git rev-parse HEAD gave {commit!r}, not a 40-hex commit")
    return commit


def current_commit() -> str:
    """The commit of this feeds checkout: the default ``commit_of`` of every recipe."""
    return recipe_commit(REPO_ROOT)


def stamp(built_by: dict[str, Any], commit_of: CommitOf) -> tuple[dict[str, Any], str]:
    """``built_by`` with ``recipe_commit`` set from the build, and that commit.

    A caller-supplied ``recipe_commit`` is accepted only when it equals the build's.
    """
    commit = commit_of()
    if not COMMIT.fullmatch(commit):
        raise PackageRefused(f"the recipe commit {commit!r} is not a 40-hex sha")
    supplied = built_by.get("recipe_commit")
    if supplied is not None and supplied != commit:
        raise PackageRefused(
            f'built_by["recipe_commit"] is {supplied!r} but this build runs at {commit}; '
            "omit it and the build records its own"
        )
    return {**built_by, "recipe_commit": commit}, commit


def require_hvantk_release(built_by: dict[str, Any]) -> str:
    """The ``hvantk_release`` of ``built_by``, which a recipe built on hvantk must record."""
    release = built_by.get("hvantk_release")
    if not isinstance(release, str) or not release.strip():
        raise PackageRefused('built_by["hvantk_release"] is required (the hvantk commit that ran)')
    return release


def hvantk_sidecar(table: Path) -> tuple[str, date]:
    """``(builder_commit, build date)`` from hvantk's ``<table>.provenance.json``.

    The commit must be a full 40-hex sha: hvantk records ``hvantk <version>`` when it
    finds neither a git checkout nor ``HVANTK_BUILDER_COMMIT``, and that run is refused.
    """
    path = table.with_name(table.name + ".provenance.json")
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise PackageRefused(
            f"hvantk's provenance sidecar {path.name} is missing or unreadable next to "
            f"{table.name}: {exc}"
        ) from exc
    commit = loaded.get("builder_commit") if isinstance(loaded, dict) else None
    if not isinstance(commit, str) or not COMMIT.fullmatch(commit):
        raise PackageRefused(
            f"{path.name}: builder_commit is {commit!r}, not a 40-hex commit "
            "(run the job with HVANTK_BUILDER_COMMIT set to the hvantk checkout's commit)"
        )
    stamp_text = loaded.get("build_timestamp")
    try:
        made = datetime.fromisoformat(str(stamp_text)).date()
    except ValueError as exc:
        raise PackageRefused(f"{path.name}: build_timestamp {stamp_text!r} is not ISO") from exc
    return commit, made
