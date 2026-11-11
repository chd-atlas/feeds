"""Vendor the atlas's exported JSON Schemas at one atlas commit.

Usage: sync_atlas_schemas.py <atlas_checkout> [--commit SHA]

Reads ``schemas/tables/*.json`` and ``schemas/package_manifest.schema.json``
from the commit's tree (``git show``), never from the working tree, writes them
under ``schemas/atlas/``, deletes vendored files the commit no longer has, and
records the full commit sha in ``schemas/atlas/ATLAS_COMMIT``.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "schemas" / "atlas"


def _git(atlas: Path, *args: str) -> bytes:
    return subprocess.run(["git", "-C", str(atlas), *args], check=True, capture_output=True).stdout


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("atlas", type=Path)
    parser.add_argument("--commit", default=None)
    args = parser.parse_args(argv)
    sha = _git(args.atlas, "rev-parse", "--verify", f"{args.commit or 'HEAD'}^{{commit}}")
    sha_text = sha.decode().strip()
    tables = [
        line
        for line in _git(args.atlas, "ls-tree", "--name-only", sha_text, "schemas/tables/")
        .decode()
        .splitlines()
        if line.endswith(".json")
    ]
    wanted = {
        Path(line).relative_to("schemas"): line
        for line in [*tables, "schemas/package_manifest.schema.json"]
    }
    contents = {rel: _git(args.atlas, "show", f"{sha_text}:{path}") for rel, path in wanted.items()}
    (DEST / "tables").mkdir(parents=True, exist_ok=True)
    for stale in [*DEST.glob("tables/*.json"), *DEST.glob("*.json")]:
        if stale.relative_to(DEST) not in contents:
            stale.unlink()
    for rel, data in contents.items():
        (DEST / rel).write_bytes(data)
    (DEST / "ATLAS_COMMIT").write_text(sha_text + "\n")
    print(f"vendored {len(contents)} schemas at {sha_text}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
