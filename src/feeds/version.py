"""The version of a package: the upstream release, else the input's hash prefix."""

from __future__ import annotations

import re

# The atlas's models/package.py SEGMENT: a version is a directory name.
SEGMENT = r"^(?:[A-Za-z0-9_-][A-Za-z0-9._-]*|\.[A-Za-z0-9_-][A-Za-z0-9._-]*|\.\.[A-Za-z0-9._-]+)$"
_SHA256 = re.compile(r"[0-9a-f]{64}")


def package_version(*, upstream_release: str | None, input_sha256: str, revision: int = 1) -> str:
    """``<upstream release>`` or ``<first 8 hex of input_sha256>``, plus ``-r<n>`` for n >= 2."""
    if not _SHA256.fullmatch(input_sha256):
        raise ValueError("input_sha256 must be 64 lowercase hex characters")
    if revision < 1:
        raise ValueError(f"revision must be at least 1, got {revision}")
    base = upstream_release if upstream_release else input_sha256[:8]
    version = base if revision == 1 else f"{base}-r{revision}"
    if not re.fullmatch(SEGMENT, version):
        raise ValueError(f"version '{version}' is not a valid path segment")
    return version
