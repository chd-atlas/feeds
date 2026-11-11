"""Read a recipe's ``source.yaml``: what the source is called, where it lives, its terms."""

from __future__ import annotations

from pathlib import Path

import yaml

REQUIRED = ("name", "url", "licence", "redistribution", "licence_evidence")


class SourceRefused(Exception):
    """``source.yaml`` is unreadable, is not a mapping, or lacks a required key."""


def load_source(path: Path) -> dict[str, str]:
    """The five keys of spec section 4.3, each a non-blank string.

    A key that is missing, not a string or blank after stripping is refused by
    name, so a package is never cut with an empty licence or evidence line.
    """
    try:
        loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise SourceRefused(f"{path}: cannot read: {exc}") from exc
    if not isinstance(loaded, dict):
        raise SourceRefused(f"{path}: expected a mapping of keys, got {type(loaded).__name__}")
    out: dict[str, str] = {}
    for key in REQUIRED:
        value = loaded.get(key)
        if not isinstance(value, str) or not value.strip():
            raise SourceRefused(f"{path}: '{key}' is missing or blank")
        out[key] = value.strip()
    return out
