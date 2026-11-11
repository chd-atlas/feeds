from pathlib import Path

import pytest

from feeds.source import REQUIRED, SourceRefused, load_source

_GOOD = {
    "name": "Demo",
    "url": "https://example.org/",
    "licence": "CC0-1.0",
    "redistribution": "permitted",
    "licence_evidence": "  all data is released under CC0  ",
}


def _write(path: Path, entries: dict[str, str | None]) -> Path:
    path.write_text(
        "".join(f"{k}: {v!r}\n" if v is not None else f"{k}:\n" for k, v in entries.items()),
        encoding="utf-8",
    )
    return path


def test_the_five_keys_are_returned_stripped(tmp_path: Path) -> None:
    got = load_source(_write(tmp_path / "source.yaml", _GOOD))
    assert got == _GOOD | {"licence_evidence": "all data is released under CC0"}


@pytest.mark.parametrize("key", REQUIRED)
@pytest.mark.parametrize("bad", [None, "", "   "], ids=["missing-value", "empty", "blank"])
def test_a_missing_or_blank_key_is_refused_by_name(
    tmp_path: Path, key: str, bad: str | None
) -> None:
    path = _write(tmp_path / "source.yaml", _GOOD | {key: bad})
    with pytest.raises(SourceRefused, match=f"'{key}'"):
        load_source(path)
    absent = {k: v for k, v in _GOOD.items() if k != key}
    with pytest.raises(SourceRefused, match=f"'{key}'"):
        load_source(_write(tmp_path / "absent.yaml", absent))


def test_a_file_that_is_not_a_mapping_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "source.yaml"
    path.write_text("- a\n- b\n", encoding="utf-8")
    with pytest.raises(SourceRefused, match="mapping"):
        load_source(path)
