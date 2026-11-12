"""The Orphadata ``<JDBOR>`` header reader (ported from the atlas converter's tests)."""

from pathlib import Path

import pytest

from feeds.orphadata import jdbor_release
from feeds.package import PackageRefused

_JDBOR = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<JDBOR date="{date}" version="{version}" copyright="Orphanet (c) 2026">\n'
    "  <Availability/>\n</JDBOR>\n"
)


def test_the_version_and_date_are_read_and_a_wrong_root_or_blank_version_refused(
    tmp_path: Path,
) -> None:
    product = tmp_path / "en_product4.xml"
    product.write_text(
        _JDBOR.format(date="2027-01-02 03:04:05", version="1.3.50 / 4.1.9 [2026-12-01]")
    )
    assert jdbor_release(product) == ("1.3.50 / 4.1.9 [2026-12-01]", "2027-01-02 03:04:05")

    product.write_text('<?xml version="1.0"?>\n<JDBOR date="2027-01-02 03:04:05"/>\n')
    with pytest.raises(PackageRefused, match="no version"):
        jdbor_release(product)
    product.write_text(_JDBOR.format(date="  ", version="1.3.50"))
    with pytest.raises(PackageRefused, match="no version or no date"):
        jdbor_release(product)

    # Some other XML, whose attributes would otherwise be cited as Orphadata's release.
    product.write_text(
        _JDBOR.replace("JDBOR", "HPODisorderSetStatusList").format(
            date="2027-01-02 03:04:05", version="1.3.50"
        )
    )
    with pytest.raises(PackageRefused, match="not <JDBOR>"):
        jdbor_release(product)

    product.write_text("not xml")
    with pytest.raises(PackageRefused, match="not readable XML"):
        jdbor_release(product)
