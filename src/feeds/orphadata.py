"""The header every Orphadata XML file states on its root ``<JDBOR>`` element."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from feeds.package import PackageRefused


def jdbor_release(path: Path) -> tuple[str, str]:
    """The ``version`` and ``date`` an Orphadata file states on its root ``<JDBOR>`` element.

    Moved from the atlas's ``scripts/convert_orphanet.py``. Read from the first start
    event only, so a 50 MB file costs one chunk. Refuses a file whose root is not
    ``JDBOR`` or that states either one blank: a release recorded as nothing would
    cite as nothing.
    """
    try:
        for _, element in ET.iterparse(path, events=("start",)):
            if element.tag != "JDBOR":
                raise PackageRefused(f"{path.name}: root element is <{element.tag}>, not <JDBOR>")
            version = element.get("version", "").strip()
            date = element.get("date", "").strip()
            if not version or not date:
                raise PackageRefused(f"{path.name}: <JDBOR> states no version or no date")
            return version, date
    except ET.ParseError as exc:
        raise PackageRefused(f"{path.name}: not readable XML: {exc}") from exc
    raise PackageRefused(f"{path.name}: no XML element")
