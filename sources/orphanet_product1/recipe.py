"""Orphadata product 1 (``en_product1.xml``) as the ``orphanet_product1`` package.

Moved from the atlas's ``scripts/convert_orphanet.py`` (``_omim_rows``), whose reasoning
is kept: every Orphanet-to-OMIM reference, whatever its relation or validation status.
The atlas's rule keeps only validated ``E``/``BTNT``; the package keeps the rest so the
refusals stay visible. References of any other source (ICD-10, ...) are not rows and not
drops: they are not this table's subject.

The relation is the code at the front of Orphanet's name ("BTNT (ORPHAcode is broader
than ...)" -> "BTNT"). Licence, citation and data version are as for product 4 (see
``sources/orphanet_product4/recipe.py``): the file's ``<JDBOR>`` header is read by
``feeds.orphadata.jdbor_release`` and written into the manifest beside the sha256.

The build accepts only the pinned download (the sha256 the atlas's
``mirrors/sources.yaml`` records) with ``retrieved_on`` 2026-10-25 (plan Decision 3).
``input_counts.rows`` is the OMIM references read (plan Decision 4). ``ids`` is
``hgnc_id``: the package has no gene column.
"""

from __future__ import annotations

import hashlib
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final

from feeds.orphadata import jdbor_release
from feeds.package import PackageRefused, TableOut, write_package
from feeds.provenance import CommitOf, current_commit, stamp
from feeds.source import load_source
from feeds.version import package_version

HERE = Path(__file__).resolve().parent
SOURCE: Final = "orphanet_product1"
DOWNLOAD_URL: Final = "https://www.orphadata.com/data/xml/en_product1.xml"
PINNED_SHA256: str = "df8d562a0c6011af36a74eb4000ce81ca7d723e8031010819fb71727c0962bbb"
RETRIEVED_ON: Final = "2026-10-25"
# The pinned file's own <JDBOR> header (atlas mirrors/sources.yaml, read 2026-10-26).
FILE_VERSION: Final = "1.3.42 / 4.1.8 [2025-03-03]"
FILE_DATE: Final = "2026-06-23 07:53:50"

ORDER: Final = ("orpha", "omim", "relation", "validation", "source")


def _text(element: ET.Element, path: str) -> str:
    return (element.findtext(path) or "").strip()


def _rows(disorder: ET.Element) -> list[tuple[str, ...]]:
    """One product 1 ``Disorder``, as a row per OMIM reference."""
    orpha = f"ORPHA:{_text(disorder, 'OrphaCode')}"
    return [
        (
            orpha,
            f"OMIM:{_text(reference, 'Reference')}",
            # "BTNT (ORPHAcode is broader than ...)" -> "BTNT"
            _text(reference, "DisorderMappingRelation/Name").split(" ", 1)[0],
            _text(reference, "DisorderMappingValidationStatus/Name"),
            SOURCE,
        )
        for reference in disorder.iter("ExternalReference")
        if _text(reference, "Source") == "OMIM"
    ]


def read_rows(xml: Path) -> list[tuple[str, ...]]:
    """The OMIM reference rows in table order, sorted by ``(orpha, omim)``.

    Pure: no package is written, and the sha256 pin is not checked.
    """
    rows: list[tuple[str, ...]] = []
    try:
        for _, element in ET.iterparse(xml):
            if element.tag == "Disorder":
                rows.extend(_rows(element))
                element.clear()
    except ET.ParseError as exc:
        raise PackageRefused(f"{xml.name}: not readable XML: {exc}") from exc
    rows.sort(key=lambda row: (row[0], row[1]))
    return rows


def build(
    xml: Path,
    *,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Write the ``orphanet_product1`` package from the pinned download into the empty ``out``.

    ``revision`` has no default (F13): 1 for the first build of the input, 2, 3, ...
    for a recipe change or a rebuild whose ``manifest.json`` bytes differ.
    """
    built_by, _ = stamp(built_by, commit_of)
    sha256 = hashlib.sha256(xml.read_bytes()).hexdigest()
    if sha256 != PINNED_SHA256:
        raise PackageRefused(
            f"{xml.name} has sha256 {sha256}, not the pinned {PINNED_SHA256}; "
            "this recipe packages the 2026-10-25 download only"
        )
    file_version, file_date = jdbor_release(xml)
    rows = read_rows(xml)
    source = load_source(HERE / "source.yaml")
    cells: list[Sequence[str | None]] = list(rows)
    return write_package(
        out,
        package=SOURCE,
        version=package_version(upstream_release=None, input_sha256=sha256, revision=revision),
        ids="hgnc_id",
        coverage="every_gene",
        tables=[TableOut("orphanet_omim", 1, ORDER, cells, None)],
        inputs=[
            {
                "name": xml.name,
                "url": DOWNLOAD_URL,
                "sha256": sha256,
                "retrieved_on": RETRIEVED_ON,
                "upstream_release": file_version,
            }
        ],
        input_counts={"rows": len(rows)},
        dropped=[],
        built_by=built_by,
        source={
            "name": source["name"],
            "version": f"sha256:{sha256}",
            "file_version": file_version,
            "file_date": file_date,
            "retrieved_on": RETRIEVED_ON,
            "url": source["url"],
            "licence": source["licence"],
            "redistribution": source["redistribution"],
        },
        licence_evidence=source["licence_evidence"],
    )
