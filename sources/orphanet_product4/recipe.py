"""Orphadata product 4 (``en_product4.xml``) as the ``orphanet_product4`` package.

Moved from the atlas's ``scripts/convert_orphanet.py`` (``_phenotype_rows``), whose
reasoning is kept, with one change of kind: the converter kept only HPO ids under
HP:0001626 (cardiovascular); the package carries EVERY ``HPODisorderAssociation``, and
the atlas filters at ingest (owner decision F23). There is no HPO file and no HPO
filter here.

LICENCE AND DATA VERSION
------------------------
The file declares ``CC-BY-4.0`` in its own ``<Licence>`` block, and Orphadata's legal
notice says the same of its Orphadata Science files, asking to be cited as "Orphadata
Science: Free access data from Orphanet. (c) INSERM 1999. Available on
http://sciences.orphadata.com/. Data version [date]". That data version is the root
``<JDBOR date=... version=...>`` header, read by ``feeds.orphadata.jdbor_release`` and
written into the manifest's ``source.file_version`` / ``file_date`` beside the sha256,
so the checksum, the version and the date cannot come from three different downloads.
The band (``HPOFrequency/Name``) is Orphanet's claim and is carried verbatim, never
mapped to a number.

EXACT DUPLICATES
----------------
Product 4 annotates some (orpha, hpo) pairs twice. Where every table column repeats the
row is collapsed into one ``dropped`` entry (count, up to five
``ORPHA:<code> HP:<id>`` examples); where only the band differs both rows are kept,
which is why the band is in the table key ``(orpha, hpo, frequency)``. "Identical" means
in the carried columns only: for 1 of the 24 (ORPHA:2806 / HP:0010856, association ids
648555 and 648556) the file's ``DiagnosticCriteria``, which the table does not carry,
differs. Measured by
this recipe on the pinned file (see the README): 116,664 annotations, 116,640 rows, 24
collapsed. ``input_counts.rows`` is the annotations read (plan Decision 4).

The build accepts only the pinned download: the sha256 below is the one the atlas's
``mirrors/sources.yaml`` records for this source, and ``retrieved_on`` is the atlas's
2026-10-25 (plan Decision 3). ``ids`` is ``hgnc_id``: the package has no gene column.
"""

from __future__ import annotations

import hashlib
import xml.etree.ElementTree as ET
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final

from feeds.orphadata import jdbor_release
from feeds.package import PackageRefused, TableOut, write_package
from feeds.provenance import CommitOf, current_commit, stamp
from feeds.source import load_source
from feeds.version import package_version

HERE = Path(__file__).resolve().parent
SOURCE: Final = "orphanet_product4"
DOWNLOAD_URL: Final = "https://www.orphadata.com/data/xml/en_product4.xml"
PINNED_SHA256: str = "4f44e8a61201399911aa1ba44a293c0ccaa5ce11272c47d40862255cb72f6b32"
RETRIEVED_ON: Final = "2026-10-25"
# The pinned file's own <JDBOR> header (atlas mirrors/sources.yaml, read 2026-10-26).
FILE_VERSION: Final = "1.3.42 / 4.1.8 [2025-03-03]"
FILE_DATE: Final = "2026-06-23 07:57:18"

ORDER: Final = ("orpha", "disorder_name", "hpo", "hpo_label", "frequency", "source")
_SORT_AT: Final = (0, 2, 4)  # (orpha, hpo, frequency)
_EXAMPLES: Final = 5


def _text(element: ET.Element, path: str) -> str:
    return (element.findtext(path) or "").strip()


def _rows(element: ET.Element) -> list[tuple[str, ...]]:
    """One ``HPODisorderSetStatus``, as a row per annotation, none filtered."""
    disorder = element.find("Disorder")
    if disorder is None:
        return []
    orpha = f"ORPHA:{_text(disorder, 'OrphaCode')}"
    name = _text(disorder, "Name")
    return [
        (
            orpha,
            name,
            _text(association, "HPO/HPOId"),
            _text(association, "HPO/HPOTerm"),
            # Verbatim, never mapped to a number: the band is Orphanet's claim.
            _text(association, "HPOFrequency/Name"),
            SOURCE,
        )
        for association in disorder.iter("HPODisorderAssociation")
    ]


def _read(xml: Path) -> tuple[list[tuple[str, ...]], list[dict[str, Any]], int]:
    """``(rows, dropped, annotations read)``."""
    annotations: list[tuple[str, ...]] = []
    try:
        for _, element in ET.iterparse(xml):
            if element.tag == "HPODisorderSetStatus":
                annotations.extend(_rows(element))
                element.clear()
    except ET.ParseError as exc:
        raise PackageRefused(f"{xml.name}: not readable XML: {exc}") from exc
    counts = Counter(annotations)
    rows = sorted(counts, key=lambda row: tuple(row[i] for i in _SORT_AT))
    collapsed = len(annotations) - len(rows)
    dropped: list[dict[str, Any]] = []
    if collapsed:
        pairs = sorted({f"{row[0]} {row[2]}" for row, n in counts.items() if n > 1})
        dropped.append(
            {
                "reason": f"annotation identical to another in every table column in {xml.name}",
                "count": collapsed,
                "examples": pairs[:_EXAMPLES],
            }
        )
    return rows, dropped, len(annotations)


def read_rows(xml: Path) -> tuple[list[tuple[str, ...]], list[dict[str, Any]]]:
    """``(rows, dropped)``: the rows in table order sorted by the key, and the collapse.

    Pure: no package is written, and the sha256 pin is not checked.
    """
    rows, dropped, _ = _read(xml)
    return rows, dropped


def build(
    xml: Path,
    *,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Write the ``orphanet_product4`` package from the pinned download into the empty ``out``.

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
    rows, dropped, read = _read(xml)
    source = load_source(HERE / "source.yaml")
    cells: list[Sequence[str | None]] = list(rows)
    return write_package(
        out,
        package=SOURCE,
        version=package_version(upstream_release=None, input_sha256=sha256, revision=revision),
        ids="hgnc_id",
        coverage="every_gene",
        tables=[TableOut("orphanet_phenotypes", 1, ORDER, cells, None)],
        inputs=[
            {
                "name": xml.name,
                "url": DOWNLOAD_URL,
                "sha256": sha256,
                "retrieved_on": RETRIEVED_ON,
                "upstream_release": file_version,
            }
        ],
        input_counts={"rows": read},
        dropped=dropped,
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
