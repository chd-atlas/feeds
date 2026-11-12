"""The Mondo SSSOM mapping file as the ``mondo_sssom`` package.

Moved from the atlas's ``scripts/convert_orphanet.py`` (``_xref_rows`` and the SSSOM
reading in ``main``), whose rules are kept: a row is a ``skos:exactMatch`` from a
``MONDO:`` subject to an ``Orphanet:`` or ``OMIM:`` object; ``Orphanet:<n>`` is written
``ORPHA:<n>`` as Orphanet spells it; rows are de-duplicated and sorted by
``(mondo, target)``. The ``#`` lines before the header are SSSOM metadata and are skipped.

One change: the converter kept only MONDO terms an atlas validity mirror names. The
package carries every MONDO subject, and the atlas filters at ingest (owner decision
F23), so this recipe reads no atlas data. Other predicates and other targets are not
rows and not drops: they are not this table's subject (plan Decision 4).

The build accepts only the pinned download (the sha256 the atlas's
``mirrors/sources.yaml`` records, Mondo release ``v2026-10-06``) with ``retrieved_on``
2026-10-25 (plan Decision 3). ``input_counts.rows`` is the carried matches read; a pair
repeated in the file is written once and the repeat is listed under ``dropped``.
``ids`` is ``hgnc_id``: the package has no gene column.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any, Final

from feeds.package import PackageRefused, TableOut, write_package
from feeds.provenance import CommitOf, current_commit, stamp
from feeds.source import load_source
from feeds.version import package_version

HERE = Path(__file__).resolve().parent
SOURCE: Final = "mondo_sssom"
PINNED_SHA256: str = "ec15afdb32fab3d3cb679d357e25ef4fb38b5e2c468ff4ac5456da10e2f6c80b"
RETRIEVED_ON: Final = "2026-10-25"
UPSTREAM_RELEASE: Final = "v2026-10-06"
VERSION: Final = "2026-10-06"
DOWNLOAD_URL: Final = (
    "https://raw.githubusercontent.com/monarch-initiative/mondo/v2026-10-06"
    "/src/ontology/mappings/mondo.sssom.tsv"
)

ORDER: Final = ("mondo", "target", "source")
_UNWRITABLE: Final = ('"', "\t", "\n", "\r")


def _body(path: Path) -> Iterator[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        for line in handle:
            if not line.startswith("#"):
                yield line


def read_rows(sssom: Path) -> tuple[list[tuple[str, ...]], int]:
    """The distinct carried rows sorted by ``(mondo, target)``, and the matches read.

    Pure: no package is written, and the sha256 pin is not checked.
    """
    # QUOTE_NONE, unlike the converter's default quoting: SSSOM cells are unquoted (the
    # pinned file holds 0 ``"``), and default quoting could merge records across a stray
    # quote where QUOTE_NONE cannot. A quote in a carried cell is refused below.
    reader = csv.DictReader(_body(sssom), delimiter="\t", quoting=csv.QUOTE_NONE)
    pairs: list[tuple[str, str]] = []
    for record in reader:
        subject, predicate, obj = (
            record.get("subject_id"),
            record.get("predicate_id"),
            record.get("object_id"),
        )
        if subject is None or predicate is None or obj is None:
            raise PackageRefused(f"{sssom.name}: no subject_id, predicate_id or object_id column")
        if predicate != "skos:exactMatch" or not subject.startswith("MONDO:"):
            continue
        if obj.startswith("Orphanet:"):
            target = "ORPHA:" + obj.split(":", 1)[1]
        elif obj.startswith("OMIM:"):
            target = obj
        else:
            continue
        if any(ch in cell for cell in (subject, target) for ch in _UNWRITABLE):
            raise PackageRefused(
                f"{sssom.name}: {subject} -> {target!r} holds a quote, tab or newline"
            )
        pairs.append((subject, target))
    rows: list[tuple[str, ...]] = [(m, t, SOURCE) for m, t in sorted(set(pairs))]
    return rows, len(pairs)


def build(
    sssom: Path,
    *,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Write the ``mondo_sssom`` package from the pinned download into the empty ``out``.

    ``revision`` has no default (F13): 1 for the first build of the input, 2, 3, ...
    for a recipe change or a rebuild whose ``manifest.json`` bytes differ.
    """
    built_by, _ = stamp(built_by, commit_of)
    sha256 = hashlib.sha256(sssom.read_bytes()).hexdigest()
    if sha256 != PINNED_SHA256:
        raise PackageRefused(
            f"{sssom.name} has sha256 {sha256}, not the pinned {PINNED_SHA256}; "
            "this recipe packages the 2026-10-25 download only"
        )
    rows, matches = read_rows(sssom)
    source = load_source(HERE / "source.yaml")
    cells: list[Sequence[str | None]] = list(rows)
    dropped: list[dict[str, Any]] = []
    if matches != len(rows):
        dropped.append({"reason": "pair repeated in the file", "count": matches - len(rows)})
    return write_package(
        out,
        package=SOURCE,
        version=package_version(upstream_release=VERSION, input_sha256=sha256, revision=revision),
        ids="hgnc_id",
        coverage="every_gene",
        tables=[TableOut("mondo_xrefs", 1, ORDER, cells, None)],
        inputs=[
            {
                "name": sssom.name,
                "url": DOWNLOAD_URL,
                "sha256": sha256,
                "retrieved_on": RETRIEVED_ON,
                "upstream_release": UPSTREAM_RELEASE,
            }
        ],
        input_counts={"rows": matches},
        dropped=dropped,
        built_by=built_by,
        source={
            "name": source["name"],
            "version": f"sha256:{sha256}",
            "retrieved_on": RETRIEVED_ON,
            "url": source["url"],
            "licence": source["licence"],
            "redistribution": source["redistribution"],
        },
        licence_evidence=source["licence_evidence"],
    )
