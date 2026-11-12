"""The GenCC submissions export as the ``gencc`` package.

Moved from the atlas's ``scripts/convert_gencc.py``, whose reasoning is kept.

GenCC is a pass-through aggregator, not a primary curator. Its own terms: "The GenCC
does not independently verify the submitted information." It harmonises ~19 member
databases' vocabularies (gene, disease, mode of inheritance, classification) onto one
scale and publishes every submitter's verdict side by side, with none adjudicated. So
the table's key is (gene, disease, moi, submitter, sgc_id) -- the submitter is part of
the key on purpose: two submitters disagreeing about one gene-disease pair is the
normal case here, not an error to resolve.

``sgc_id`` is in the key too, and every row GenCC publishes is written: a submitter
can also disagree with *itself*. Measured by the atlas converter on 2026-08-03 against
the real export (30,410 rows, ``?format=new``), (gene, disease, moi, submitter) alone
still collides on 133 groups, and 70 of those 133 carry more than one distinct
``classification_title`` -- e.g. Ambry Genetics submitting both ``Limited``
(SGC-102815) and ``Strong`` (SGC-104042) for HGNC:20226 / MONDO:0859332 / Autosomal
recessive. ``version_number`` is identical within every one of the 133 groups, so
GenCC marks neither row as superseding the other, and nothing else in the export says
which is current. There is therefore no principled way to keep one and drop the other
without inventing a verdict GenCC itself does not assert, so every row is kept and
``sgc_id`` -- GenCC's own row id -- makes the key unique instead.

``?format=new`` is mandatory when downloading: GenCC's download page states the legacy
TSV format (keyed on a UUID rather than an ``sgc_id``/``version_number`` pair) "will be
removed and unavailable after September 30, 2026."
(``https://search.thegencc.org/download/action/submissions-export-tsv?format=new``;
the quota is 20 requests per IP per day.)

The carried columns are GenCC's harmonised ones, not the raw ``submitted_as_*`` ones
(which keep the submitter's own spelling), except ``submitted_as_date`` and
``submitted_as_public_report_url``, which have no harmonised twin.

THE EXPORT IS QUOTED TSV
------------------------
Measured on the 2026-11-11 export: 46,968 physical lines but 30,376 records, because
free-text cells the table does not carry contain quoted newlines. Splitting lines would
invent thousands of rows, so the file is parsed as quoted TSV (polars ``read_csv``,
``infer_schema_length=0``, as the converter does). The carried columns hold no ``"``
once parsed; a carried cell that does, or that holds a tab or line break, refuses
naming its ``sgc_id``, because the table is written unquoted and the atlas reads it
with polars' default quoting.

A row whose ``gene_curie`` is not ``HGNC:\\d+`` or whose ``disease_curie`` is not
``MONDO:\\d{7}`` is not a row of this table: it goes to ``dropped`` with its count.
Measured 2026-11-11: none does, so ``dropped`` is ``[]`` on real data; the guard stays
because GenCC's harmonisation is observed, not documented. ``ids`` is ``hgnc_id``.

TWO ENTRY POINTS
----------------
``build`` parses a live download; version is the input's sha256[:8] (GenCC states no
release) and ``source.version`` is ``retrieved <retrieved_on>``, the atlas's wording.
``bootstrap`` (owner decision F21) packages the atlas's committed mirror
``mirrors/gencc_submissions.tsv`` at one named atlas commit, once, because the
2026-08-03 export is no longer retrievable. It refuses any other bytes. The manifest's
input URL is the raw permalink at the atlas commit.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Sequence
from datetime import date
from pathlib import Path
from typing import Any, Final

import polars as pl

from feeds.package import PackageRefused, TableOut, write_package
from feeds.provenance import CommitOf, current_commit, stamp
from feeds.source import load_source
from feeds.version import package_version

HERE = Path(__file__).resolve().parent
SOURCE: Final = "gencc"
DOWNLOAD_URL: Final = (
    "https://search.thegencc.org/download/action/submissions-export-tsv?format=new"
)

# Upstream header -> table column, in the table's column order.
_COLUMNS: Final[dict[str, str]] = {
    "sgc_id": "sgc_id",
    "gene_curie": "gene",
    "gene_symbol": "gene_symbol",
    "disease_curie": "disease",
    "disease_title": "disease_label",
    "moi_title": "moi",
    "classification_title": "classification",
    "submitter_title": "submitter",
    "submitted_as_date": "submitted_on",
    "submitted_as_public_report_url": "report_url",
}
ORDER: Final = tuple(_COLUMNS.values())
_SORT_AT = tuple(ORDER.index(c) for c in ("gene", "disease", "moi", "submitter", "sgc_id"))
_HGNC = re.compile(r"HGNC:\d+")
_MONDO = re.compile(r"MONDO:\d{7}")
_UNWRITABLE = re.compile(r'["\t\r\n]')
_EXAMPLES: Final = 5

# The atlas's committed mirror at the bootstrap commit (atlas 4b1f8f6).
MIRROR_NAME: Final = "gencc_submissions.tsv"
MIRROR_SHA256: Final = "88c77272f0fd7bf807387ca7b5ae549d15183bd57e45a97f26cafc4d4e5c7808"
MIRROR_RETRIEVED_ON: Final = date(2026, 8, 3)
# The atlas commit whose mirror this is; the manifest URL is a permanent permalink to it.
ATLAS_COMMIT: Final = "4b1f8f68cbc30beaf932fbbc58e3617ebe401b6d"
_ATLAS_RAW = "https://raw.githubusercontent.com/chd-atlas/chd-atlas.github.io"


def _check_cells(rows: Sequence[Sequence[str]], origin: str) -> None:
    for row in rows:
        for name, cell in zip(ORDER, row, strict=True):
            if _UNWRITABLE.search(cell):
                raise PackageRefused(
                    f"{origin}: {row[0]} column {name} holds a quote, tab or line break "
                    f"({cell[:40]!r}); the table is written unquoted"
                )


def read_rows(tsv: Path) -> tuple[list[tuple[str, ...]], list[dict[str, Any]]]:
    """``(rows, dropped)``: the rows in table order, sorted by the key, and the drops.

    Pure: no package is written. Empty cells are ``""``.
    """
    rows, dropped, _ = _read(tsv)
    return rows, dropped


def _read(tsv: Path) -> tuple[list[tuple[str, ...]], list[dict[str, Any]], set[str]]:
    """``read_rows`` plus the distinct ``gene_curie`` values of every record read."""
    frame = pl.read_csv(tsv, separator="\t", infer_schema_length=0)
    missing = sorted(set(_COLUMNS) - set(frame.columns))
    if missing:
        raise PackageRefused(f"{tsv.name}: expected columns absent: {missing}")
    records = frame.select(list(_COLUMNS)).fill_null("").rows()
    kept: list[tuple[str, ...]] = []
    outside: list[str] = []
    genes_read = {record[1] for record in records if record[1]}
    for record in records:
        row = tuple(record)
        if _HGNC.fullmatch(row[1]) and _MONDO.fullmatch(row[3]):
            kept.append(row)
        else:
            outside.append(row[0])
    _check_cells(kept, tsv.name)
    kept.sort(key=lambda row: tuple(row[i] for i in _SORT_AT))
    dropped: list[dict[str, Any]] = []
    if outside:
        dropped.append(
            {
                "reason": "gene_curie is not HGNC:<digits> or disease_curie is not "
                "MONDO:<7 digits>",
                "count": len(outside),
                "examples": sorted(outside)[:_EXAMPLES],
            }
        )
    return kept, dropped, genes_read


def _write(
    rows: list[tuple[str, ...]],
    *,
    dropped: list[dict[str, Any]],
    genes_read: set[str],
    inputs: list[dict[str, Any]],
    retrieved_on: date,
    revision: int,
    sha256: str,
    built_by: dict[str, Any],
    out: Path,
) -> Path:
    source = load_source(HERE / "source.yaml")
    cells: list[Sequence[str | None]] = list(rows)
    return write_package(
        out,
        package=SOURCE,
        version=package_version(upstream_release=None, input_sha256=sha256, revision=revision),
        ids="hgnc_id",
        coverage="every_gene",
        tables=[TableOut("gencc_submissions", 1, ORDER, cells, "gene")],
        inputs=inputs,
        # Decision 4: every data record read, kept or dropped. ``genes`` is the distinct
        # non-empty gene_curie over those same records, so a dropped record's gene cell
        # (e.g. an Ensembl id) is counted; rows_in == table rows + sum(dropped counts).
        input_counts={
            "rows": len(rows) + sum(entry["count"] for entry in dropped),
            "genes": len(genes_read),
        },
        dropped=dropped,
        built_by=built_by,
        source={
            "name": source["name"],
            "version": f"retrieved {retrieved_on.isoformat()}",
            "retrieved_on": retrieved_on.isoformat(),
            "url": source["url"],
            "licence": source["licence"],
            "redistribution": source["redistribution"],
        },
        licence_evidence=source["licence_evidence"],
    )


def build(
    tsv: Path,
    *,
    retrieved_on: date,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Write the ``gencc`` package from a live download into the empty ``out``.

    ``revision`` has no default (F13): 1 for the first build of the input, 2, 3, ...
    for a recipe change or a rebuild whose ``manifest.json`` bytes differ.
    """
    built_by, _ = stamp(built_by, commit_of)
    rows, dropped, genes_read = _read(tsv)
    sha256 = hashlib.sha256(tsv.read_bytes()).hexdigest()
    return _write(
        rows,
        dropped=dropped,
        genes_read=genes_read,
        inputs=[
            {
                "name": tsv.name,
                "url": DOWNLOAD_URL,
                "sha256": sha256,
                "retrieved_on": retrieved_on.isoformat(),
                "upstream_release": None,
            }
        ],
        retrieved_on=retrieved_on,
        revision=revision,
        sha256=sha256,
        built_by=built_by,
        out=out,
    )


def bootstrap(
    mirror: Path,
    *,
    atlas_commit: str,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Package the atlas's committed mirror at ``atlas_commit`` (decision F21), once.

    Refuses any ``atlas_commit`` but ``ATLAS_COMMIT``, any file whose sha256 is not
    ``MIRROR_SHA256`` and any header that is not the table's column order.
    """
    built_by, _ = stamp(built_by, commit_of)
    if atlas_commit != ATLAS_COMMIT:
        raise PackageRefused(
            f"atlas_commit {atlas_commit!r} is not the pinned {ATLAS_COMMIT}; the mirror's "
            "permalink must name the commit it was read at"
        )
    data = mirror.read_bytes()
    sha256 = hashlib.sha256(data).hexdigest()
    if sha256 != MIRROR_SHA256:
        raise PackageRefused(
            f"{mirror.name} has sha256 {sha256}, not the pinned {MIRROR_SHA256}; "
            "the bootstrap packages the atlas's mirror at one commit only"
        )
    lines = data.decode("utf-8").splitlines()
    if tuple(lines[0].split("\t")) != ORDER:
        raise PackageRefused(
            f"{mirror.name}: header {lines[0]!r} is not the table's column order {list(ORDER)}"
        )
    rows: list[tuple[str, ...]] = []
    for number, line in enumerate(lines[1:], start=2):
        cells = tuple(line.split("\t"))
        if len(cells) != len(ORDER):
            raise PackageRefused(f"{mirror.name} line {number}: {len(cells)} cells")
        rows.append(cells)
    _check_cells(rows, mirror.name)
    return _write(
        rows,
        dropped=[],
        genes_read={row[1] for row in rows},
        inputs=[
            {
                "name": MIRROR_NAME,
                "url": f"{_ATLAS_RAW}/{atlas_commit}/mirrors/{MIRROR_NAME}",
                "sha256": sha256,
                "retrieved_on": MIRROR_RETRIEVED_ON.isoformat(),
                "upstream_release": None,
            }
        ],
        retrieved_on=MIRROR_RETRIEVED_ON,
        revision=revision,
        sha256=sha256,
        built_by=built_by,
        out=out,
    )
