"""The ClinGen gene-disease validity bulk CSV as the ``clingen`` package.

Moved from the atlas's ``scripts/convert_clingen.py``, whose reasoning is kept.
The source file carries no version tag -- only a same-day ``FILE CREATED`` banner
stamp -- so the banner date is the package version and the banner text is the
``source.version`` the atlas pins. Fetching at build time instead would drift the
corpus between builds with nothing to diff.

THE PREAMBLE IS A TRAP, AND THE TRAP IS SILENT
----------------------------------------------
Lines 1-3 are a title block, line 4 and line 6 are ``+++`` separator rows, and
line 5 is the real header. A reader with default arguments takes line 1 as the
header and produces ten columns named after the title block -- no exception, no
empty table, just every column name wrong. So the header is taken as the line
after the four-line banner, the banner's second line must say
``FILE CREATED: YYYY-MM-DD``, an expected column that is absent refuses naming it,
and a row whose ``GENE SYMBOL`` starts with ``+++`` is layout: neither counted
nor dropped. ``ids`` is ``hgnc_id`` (the CSV carries the HGNC id itself) and the
package covers every gene ClinGen has curated; the atlas does not registry-filter
this table on ingest because its gene registry is derived from it.

TWO ENTRY POINTS
----------------
``build`` parses a live download. ``bootstrap`` (owner decision F21) packages the
atlas's committed mirror ``mirrors/clingen_gene_validity.tsv`` at one named atlas
commit, once, because the 2026-08-03 download is no longer retrievable: ClinGen
serves only the current file. It refuses any other bytes. The two are kept
apart so that the live parser is never asked to read a mirror and a mirror is never
passed off as a download; the manifest's input URL is the raw permalink at the atlas
commit, and its ``upstream_release`` is the banner text of the file the mirror was
made from.

Measured 2026-11-11, by running ``build`` on that day's download: it stops in
``write_package`` with exactly one primary-key error, a repeated (gene, disease, moi)
key for FOXE3 ('HGNC:3808', 'MONDO:0019625', 'AD'). ``read_rows`` is pure, so row
parity with the atlas converter can be measured without packaging.
"""

from __future__ import annotations

import csv
import hashlib
import re
from collections.abc import Sequence
from datetime import date
from pathlib import Path
from typing import Any, Final

from feeds.package import PackageRefused, TableOut, write_package
from feeds.provenance import CommitOf, current_commit, stamp
from feeds.source import load_source
from feeds.version import package_version

HERE = Path(__file__).resolve().parent
SOURCE: Final = "clingen"
DOWNLOAD_URL: Final = "https://search.clinicalgenome.org/kb/gene-validity/download"

# Lines 1-3 are a title block; line 4 is a separator; line 5 is the header.
_PREAMBLE_ROWS = 4
_SEPARATOR = "+++"
_BANNER = re.compile(r"FILE CREATED: (\d{4}-\d{2}-\d{2})")

# Upstream header -> table column, in the table's column order.
_COLUMNS: Final[dict[str, str]] = {
    "GENE ID (HGNC)": "gene",
    "GENE SYMBOL": "gene_symbol",
    "DISEASE ID (MONDO)": "disease",
    "DISEASE LABEL": "disease_label",
    "MOI": "moi",
    "SOP": "sop",
    "CLASSIFICATION": "classification",
    "CLASSIFICATION DATE": "classification_date",
    "GCEP": "gcep",
    "ONLINE REPORT": "report_url",
}
ORDER: Final = tuple(_COLUMNS.values())
_SORT_AT = (ORDER.index("gene"), ORDER.index("disease"), ORDER.index("moi"))

# The atlas's committed mirror at the bootstrap commit (atlas 4b1f8f6).
MIRROR_NAME: Final = "clingen_gene_validity.tsv"
MIRROR_SHA256: Final = "b4c89cb6b0bbb5016894661ef9add7b432b91240118e5b89e55fd3dc7e6ec1d1"
MIRROR_BANNER: Final = "FILE CREATED: 2026-08-03"
MIRROR_RETRIEVED_ON: Final = date(2026, 8, 3)
# The atlas commit whose mirror this is; the manifest URL is a permanent permalink to it.
ATLAS_COMMIT: Final = "4b1f8f68cbc30beaf932fbbc58e3617ebe401b6d"
_ATLAS_RAW = "https://raw.githubusercontent.com/chd-atlas/chd-atlas.github.io"


def read_rows(csv_path: Path) -> tuple[str, list[tuple[str, ...]]]:
    """``(banner text, rows)``: the ``FILE CREATED`` line and the rows in table order.

    Rows are sorted by ``(gene, disease, moi)``. Pure: no package is written.
    """
    with csv_path.open(encoding="utf-8", newline="") as handle:
        records = list(csv.reader(handle))
    if len(records) <= _PREAMBLE_ROWS:
        raise PackageRefused(f"{csv_path.name}: shorter than the four-line banner and a header")
    banner = records[1][0] if records[1] else ""
    if not _BANNER.fullmatch(banner):
        raise PackageRefused(
            f"{csv_path.name}: line 2 is {banner!r}, not 'FILE CREATED: YYYY-MM-DD'"
        )
    header = records[_PREAMBLE_ROWS]
    missing = sorted(set(_COLUMNS) - set(header))
    if missing:
        raise PackageRefused(f"{csv_path.name}: expected columns absent: {missing}")
    at = {name: header.index(name) for name in _COLUMNS}
    symbol_at = at["GENE SYMBOL"]
    rows: list[tuple[str, ...]] = []
    for record in records[_PREAMBLE_ROWS + 1 :]:
        if not record:
            continue
        if len(record) != len(header):
            raise PackageRefused(
                f"{csv_path.name}: a row has {len(record)} cells, the header {len(header)}"
            )
        if record[symbol_at].startswith(_SEPARATOR):
            continue
        rows.append(tuple(record[at[name]] for name in _COLUMNS))
    rows.sort(key=lambda row: tuple(row[i] for i in _SORT_AT))
    return banner, rows


def _write(
    rows: list[tuple[str, ...]],
    *,
    banner: str,
    inputs: list[dict[str, Any]],
    retrieved_on: date,
    revision: int,
    sha256: str,
    built_by: dict[str, Any],
    out: Path,
) -> Path:
    source = load_source(HERE / "source.yaml")
    release = _BANNER.fullmatch(banner)
    if release is None:
        raise PackageRefused(f"banner {banner!r} is not 'FILE CREATED: YYYY-MM-DD'")
    cells: list[Sequence[str | None]] = list(rows)
    return write_package(
        out,
        package=SOURCE,
        version=package_version(
            upstream_release=release.group(1), input_sha256=sha256, revision=revision
        ),
        ids="hgnc_id",
        coverage="every_gene",
        tables=[TableOut("clingen_validity", 1, ORDER, cells, "gene")],
        inputs=inputs,
        input_counts={"rows": len(rows), "genes": len({row[0] for row in rows})},
        dropped=[],
        built_by=built_by,
        source={
            "name": source["name"],
            "version": banner,
            "retrieved_on": retrieved_on.isoformat(),
            "url": source["url"],
            "licence": source["licence"],
            "redistribution": source["redistribution"],
        },
        licence_evidence=source["licence_evidence"],
    )


def build(
    csv_path: Path,
    *,
    retrieved_on: date,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Write the ``clingen`` package from a live download into the empty ``out``.

    ``revision`` has no default (F13): 1 for the first build of the input, 2, 3, ...
    for a recipe change or a rebuild whose ``manifest.json`` bytes differ.
    """
    built_by, _ = stamp(built_by, commit_of)
    banner, rows = read_rows(csv_path)
    sha256 = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    return _write(
        rows,
        banner=banner,
        inputs=[
            {
                "name": csv_path.name,
                "url": DOWNLOAD_URL,
                "sha256": sha256,
                "retrieved_on": retrieved_on.isoformat(),
                "upstream_release": banner,
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
    return _write(
        rows,
        banner=MIRROR_BANNER,
        inputs=[
            {
                "name": MIRROR_NAME,
                "url": f"{_ATLAS_RAW}/{atlas_commit}/mirrors/{MIRROR_NAME}",
                "sha256": sha256,
                "retrieved_on": MIRROR_RETRIEVED_ON.isoformat(),
                "upstream_release": MIRROR_BANNER,
            }
        ],
        retrieved_on=MIRROR_RETRIEVED_ON,
        revision=revision,
        sha256=sha256,
        built_by=built_by,
        out=out,
    )
