"""The HGNC package: every approved gene, in the atlas's ``genes`` v2 table.

Moved from the atlas's ``scripts/convert_hgnc.py`` (its docstring has the
measurements behind each rule). Two inputs: hvantk's ``hgnc:lookup`` table
exported as gzipped TSV (the rows), and HGNC's complete set
(``hgnc_complete_set.txt``, the input that is hashed and counted).

Rules, each of which guards a silent wrong claim:

* ``mane_select`` is HGNC's ``ENST...|NM_...`` pair; the column holds the RefSeq
  half, chosen by its ``NM_`` prefix and not by position, so a pair without one
  yields null rather than an Ensembl id in a column read as RefSeq.
* ``uniprot`` is written verbatim, a list included (the atlas's decision 3),
  except for the genes in ``UNIPROT_PINS``: HGNC orders a list alphabetically
  and not by which protein is the gene's own, so position decides nothing. A pin
  whose accession HGNC no longer lists for that gene refuses the build, and so
  does a pinned gene missing from the table.
* The rows hvantk dropped (it keeps ``Approved`` genes only) must equal the
  complete set's rows with another status, gene for gene. They are published in
  ``dropped`` with their count, or the build refuses. Measured 2026-11-11 (review
  of this recipe): the live set (45,387 rows) and the 2026-10-06 archive (45,233
  rows) hold only Approved rows, so ``dropped`` is ``[]`` on real data and this
  guard is exercised by fixtures alone.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Mapping, Sequence
from datetime import date
from pathlib import Path
from typing import Any

import polars as pl

from feeds.package import PackageRefused, TableOut, write_package
from feeds.provenance import (
    CommitOf,
    current_commit,
    hvantk_sidecar,
    require_hvantk_release,
    stamp,
)
from feeds.source import load_source
from feeds.version import package_version

HERE = Path(__file__).resolve().parent

# The hvantk fields this recipe reads, in the `fields` list the sbatch passes to
# `hvantk reprocess`. At hvantk 8683258d the export flattens the struct-array column
# `lsdb` instead of refusing it; the selection stays so the package input carries only
# what this recipe reads. test_recipe checks the sbatch carries exactly these.
HVANTK_FIELDS = (
    "hgnc_id",
    "gene_symbol",
    "gene_name",
    "alias_symbols",
    "prev_symbols",
    "ensembl_gene_id",
    "entrez_id",
    "location",
    "uniprot_ids",
    "mane_select",
)

# The hvantk fields that are arrays. hvantk's TSV writer joins an array's elements with
# `|` after writing `%` as `%25` and `|` as `%7C` inside each element, so a reader splits
# on `|` and percent-decodes each element (hvantk/core/io/_tsv.py at 8683258d).
HVANTK_ARRAY_FIELDS = ("alias_symbols", "prev_symbols", "uniprot_ids", "mane_select")

# `genes` v2 columns in the atlas's order, each with the hvantk field it copies.
GENES_COLUMNS = (
    "hgnc_id",
    "symbol",
    "name",
    "aliases",
    "prev_symbols",
    "ensembl_gene",
    "ncbi_gene",
    "locus",
    "uniprot",
    "mane_select",
)
_FROM_HVANTK = dict(
    zip(
        GENES_COLUMNS,
        (
            "hgnc_id",
            "gene_symbol",
            "gene_name",
            "alias_symbols",
            "prev_symbols",
            "ensembl_gene_id",
            "entrez_id",
            "location",
            "uniprot_ids",
            "mane_select",
        ),
        strict=True,
    )
)

# The gene's own protein, for every gene whose `uniprot_ids` lists more than one
# accession and whose choice a human has made.
UNIPROT_PINS: Mapping[str, str] = {
    # RBM10. UniProt REST, 2026-10-25 (as recorded in the atlas's convert_hgnc.py):
    # P98175 is RBM10_HUMAN, "RNA-binding protein 10", reviewed, 930 aa, gene
    # RBM10. P0DW28 is MNS60_HUMAN, "Ribosome biogenesis inhibitor MINAS-60",
    # reviewed, 130 aa, also gene RBM10 -- "the product of an alternative open
    # reading frame (alt-ORF) of transcripts coding for the RBM10 (AC P98175)
    # protein". HGNC lists both, "P0DW28|P98175", so first-of-list is the
    # microprotein.
    "HGNC:9896": "P98175",
}

HPC_URL_FORMAT = "https://github.com/chd-atlas/feeds/tree/{commit}/sources/hgnc/hpc"

# The complete-set column each hvantk field is renamed from (HGNC_GENE_FIELDS in
# hvantk/skills/hgnc/shared/constants.py at 8683258d); the other fields keep their name.
_RAW_COLUMN = {
    "gene_symbol": "symbol",
    "gene_name": "name",
    "alias_symbols": "alias_symbol",
    "prev_symbols": "prev_symbol",
}
_MISMATCH_EXAMPLES = 5

_APPROVED = "Approved"
_ISO_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
_ARCHIVE_NAME = re.compile(r"hgnc_complete_set_(\d{4}-\d{2}-\d{2})\.txt")
_EXAMPLES = 5


_PERCENT = re.compile(r"%(.?.?)")
_DECODED = {"25": "%", "7C": "|"}


def _decode_cell(cell: str | None, column: str, gene: str | None) -> str | None:
    """Percent-decode each ``|``-separated element of an array cell, refusing ambiguity.

    A ``%`` not followed by ``25`` or ``7C`` is malformed; an element that decodes to
    something holding ``|`` would make the atlas's ``|``-joined cell ambiguous.
    """
    if cell is None:
        return None
    elements = []
    for element in cell.split("|"):
        for escape in _PERCENT.findall(element):
            if escape not in _DECODED:
                raise PackageRefused(
                    f"{gene}, column {column}: '%{escape}' is not a hvantk escape (%25, %7C)"
                )
        decoded = _PERCENT.sub(lambda m: _DECODED[m.group(1)], element)
        if "|" in decoded:
            raise PackageRefused(
                f"{gene}, column {column}: element {element!r} decodes to text holding '|', "
                "which would make the |-joined cell ambiguous"
            )
        elements.append(decoded)
    return "|".join(elements)


def _read_hvantk(path: Path) -> pl.DataFrame:
    """The hvantk export as strings, unquoted and untyped: no cell is retyped."""
    frame = pl.read_csv(path, separator="\t", quote_char=None, infer_schema=False)
    missing = [f for f in HVANTK_FIELDS if f not in frame.columns]
    if missing:
        raise PackageRefused(
            f"{path.name}: the hvantk table lacks column(s) {', '.join(missing)}; "
            f"has {', '.join(frame.columns)}"
        )
    for column in HVANTK_ARRAY_FIELDS:
        frame = frame.with_columns(
            pl.Series(
                column,
                [
                    _decode_cell(cell, column, gene)
                    for cell, gene in zip(frame[column], frame["hgnc_id"], strict=True)
                ],
                dtype=pl.String,
            )
        )
    return frame


def _read_raw(path: Path) -> pl.DataFrame:
    """The complete set with polars' default quoting, as convert_hgnc.py read it.

    Its multi-value cells are quoted (``"P0DW28|P98175"``), unlike the hvantk export.
    """
    frame = pl.read_csv(path, separator="\t", infer_schema_length=0)
    needed = ["hgnc_id", "status", *(_RAW_COLUMN.get(f, f) for f in HVANTK_FIELDS)]
    missing = sorted({c for c in needed if c not in frame.columns})
    if missing:
        raise PackageRefused(f"{path.name}: the complete set lacks column(s) {', '.join(missing)}")
    return frame


def _tokens(cell: str | None) -> str | None:
    """A raw complete-set cell as hvantk tokenizes a pipe-separated field.

    hvantk/skills/hgnc/builder.py ``_split_pipe_tokens`` at 8683258d is
    ``expr.replace('"', "").split("\\|").map(strip).filter(!= "")``, giving an empty
    array for an empty or missing cell. Returned ``|``-joined, or None for no token.
    """
    if cell is None or cell == "":
        return None
    parts = [t.strip() for t in cell.replace('"', "").split("|")]
    return "|".join(t for t in parts if t != "") or None


def _cross_check(frame: pl.DataFrame, raw: pl.DataFrame) -> None:
    """Refuse when a cell of the hvantk table is not what the complete set says.

    hvantk is a dependency, not an authority: its 8683258d export once double-encoded
    35 non-ASCII aliases (``RORα`` as ``RORÎ±``) and every other check passed. Each
    cell the recipe reads (after percent-decoding, before the recipe's own rules) is
    compared with the same gene's raw field: array fields normalised as hvantk
    documents (``_tokens``), scalars as the raw text.
    """
    approved = raw.filter(pl.col("status").fill_null("") == _APPROVED)
    by_gene = {r["hgnc_id"]: r for r in approved.iter_rows(named=True)}
    found: list[str] = []
    count = 0
    for record in frame.iter_rows(named=True):
        theirs = by_gene.get(record["hgnc_id"])
        if theirs is None:
            continue  # `_dropped` has already refused a gene the complete set lacks
        for field in HVANTK_FIELDS:
            wanted = theirs[_RAW_COLUMN.get(field, field)]
            if field in HVANTK_ARRAY_FIELDS:
                wanted = _tokens(wanted)
            if record[field] != wanted:
                count += 1
                if len(found) < _MISMATCH_EXAMPLES:
                    found.append(
                        f"{record['hgnc_id']} {field}: hvantk {record[field]!r}, "
                        f"complete set {wanted!r}"
                    )
    if count:
        raise PackageRefused(
            f"{count} cell(s) of the hvantk table differ from the complete set "
            f"(showing {len(found)}): " + "; ".join(found)
        )


def _refseq(cell: str | None) -> str | None:
    """The first ``NM_`` token of HGNC's pipe-joined MANE Select pair, else null."""
    if cell is None:
        return None
    return next((t for t in cell.split("|") if t.startswith("NM_")), None)


def _uniprot(gene: str, symbol: str | None, cell: str | None) -> str | None:
    pinned = UNIPROT_PINS.get(gene)
    if pinned is None:
        return cell
    listed = cell.split("|") if cell is not None else []
    if pinned not in listed:
        raise PackageRefused(
            f"{gene} ({symbol}): the pinned UniProt accession {pinned} is not among those "
            f"HGNC lists ({cell}); re-read UniProt and update UNIPROT_PINS"
        )
    return pinned


def _dropped(raw: pl.DataFrame, kept: set[str], rows_kept: int) -> list[dict[str, Any]]:
    """One entry per status the complete set holds besides Approved, checked against hvantk.

    The table must hold exactly the Approved ids of the complete set: a gene
    missing from it, one it holds that is not Approved, or a repeat would each
    publish a package that is not "every approved gene".
    """
    status = raw["status"].fill_null("")
    ids = raw["hgnc_id"]
    approved = {i for i, s in zip(ids, status, strict=True) if s == _APPROVED}
    approved.discard(None)
    if len(kept) != rows_kept:
        raise PackageRefused(
            f"the hvantk table repeats an hgnc_id ({rows_kept} rows, {len(kept)} ids)"
        )
    absent = sorted(approved - kept)
    extra = sorted(kept - approved)
    if absent or extra:
        raise PackageRefused(
            "the hvantk table is not the complete set's Approved genes: "
            f"{len(absent)} missing (e.g. {absent[:3]}), {len(extra)} not Approved there "
            f"(e.g. {extra[:3]})"
        )
    entries: list[dict[str, Any]] = []
    # dict.fromkeys keeps raw-file order and no set decides it; `sorted` then fixes the output.
    for label in sorted(x for x in dict.fromkeys(status.to_list()) if x != _APPROVED):
        gone = sorted({i for i, s in zip(ids, status, strict=True) if s == label and i is not None})
        shown = label or "(blank)"
        entries.append(
            {
                "reason": f"status is {shown}, not Approved; hvantk keeps Approved genes only",
                "count": int((status == label).sum()),
                "examples": gone[:_EXAMPLES],
            }
        )
    if sum(e["count"] for e in entries) != raw.height - rows_kept:
        raise PackageRefused(
            f"{raw.height - rows_kept} complete-set rows are absent from the hvantk table but "
            f"{sum(e['count'] for e in entries)} are not Approved; they must be the same rows"
        )
    return entries


def build(
    hvantk_tsv: Path,
    raw: Path,
    *,
    raw_url: str,
    retrieved_on: date,
    release: str,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Write the ``hgnc`` package into the empty directory ``out`` and return it.

    ``release`` is the complete set's archive date (e.g. ``2026-10-06``) and is
    the package version at revision 1. ``revision`` has no default: a published
    version is permanent (F13), so the operator states whether this build is the
    first for its input or a recipe change (``-r2``, ...) or a rebuild whose
    ``manifest.json`` bytes differ.

    ``built_by`` must carry ``hvantk_release``, the full commit of the hvantk
    checkout that ran, and it must equal the ``builder_commit`` of hvantk's
    sidecar ``<hvantk_tsv>.provenance.json``. ``recipe_commit`` is set by the
    build from ``commit_of`` (the clean feeds checkout; tests stub it).
    """
    if not _ISO_DATE.fullmatch(release):
        raise PackageRefused(f"release '{release}' is not an ISO date (YYYY-MM-DD)")
    named = _ARCHIVE_NAME.fullmatch(raw.name)
    if named and named.group(1) != release:
        raise PackageRefused(
            f"{raw.name} is the {named.group(1)} archive but release is '{release}'"
        )
    built_by, recipe_commit = stamp(built_by, commit_of)
    release_commit = require_hvantk_release(built_by)
    builder_commit, export_made_on = hvantk_sidecar(hvantk_tsv)
    if builder_commit != release_commit:
        raise PackageRefused(
            f"built_by hvantk_release is {release_commit!r} but hvantk's sidecar records "
            f"builder_commit {builder_commit}"
        )
    source = load_source(HERE / "source.yaml")
    frame = _read_hvantk(hvantk_tsv)
    raw_frame = _read_raw(raw)
    sha256 = hashlib.sha256(raw.read_bytes()).hexdigest()

    hgnc_ids = frame["hgnc_id"].to_list()
    if any(i is None for i in hgnc_ids):
        raise PackageRefused("the hvantk table has a row with no hgnc_id")
    dropped = _dropped(raw_frame, set(hgnc_ids), frame.height)

    wanted = set(UNIPROT_PINS) - set(hgnc_ids)
    if wanted:
        raise PackageRefused(
            f"pinned gene(s) {', '.join(sorted(wanted))} are not in the hvantk table; "
            "re-read UniProt and update UNIPROT_PINS"
        )

    _cross_check(frame, raw_frame)

    rows: list[Sequence[str | None]] = []
    for record in frame.sort("hgnc_id").iter_rows(named=True):
        cells = {column: record[field] for column, field in _FROM_HVANTK.items()}
        cells["uniprot"] = _uniprot(record["hgnc_id"], record["gene_symbol"], record["uniprot_ids"])
        cells["mane_select"] = _refseq(record["mane_select"])
        rows.append([cells[c] for c in GENES_COLUMNS])

    return write_package(
        out,
        package="hgnc",
        version=package_version(upstream_release=release, input_sha256=sha256, revision=revision),
        ids="hgnc_id",
        coverage="every_gene",
        tables=[TableOut("genes", 2, GENES_COLUMNS, rows, "hgnc_id")],
        inputs=[
            {
                "name": raw.name,
                "url": raw_url,
                "sha256": sha256,
                "retrieved_on": retrieved_on.isoformat(),
                "upstream_release": release,
            },
            {
                "name": hvantk_tsv.name,
                # Not a download location: hvantk's export is a cluster intermediate and is
                # not published. The URL pins the job code, at the commit that built this.
                "url": HPC_URL_FORMAT.format(commit=recipe_commit),
                "sha256": hashlib.sha256(hvantk_tsv.read_bytes()).hexdigest(),
                "retrieved_on": export_made_on.isoformat(),
                "upstream_release": None,
            },
        ],
        input_counts={
            "rows": raw_frame.height,
            "genes": raw_frame["hgnc_id"].drop_nulls().n_unique(),
        },
        dropped=dropped,
        built_by=built_by,
        source={
            "name": source["name"],
            "version": release,
            "retrieved_on": retrieved_on.isoformat(),
            "url": source["url"],
            "licence": source["licence"],
            "redistribution": source["redistribution"],
        },
        licence_evidence=source["licence_evidence"],
    )
