"""Farah et al. 2024 developing human heart, per-heart pseudobulk, as the ``farah2024_hoc`` package.

Moved from the atlas's ``scripts/convert_farah2024.py``; the measurements named
below are that script's docstring and tests', taken by their author on the dates
given, and were not re-taken here except where a line says so.

Input is the TSV ``hpc/farah_extract.py`` writes on the cluster: one row per
(gene, heart, cell type, region) with the group's summed counts, the number of
its cells with a non-zero count, its cell count, and its library size (the
``sum`` layer summed over every gene the source measured, not only the registry's).

Per (gene, cell_type, region) the package carries the median across hearts of
counts-per-million, with quartiles when at least three hearts support them, the
number of hearts, the cells pooled, and the share of those cells expressing the
gene. A heart is ``{age}pcw-{replicate}``: two hearts of one gestational age are
two donors.

THE SLICE
---------
The package is a registry slice: ``coverage`` is ``{"registry_slice":
{"atlas_commit": ..., "genes": ...}}``, the atlas's ``mirrors/genes.tsv`` at that
commit. The source names its genes by symbol and carries no Ensembl id, so the
slice keeps a source symbol only when it is the **approved symbol** of a
registry row (``ids: symbol``, ``gene`` = the source's symbol). An alias or
previous symbol cannot be told from an unrelated gene's approved symbol: the
atlas converter measured (2026-10-28, against its ``mirrors/genes.tsv``) that
``HGF`` is SOS1's alias and hepatocyte growth factor's symbol, ``PKD1`` is
PRKD1's alias and polycystin 1's, and ``POP1`` is POPDC1's alias and the
ribonuclease P subunit's. A source symbol that is not a registry gene's approved
symbol is left out and counted in ``dropped``. A registry gene whose approved
symbol the source lacks has no rows (CRIPTO: the source uses its previous symbol
TDGF1).

WHAT STOPS THE BUILD
--------------------
A source symbol with leading or trailing whitespace, a heart counted
twice in one group, a heart's group with fewer than 10 cells, and any malformed
input row refuse the build before a byte is written. The 10-cell rule is refused
rather than dropped: hvantk drops such a group at the job's ``--min-cells``
(``matrix_utils.py:205`` at hvantk commit 96ceafb4, per the atlas converter), so
one arriving means the job ran without it and its whole output is suspect. It is
checked on every row, whether or not the row's symbol is in the slice.

VERSION
-------
The upstream file is the 3.3 GB h5ad, which stays on the cluster, so this recipe
cannot hash it. ``h5ad_sha256`` is the ``input sha256:`` line the cluster job
prints, and the version is its first 8 hex digits. The pseudobulk file's own
sha256 is recorded as a second ``inputs`` entry.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from datetime import date
from pathlib import Path
from typing import Any, Final, NamedTuple

from feeds.package import PackageRefused, TableOut, write_package
from feeds.source import load_source
from feeds.version import package_version

HERE = Path(__file__).resolve().parent

SOURCE: Final = "farah2024_hoc"
DATASET: Final = "hoc"
H5AD_NAME: Final = "farah_2024.h5ad"
# The h5ad's build date (its uns/hvantk_metadata/build_date, per the atlas's
# mirrors/sources.yaml); the raw download date is not recorded.
H5AD_RETRIEVED_ON: Final = date(2026, 4, 18)
PSEUDOBULK_NAME: Final = "pseudobulk.tsv"
# Not a download location: the pseudobulk is a cluster-local intermediate and is not
# published. The URL pins the code that produced it, at the recipe commit that ran.
PSEUDOBULK_URL_FORMAT: Final = (
    "https://github.com/chd-atlas/feeds/tree/{commit}/sources/farah2024_hoc/hpc"
)

UNIT: Final = "cpm"
# The next four values are copied from the atlas's src/chd_atlas/tables.py. The
# atlas enforces them when it ingests a package (CTP009 holds quartiles to
# MIN_HEARTS_FOR_QUARTILES, CTP010 holds n_cells to the donors); this recipe
# applies them first so a violation stops here rather than at ingest.
CELL_TYPE_REGIONS: Final = frozenset({"LA", "RA", "LV", "RV", "IVS"})
POOLED_REGION: Final = "all"
MIN_CELLS_PER_HEART: Final = 10
MIN_HEARTS_FOR_QUARTILES: Final = 3
REGIONS: Final = CELL_TYPE_REGIONS | {POOLED_REGION}

COLUMNS: Final = (
    "dataset",
    "gene",
    "cell_type",
    "region",
    "median_abundance",
    "unit",
    "q25",
    "q75",
    "n_donors",
    "n_cells",
    "fraction_expressing",
    "source",
)
INPUT_COLUMNS: Final = (
    "symbol",
    "heart",
    "cell_type",
    "region",
    "counts",
    "count_nonzero",
    "n_cells",
    "library_size",
)
_SHA256: Final = re.compile(r"[0-9a-f]{64}")
_COMMIT: Final = re.compile(r"[0-9a-f]{40}")
_EXAMPLES: Final = 5


class PseudobulkRow(NamedTuple):
    symbol: str
    heart: str
    cell_type: str
    region: str
    counts: float
    count_nonzero: int
    n_cells: int
    library_size: float


def quantile_linear(sorted_values: Sequence[float], percentile: int) -> float:
    """The R type-7 / numpy-default quantile of an already-sorted sequence.

    Copied from the atlas's ``scripts/convert_cardoso_moreira.py``. Written out
    rather than imported so the estimator the package carries is readable in the
    file that produced the numbers. ``sorted_values`` must be ascending.
    """
    n = len(sorted_values)
    if n == 1:
        return float(sorted_values[0])
    position = (n - 1) * (percentile / 100.0)
    lower = int(position)
    upper = min(lower + 1, n - 1)
    weight = position - lower
    return float(sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight)


def _finite(text: str, what: str, where: str) -> float:
    try:
        value = float(text)
    except ValueError:
        raise PackageRefused(f"{where}: {what} {text!r} is not a number") from None
    if not math.isfinite(value):
        raise PackageRefused(f"{where}: {what} {text!r} is not finite")
    return value


def _tsv_lines(path: Path) -> list[list[str]]:
    """Every line of ``path`` split on tabs, with no quoting and no retyping."""
    text = path.read_text(encoding="utf-8")
    if text.endswith("\n"):
        text = text[:-1]
    return [line.split("\t") for line in text.split("\n")]


def read_pseudobulk(path: Path) -> list[PseudobulkRow]:
    """Every row of the extract, or PackageRefused naming the first bad one."""
    lines = _tsv_lines(path)
    header = lines[0]
    if tuple(header) != INPUT_COLUMNS:
        raise PackageRefused(f"{path}: header is {header!r}, expected {list(INPUT_COLUMNS)}")
    rows: list[PseudobulkRow] = []
    for number, record in enumerate(lines[1:], start=2):
        where = f"{path}:{number}"
        if len(record) != len(INPUT_COLUMNS):
            raise PackageRefused(f"{where}: {len(record)} fields, expected {len(INPUT_COLUMNS)}")
        symbol, heart, cell_type, region, counts, nonzero, cells, library = record
        if not (symbol and heart and cell_type):
            raise PackageRefused(f"{where}: blank symbol, heart or cell_type")
        if symbol != symbol.strip():
            raise PackageRefused(f"{where}: symbol {symbol!r} has leading or trailing whitespace")
        if region not in REGIONS:
            raise PackageRefused(f"{where}: region {region!r} is not one of {sorted(REGIONS)}")
        count_value = _finite(counts, "counts", where)
        library_value = _finite(library, "library_size", where)
        nonzero_value = _finite(nonzero, "count_nonzero", where)
        cells_value = _finite(cells, "n_cells", where)
        if count_value < 0:
            raise PackageRefused(f"{where}: counts {counts!r} is negative")
        if library_value <= 0:
            raise PackageRefused(f"{where}: library_size {library!r} must be positive")
        if nonzero_value != int(nonzero_value) or cells_value != int(cells_value):
            raise PackageRefused(f"{where}: count_nonzero and n_cells must be whole numbers")
        if not 0 <= nonzero_value <= cells_value or cells_value < 1:
            raise PackageRefused(
                f"{where}: count_nonzero {nonzero!r} must lie within 0..n_cells {cells!r}"
            )
        rows.append(
            PseudobulkRow(
                symbol,
                heart,
                cell_type,
                region,
                count_value,
                int(nonzero_value),
                int(cells_value),
                library_value,
            )
        )
    return rows


def read_registry(path: Path) -> tuple[frozenset[str], int]:
    """The approved symbols of the atlas's ``genes.tsv``, and its row count.

    Only ``hgnc_id`` and ``symbol`` are read; aliases and previous symbols are
    never part of the slice. A registry listing one id or one symbol twice is
    refused, as is a blank cell in either column.
    """
    lines = _tsv_lines(path)
    header = lines[0]
    for column in ("hgnc_id", "symbol"):
        if column not in header:
            raise PackageRefused(f"{path}: no '{column}' column in the registry")
    at_id, at_symbol = header.index("hgnc_id"), header.index("symbol")
    ids: set[str] = set()
    symbols: set[str] = set()
    for number, record in enumerate(lines[1:], start=2):
        if len(record) != len(header):
            raise PackageRefused(f"{path}:{number}: {len(record)} fields, expected {len(header)}")
        gene, symbol = record[at_id], record[at_symbol]
        if not gene.strip() or not symbol.strip():
            raise PackageRefused(f"{path}:{number}: blank hgnc_id or symbol")
        if gene in ids or symbol in symbols:
            raise PackageRefused(f"{path}:{number}: {gene} / {symbol} is listed twice")
        ids.add(gene)
        symbols.add(symbol)
    return frozenset(symbols), len(ids)


def slice_symbols(symbols: Sequence[str], approved: frozenset[str]) -> set[str]:
    """The source symbols that are a registry gene's approved symbol, exactly as written.

    ``read_pseudobulk`` has already refused a padded symbol, so none is stripped here.
    """
    return {symbol for symbol in symbols if symbol in approved}


def summarise(rows: Sequence[PseudobulkRow], kept: set[str]) -> list[dict[str, object]]:
    """One row per (gene, cell_type, region), sorted to the table's sort key.

    Rows whose symbol is not in ``kept`` are skipped. A heart's group with fewer
    than ``MIN_CELLS_PER_HEART`` cells stops the build, whatever its symbol.
    """
    cpm: dict[tuple[str, str, str], dict[str, float]] = defaultdict(dict)
    nonzero: dict[tuple[str, str, str], int] = defaultdict(int)
    cells: dict[tuple[str, str, str], int] = defaultdict(int)
    for row in rows:
        if row.n_cells < MIN_CELLS_PER_HEART:
            raise PackageRefused(
                f"group {row.heart}/{row.cell_type}/{row.region} has n_cells {row.n_cells}, "
                f"below {MIN_CELLS_PER_HEART}: hvantk drops such groups at --min-cells "
                f"{MIN_CELLS_PER_HEART}, so this pseudobulk was made without it"
            )
        if row.symbol not in kept:
            continue
        key = (row.symbol, row.cell_type, row.region)
        if row.heart in cpm[key]:
            raise PackageRefused(f"heart {row.heart} appears twice for {key}")
        cpm[key][row.heart] = row.counts / row.library_size * 1e6
        nonzero[key] += row.count_nonzero
        cells[key] += row.n_cells

    out: list[dict[str, object]] = []
    for key in sorted(cpm):
        gene, cell_type, region = key
        values = sorted(cpm[key].values())
        enough = len(values) >= MIN_HEARTS_FOR_QUARTILES
        out.append(
            {
                "dataset": DATASET,
                "gene": gene,
                "cell_type": cell_type,
                "region": region,
                "median_abundance": quantile_linear(values, 50),
                "unit": UNIT,
                "q25": quantile_linear(values, 25) if enough else None,
                "q75": quantile_linear(values, 75) if enough else None,
                "n_donors": len(values),
                "n_cells": cells[key],
                "fraction_expressing": nonzero[key] / cells[key],
                "source": SOURCE,
            }
        )
    return out


def _cell(value: object) -> str | None:
    return None if value is None else str(value)


def _table(summary: Sequence[Mapping[str, object]]) -> TableOut:
    cells: list[Sequence[str | None]] = [[_cell(row[c]) for c in COLUMNS] for row in summary]
    return TableOut("cell_type_profiles", 1, COLUMNS, cells, "gene", shard=DATASET)


def build(
    pseudobulk: Path,
    registry: Path,
    *,
    atlas_commit: str,
    h5ad_sha256: str,
    pseudobulk_made_on: date,
    built_by: dict[str, Any],
    out: Path,
) -> Path:
    """Write the ``farah2024_hoc`` package into the empty directory ``out`` and return it.

    ``registry`` is the atlas's ``mirrors/genes.tsv`` at ``atlas_commit``.
    ``h5ad_sha256`` is the h5ad's digest as the cluster job printed it.
    ``pseudobulk_made_on`` is the job's date (its log's ``date:`` line, or ``sacct``);
    ``built_by["recipe_commit"]`` pins the pseudobulk input's URL.
    """
    if not _SHA256.fullmatch(h5ad_sha256):
        raise PackageRefused("h5ad_sha256 must be 64 lowercase hex characters")
    if not _COMMIT.fullmatch(atlas_commit):
        raise PackageRefused("atlas_commit must be a 40-character lowercase hex sha")
    recipe_commit = built_by.get("recipe_commit")
    if not isinstance(recipe_commit, str) or not _COMMIT.fullmatch(recipe_commit):
        raise PackageRefused('built_by["recipe_commit"] must be a 40-character lowercase hex sha')
    source = load_source(HERE / "source.yaml")
    pseudobulk_sha256 = hashlib.sha256(pseudobulk.read_bytes()).hexdigest()
    rows = read_pseudobulk(pseudobulk)
    approved, registry_genes = read_registry(registry)
    kept = slice_symbols([r.symbol for r in rows], approved)
    summary = summarise(rows, kept)
    if not summary:
        raise PackageRefused("no pseudobulk row is a registry gene's approved symbol")

    outside = [r.symbol for r in rows if r.symbol not in kept]
    dropped: list[dict[str, Any]] = []
    if outside:
        dropped.append(
            {
                "reason": "symbol is not an approved symbol of a registry gene",
                "count": len(outside),
                "examples": sorted(set(outside))[:_EXAMPLES],
            }
        )
    version = package_version(
        upstream_release=None,
        input_sha256=h5ad_sha256,
        # A recipe change that alters this package's bytes bumps `revision` (-r2, ... per spec F13).
        revision=1,
    )
    return write_package(
        out,
        package=SOURCE,
        version=version,
        ids="symbol",
        coverage={"registry_slice": {"atlas_commit": atlas_commit, "genes": registry_genes}},
        tables=[_table(summary)],
        inputs=[
            {
                "name": H5AD_NAME,
                "url": source["url"],
                "sha256": h5ad_sha256,
                "retrieved_on": H5AD_RETRIEVED_ON.isoformat(),
                "upstream_release": None,
            },
            {
                "name": PSEUDOBULK_NAME,
                "url": PSEUDOBULK_URL_FORMAT.format(commit=recipe_commit),
                "sha256": pseudobulk_sha256,
                "retrieved_on": pseudobulk_made_on.isoformat(),
                "upstream_release": None,
            },
        ],
        input_counts={"rows": len(rows), "genes": len({r.symbol for r in rows})},
        dropped=dropped,
        built_by=built_by,
        source={
            "name": source["name"],
            "version": f"sha256:{h5ad_sha256}",
            "retrieved_on": H5AD_RETRIEVED_ON.isoformat(),
            "url": source["url"],
            "licence": source["licence"],
            "redistribution": source["redistribution"],
        },
        licence_evidence=source["licence_evidence"],
    )
