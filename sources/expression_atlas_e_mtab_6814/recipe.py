"""The bulk developmental expression profiles as the ``expression_atlas_e_mtab_6814`` package.

Moved from the atlas's ``scripts/convert_cardoso_moreira.py``; the measurements
below are that docstring's unless a line says otherwise.

Source: **E-MTAB-6814**, the EMBL-EBI Expression Atlas reprocessing of
Cardoso-Moreira et al. 2019 (*Nature*) -- human developmental transcriptomes
across seven organs.

**The published values are the Atlas's, not the paper's, and that distinction
is load-bearing.** ``E-MTAB-6814-analysis-methods.tsv`` records the Atlas
pipeline: iRAP 1.0.1, HISAT2 2.1.0 against Ensembl release 107, FeatureCounts
1.6.2, then limma quantile normalisation within each set of biological
replicates. The paper's own supplementary tables report RPKM from a different
pipeline. The two are not interchangeable and the provenance must name which
produced a number. (A stale note survives in the atlas repository's history: the
commit adding ``rpkm`` to the unit vocabulary says the source "publishes RPKM" --
true of the paper, false of the Atlas matrix this recipe reads. ``tpm`` is what
is emitted here.)

Why TPM rather than FPKM, both of which the Atlas publishes: TPM sums to a
constant per sample, so it is the more defensible quantity to compare across
samples. FPKM's across-gene sum varies between samples.

What the Atlas gives us for free
--------------------------------
Each matrix cell is not a scalar but a five-value quantile summary of that
gene's assay group -- ``min,q25,median,q75,max`` (e.g. ``8,15,18,20,21``). So
``median_abundance``, ``q25`` and ``q75`` are read off directly rather than
computed, and ``configuration.xml`` maps each assay-group id to a
``"<stage>; <organ>"`` label plus the assays it aggregates, which yields
``tissue``, ``stage`` and ``n_samples``.

The two outputs, and why they are different populations
-------------------------------------------------------
* **profiles** -- since the recipe moved here, EVERY gene of the matrix (58,735
  Ensembl ids), under its Ensembl id (``ids: ensembl_gene_id``). The converter
  carried only the atlas's registry (185 genes, as of atlas 4b1f8f6) and joined on a registry
  column; the package has no registry, so ``chd-atlas ingest`` resolves the
  Ensembl ids against its own genes table and keeps its registry's genes.
* **profile_quantiles** -- breakpoints over the matrix's **entire** gene
  universe, because a percentile's denominator is the whole transcriptome.
  Restricting it to the registry would publish a rank against a population no
  reader means. This was already so in the converter, and still is.

Gene identity is the Ensembl id, never a symbol (the matrix's second column,
``Gene Name``, is not read).

Where it runs
-------------
On the workstation, not the cluster (plan decision 10 of the phase 3 plan): the
atlas's old ``scripts/hpc/cardoso_mirrors.sbatch`` was an 11-second standard
library pass over a public download, and is not ported. The profiles table holds
6,989,465 rows, which ``write_package`` keeps in memory and checks cell by cell
in Python; measured at 4.49 GB peak RSS on the workstation (the README has the build's
measurements).
"""

from __future__ import annotations

import csv
import hashlib
import xml.etree.ElementTree as ET
from collections.abc import Iterator, Sequence
from datetime import date
from pathlib import Path
from typing import Any, Final

from feeds.package import PackageRefused, TableOut, write_package
from feeds.provenance import CommitOf, current_commit, stamp
from feeds.source import load_source
from feeds.version import package_version

HERE = Path(__file__).resolve().parent

SOURCE: Final = "expression_atlas_e_mtab_6814"
DATASET: Final = "E-MTAB-6814"
RETRIEVED_ON: Final = date(2026, 8, 19)

MATRIX_NAME: Final = "E-MTAB-6814-tpms.tsv"
CONFIGURATION_NAME: Final = "E-MTAB-6814-configuration.xml"
_FTP: Final = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-6814/"
MATRIX_URL: Final = _FTP + MATRIX_NAME
CONFIGURATION_URL: Final = _FTP + CONFIGURATION_NAME
MATRIX_SHA256: Final = "92484a979fb8e262fe1f535b2b902ba8088310bbbec76969ba873b51d593f4be"
CONFIGURATION_SHA256: Final = "31de8dd2c78067cba2a333527d1af123e9d497443c4acdf7c5dea2da115b166f"

# The quantile estimator, named because the method must travel with the number, and
# implemented here rather than delegated to a library so the definition is auditable
# from this file alone. `linear` is the R type-7 / numpy-default rule: with the
# values sorted ascending and h = (n - 1) * p, the estimate is a linear
# interpolation between the values at floor(h) and ceil(h). Pinned against worked
# literals in the tests -- a reproducibility recipe must not inherit a quantile
# convention from a dependency's version.
QUANTILE_ESTIMATOR = "linear"
PERCENTILES = list(range(101))
UNIT = "tpm"

PROFILE_COLUMNS: Final = (
    "dataset",
    "gene",
    "tissue",
    "stage",
    "median_abundance",
    "unit",
    "q25",
    "q75",
    "n_samples",
)
QUANTILE_COLUMNS: Final = (
    "dataset",
    "tissue",
    "stage",
    "percentile",
    "value",
    "unit",
    "n_genes",
)


class AssayGroup:
    """One assay group: a (stage, tissue) cell and the assays it aggregates."""

    __slots__ = ("group_id", "stage", "tissue", "n_samples")

    def __init__(self, group_id: str, stage: str, tissue: str, n_samples: int) -> None:
        self.group_id = group_id
        self.stage = stage
        self.tissue = tissue
        self.n_samples = n_samples


def parse_assay_groups(config_xml: Path) -> dict[str, AssayGroup]:
    """Map each assay-group id to its (stage, tissue) and replicate count.

    Labels read ``"10 week post conception; heart"``. Split on the first
    semicolon only: an organ name never contains one, and a stage token might.
    A label without a semicolon REFUSES: it is a group the experiment defines on
    some axis other than stage-and-organ, and guessing which half is which would
    silently mislabel a whole column. (The converter skipped and named such a
    group; measured 2026-11-12, all 119 groups of the pinned file parse, so the
    refusal costs nothing today and cannot drop a column unseen later.)
    """
    groups: dict[str, AssayGroup] = {}
    root = ET.parse(config_xml).getroot()
    for element in root.iter("assay_group"):
        group_id = element.get("id")
        label = element.get("label")
        if not group_id or not label:
            raise PackageRefused(
                f"an assay_group has id {group_id!r} and label {label!r}; both are required"
            )
        if group_id in groups:
            raise PackageRefused(f"assay group {group_id} is defined twice")
        if ";" not in label:
            raise PackageRefused(
                f"assay group {group_id}: label {label!r} has no '; ' split between stage and organ"
            )
        stage, tissue = (part.strip() for part in label.split(";", 1))
        n_samples = sum(1 for _ in element.iter("assay"))
        groups[group_id] = AssayGroup(group_id, stage, tissue, n_samples)
    return groups


def is_absent(cell: str) -> bool:
    """A blank, ``NA`` or ``NaN`` cell: the gene was not quantified in the group.

    Distinct from a malformed cell, which ``parse_cell`` also maps to None but the
    build refuses (measured 2026-11-12: 0 of the pinned matrix's cells are malformed).
    """
    cell = cell.strip()
    return not cell or cell.upper() in {"NA", "NAN"}


def parse_cell(cell: str) -> tuple[float, float, float] | None:
    """(q25, median, q75) from the Atlas's ``min,q25,median,q75,max`` quintuple.

    Returns None for an absent measurement rather than zero (and for a
    malformed cell too, which ``build`` tells apart with ``is_absent`` and
    refuses): a gene the Atlas
    did not quantify in a group is not a gene measured at zero, and collapsing
    the two would publish a floor-level abundance the source never reported.

    A cell with a single value is accepted as a degenerate quintuple (all five
    equal), which is what the Atlas emits for a one-assay group.
    """
    cell = cell.strip()
    if not cell or cell.upper() in {"NA", "NAN"}:
        return None
    parts = cell.split(",")
    try:
        values = [float(part) for part in parts]
    except ValueError:
        return None
    if len(values) == 1:
        return (values[0], values[0], values[0])
    if len(values) != 5:
        return None
    return (values[1], values[2], values[3])


def profile_row(
    dataset: str, gene: str, group: AssayGroup, q25: float, median: float, q75: float
) -> dict[str, object]:
    """One ``profiles`` row.

    ``q25``/``q75`` are emitted as null below three replicates: quartiles of two
    points are not quartiles, and a median of two is a mean of two with none of
    the robustness the word implies.
    """
    has_quartiles = group.n_samples >= 3
    return {
        "dataset": dataset,
        "gene": gene,
        "tissue": group.tissue,
        "stage": group.stage,
        "median_abundance": median,
        "unit": UNIT,
        "q25": q25 if has_quartiles else None,
        "q75": q75 if has_quartiles else None,
        "n_samples": group.n_samples,
    }


def quantile_linear(sorted_values: Sequence[float], percentile: int) -> float:
    """The R type-7 / numpy-default quantile of an already-sorted sequence.

    Written out rather than imported so the estimator this package publishes is
    readable in the file that produced the numbers. ``sorted_values`` must be
    ascending; sorting inside would hide an O(n log n) cost per percentile.
    """
    n = len(sorted_values)
    if n == 1:
        return float(sorted_values[0])
    position = (n - 1) * (percentile / 100.0)
    lower = int(position)
    upper = min(lower + 1, n - 1)
    weight = position - lower
    return float(sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight)


def quantile_rows(
    dataset: str, group: AssayGroup, medians: list[float]
) -> Iterator[dict[str, object]]:
    """101 breakpoints for one (tissue, stage) cell, over the whole universe.

    Computed on every gene the matrix quantified in this group -- the denominator
    of a published percentile is the source's own gene universe, and ``n_genes``
    states it.
    """
    ordered = sorted(medians)
    n_genes = len(ordered)
    for percentile in PERCENTILES:
        value = quantile_linear(ordered, percentile)
        yield {
            "dataset": dataset,
            "tissue": group.tissue,
            "stage": group.stage,
            "percentile": percentile,
            "value": round(float(value), 6),
            "unit": UNIT,
            "n_genes": n_genes,
        }


def _cell(value: object) -> str | None:
    return None if value is None else str(value)


def _check_pin(path: Path, expected: str) -> str:
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise PackageRefused(
            f"{path.name} has sha256 {actual}, not the pinned {expected}; the recipe was "
            "written against that file"
        )
    return actual


def build(
    matrix: Path,
    configuration: Path,
    *,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Write the ``expression_atlas_e_mtab_6814`` package into the empty ``out``.

    ``matrix`` is ``E-MTAB-6814-tpms.tsv`` and ``configuration`` is
    ``E-MTAB-6814-configuration.xml``; both are refused unless their sha256 is
    the pinned one. ``revision`` has no default: a published version is
    permanent, so the operator states whether this is the first build of the
    input (1) or a recipe change or a rebuild whose ``manifest.json`` bytes
    differ (2, 3, ...).
    """
    built_by, _ = stamp(built_by, commit_of)
    matrix_sha = _check_pin(matrix, MATRIX_SHA256)
    configuration_sha = _check_pin(configuration, CONFIGURATION_SHA256)
    source = load_source(HERE / "source.yaml")

    groups = parse_assay_groups(configuration)
    collected: dict[str, list[float]] = {group_id: [] for group_id in groups}
    profiles: list[list[str | None]] = []
    seen: set[str] = set()
    with matrix.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        columns = {index: header[index] for index in range(2, len(header))}
        unknown = sorted(g for g in columns.values() if g not in groups)
        if unknown:
            raise PackageRefused(
                f"{matrix.name}: columns {unknown[:5]} are not assay groups of the configuration"
            )
        if len(set(columns.values())) != len(columns):
            raise PackageRefused(f"{matrix.name}: a column id is repeated in the header")
        unfilled = sorted(set(groups) - set(columns.values()))
        if unfilled:
            raise PackageRefused(
                f"configured groups {unfilled[:5]} have no column in {matrix.name}"
            )
        for record in reader:
            if len(record) != len(header):
                raise PackageRefused(
                    f"{matrix.name}: line {reader.line_num} has {len(record)} fields, "
                    f"expected {len(header)}"
                )
            gene = record[0].strip()
            if not gene:
                raise PackageRefused(f"{matrix.name}: a row has a blank gene id")
            if gene in seen:
                raise PackageRefused(f"{matrix.name}: gene {gene} appears on two rows")
            seen.add(gene)
            for index, group_id in columns.items():
                group = groups[group_id]
                if is_absent(record[index]):
                    continue
                parsed = parse_cell(record[index])
                if parsed is None:
                    raise PackageRefused(
                        f"{matrix.name}: gene {gene}, column {group_id}: cell "
                        f"{record[index]!r} is neither absent nor a quintuple"
                    )
                q25, median, q75 = parsed
                collected[group_id].append(median)
                row = profile_row(DATASET, gene, group, q25, median, q75)
                profiles.append([_cell(row[c]) for c in PROFILE_COLUMNS])

    quantiles: list[list[str | None]] = []
    for group_id, group in groups.items():
        medians = collected[group_id]
        if not medians:
            raise PackageRefused(
                f"assay group {group_id} quantified no gene; it would have no grid"
            )
        for row in quantile_rows(DATASET, group, medians):
            quantiles.append([_cell(row[c]) for c in QUANTILE_COLUMNS])

    # The schemas' sort keys; the percentile is numeric, everything else text.
    profiles.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    quantiles.sort(key=lambda r: (r[0], r[1], r[2], int(r[3] or "0")))

    version = package_version(upstream_release=None, input_sha256=matrix_sha, revision=revision)
    return write_package(
        out,
        package=SOURCE,
        version=version,
        ids="ensembl_gene_id",
        coverage="every_gene",
        tables=[
            TableOut("profiles", 1, PROFILE_COLUMNS, profiles, "gene", shard=DATASET),
            TableOut("profile_quantiles", 1, QUANTILE_COLUMNS, quantiles, None, shard=DATASET),
        ],
        inputs=[
            {
                "name": MATRIX_NAME,
                "url": MATRIX_URL,
                "sha256": matrix_sha,
                "retrieved_on": RETRIEVED_ON.isoformat(),
                "upstream_release": None,
            },
            {
                "name": CONFIGURATION_NAME,
                "url": CONFIGURATION_URL,
                "sha256": configuration_sha,
                "retrieved_on": RETRIEVED_ON.isoformat(),
                "upstream_release": None,
            },
        ],
        input_counts={"rows": len(seen), "genes": len(seen)},
        dropped=[],
        built_by=built_by,
        source={
            "name": source["name"],
            "version": f"sha256:{matrix_sha}",
            "retrieved_on": RETRIEVED_ON.isoformat(),
            "url": source["url"],
            "licence": source["licence"],
            "redistribution": source["redistribution"],
        },
        licence_evidence=source["licence_evidence"],
    )
