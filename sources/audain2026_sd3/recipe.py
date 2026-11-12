"""The Audain 2026 case-control supplement (Supplementary Data 3) as the ``audain2026_sd3`` package.

Moved from the atlas's ``scripts/convert_burden.py``; the measurements below are
that docstring's, taken by its author on 2026-08-04 against the real workbook,
and were not re-taken here (the real build's counts are in the README). The
package carries **every** symbol of the workbook, with ``gene`` the workbook's
``SYMBOL`` verbatim (``ids: symbol``); the atlas's ``chd-atlas ingest`` resolves
symbols (approved, then alias, then previous) and keeps its registry's genes.
The old converter joined approved symbols only and so lost CRIPTO (``TDGF1``),
POPDC1 (``BVES``) and CIROZ (``C1orf127``); the resolver recovers them.

WHAT WAS MEASURED, 2026-08-04, against the real workbook
--------------------------------------------------------
Supplementary Data 3 has four sheets: an Info sheet, then ``aCHD_vs_controls``
(46,342 rows), ``sCHD_vs_controls`` (46,081) and ``nsCHD_vs_controls`` (46,186).
Each carries all four carrier counts, and the three differ only in which
comparison the ``fet.*`` columns describe. That mapping was verified
arithmetically rather than assumed: recomputing the odds ratio as
(a/(N-a))/(c/(M-c)) agreed with the published value on 32,464 / 22,992 / 27,788
testable rows for the matching stratum with **zero** disagreements. Every
mismatched pairing disagreed on tens of thousands of rows. Only the zero is
quoted as a figure: it is stable at every tolerance from 1e-2 to 5e-2, whereas
the disagreement counts move with the tolerance, and an earlier draft of this
docstring quoted a range no single tolerance reproduces.

``maf`` is ``1E-3`` and ``analysis`` is the sheet's own stratum, both constant
within a sheet. Denominators are constant too -- 3,876 cases (1,471 syndromic +
2,405 non-syndromic) against 45,082 controls -- and do **not** match the
abstract's 4,747 / 52,881. The paper does not reconcile the gap: it says it
"assembled an exome dataset consisting of 57,628 samples (4747 CHD cases and
52,881 controls)" and never states 3,876, 1,471 or 2,405 anywhere. So the reason
for the difference is unstated, and this docstring does not invent one. What is
certain is which numbers the published odds ratios were computed from, and
those are the ones the package carries.

``Infinity`` is the only non-numeric value anywhere in the four statistic
columns: 448 + 187 + 292 = 927 rows, always in ``fet.odds_ratio`` and
``fet.ci_95_upper`` together, where no control carried. Those rows hold the
study's strongest results, so they are published with a null effect and
``effect_bound = unbounded_above`` rather than dropped or given an invented
ceiling. Zero rows have an inverted confidence interval.

EXCEL HAS CORRUPTED 25 GENE SYMBOLS IN THIS WORKBOOK
-----------------------------------------------------
``SYMBOL`` holds 25 distinct values matching ``\\d{4}-\\d{2}-\\d{2}`` --
``2001-03-01``, ``2001-09-01``, ``2001-12-01`` and so on. These are the classic
MARCH/SEPT/DEC autocorrection: a spreadsheet read ``MARCH1`` as a date. It is
irreversible from this file, because ``MARCH1`` and ``MARC1`` collapse to the
same cell value, which is why four (SYMBOL, csq_group) pairs are duplicated in
each sheet.

None of the 25 resolved to an in-scope gene in the atlas, so nothing it
publishes was lost. They are published in ``dropped`` ("Excel date corruption",
with the row count and up to five examples) and never guessed: if a future
workbook corrupts a symbol that *is* in scope, that gene's entire burden
evidence would vanish with a green build -- this project's characteristic
failure exactly. A silent drop is the thing to prevent; a named one is a
curator's decision.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Sequence
from datetime import date
from pathlib import Path
from typing import Any, Final

from feeds.numbers import shortest
from feeds.package import PackageRefused, TableOut, write_package
from feeds.provenance import CommitOf, current_commit, stamp
from feeds.source import load_source
from feeds.version import package_version
from feeds.xlsx import read_sheet

HERE = Path(__file__).resolve().parent

WORKBOOK_NAME: Final = "41525_2026_582_MOESM4_ESM.xlsx"
PINNED_SHA256: Final = "8ebb6fdceb38b0099d6229a1ebf4d77db9f9492a62cfafb569d2865f9645ee2b"
RETRIEVED_ON: Final = date(2026, 8, 4)

# The publication, the source registry entry, and the cohorts -- identifiers the
# atlas's `curation/` and `mirrors/` agree with, named once here rather than
# spelled inline per row.
STUDY: Final = "PMID:42230622"
SOURCE: Final = "audain2026_sd3"
CASE_COHORTS: Final = "cnchd;ddd;nottingham"
CONTROL_COHORTS: Final = "ukbb"

# sheet name -> (cohort_stratum, carrier column, denominator column). The
# carrier/denominator pair is what the arithmetic check above pinned.
SHEETS: Final[dict[str, tuple[str, str, str]]] = {
    "aCHD_vs_controls": ("all", "n_het_cases", "n_total_cases"),
    "sCHD_vs_controls": ("syndromic", "n_het_syndromic", "n_total_syndromic"),
    "nsCHD_vs_controls": ("nonsyndromic", "n_het_nonsyndromic", "n_total_nonsyndromic"),
}

# The value each sheet's own `analysis` column holds on every row (measured
# constant per sheet on the pinned workbook). Checked, so a sheet read under the
# wrong name cannot pair another comparison's statistics with this stratum's counts.
_ANALYSIS: Final[dict[str, str]] = {
    "aCHD_vs_controls": "all_cases",
    "sCHD_vs_controls": "syndromic",
    "nsCHD_vs_controls": "nonsyndromic",
}

# The paper's `csq_group` vocabulary, mapped onto the atlas's `ConsequenceClass`.
# `syn` is carried through rather than filtered: synonymous variants should show
# no enrichment, so the row is the study's own negative control and belongs
# beside the loss-of-function row it calibrates.
_CONSEQUENCE: Final[dict[str, str]] = {
    "hcLOF": "lof",
    "missC": "missense_damaging",
    "syn": "synonymous",
}

_DATE_CORRUPTED: Final = re.compile(r"\d{4}-\d{2}-\d{2}")
_EXAMPLES: Final = 5

# Every column the row builder reads, checked against each sheet's header.
_COLUMNS: Final = (
    "SYMBOL",
    "analysis",
    "csq_group",
    "maf",
    "n_het_controls",
    "n_total_controls",
    "fet.odds_ratio",
    "fet.ci_95_lower",
    "fet.ci_95_upper",
    "fet.p_value",
)

# The atlas's `burden` v1 column order (`burden_mirror.ORDER`), equal to the
# vendored schemas/atlas/tables/burden.v1.json; the writer checks every name.
ORDER: tuple[str, ...] = (
    "study",
    "gene",
    "cohort_stratum",
    "lesion_group",
    "variant_class",
    "consequence_class",
    "origin",
    "maf_max",
    "count_unit",
    "n_case_carriers",
    "n_cases",
    "n_chd_case_carriers",
    "comparator",
    "n_control_carriers",
    "n_controls",
    "expected_count",
    "effect",
    "effect_measure",
    "effect_bound",
    "ci_low",
    "ci_high",
    "pvalue",
    "pvalue_test",
    "pvalue_adjusted",
    "pvalue_adjustment",
    "case_cohorts",
    "control_cohorts",
    "method_note",
    "source",
)

# `TABLE_SCHEMAS["burden"].sort_key`, which the atlas's `validate_sort_order` enforces.
_SORT_KEY: tuple[str, ...] = (
    "gene",
    "study",
    "cohort_stratum",
    "lesion_group",
    "variant_class",
    "consequence_class",
    "origin",
    "maf_max",
)


def read_records(workbook: Path, sheet: str) -> list[dict[str, str]]:
    """One stratum's sheet as header-keyed records, refusing a sheet that lacks a needed column."""
    try:
        records = read_sheet(workbook, sheet)
    except ValueError as exc:
        raise PackageRefused(f"{workbook.name}: {exc}") from exc
    _, carriers, denominator = SHEETS[sheet]
    present = set(records[0]) if records else set()
    missing = [c for c in (*_COLUMNS, carriers, denominator) if c not in present]
    if missing:
        raise PackageRefused(f"sheet {sheet!r} has no column {', '.join(missing)}")
    expected = _ANALYSIS[sheet]
    for number, record in enumerate(records, start=2):
        if record["analysis"] != expected:
            raise PackageRefused(
                f"sheet {sheet!r} row {number}: analysis is {record['analysis']!r}, "
                f"not {expected!r}; the sheet does not hold the comparison its name says"
            )
    return records


def row(record: dict[str, str], stratum: str, carriers: str, denominator: str) -> dict[str, str]:
    """One burden row from one sheet record; ``gene`` is the workbook's symbol, verbatim."""
    odds_ratio = record["fet.odds_ratio"]
    unbounded = odds_ratio == "Infinity"
    csq = record["csq_group"]
    if csq not in _CONSEQUENCE:
        raise PackageRefused(
            f"csq_group {csq!r} (symbol {record['SYMBOL']!r}) is not one of "
            f"{', '.join(sorted(_CONSEQUENCE))}"
        )
    return {
        "study": STUDY,
        "gene": record["SYMBOL"],
        "cohort_stratum": stratum,
        "lesion_group": "",
        "variant_class": "snv_indel",
        "consequence_class": _CONSEQUENCE[csq],
        # The analysis collapsed rare variants regardless of how they were
        # inherited; it is not a de novo test.
        "origin": "any",
        "maf_max": shortest(record["maf"]),
        # The supplement's `n_het_*` columns count individuals carrying at least
        # one qualifying variant, which is why this study's denominators are its
        # sample sizes rather than twice them.
        "count_unit": "individuals",
        "n_case_carriers": record[carriers],
        "n_cases": record[denominator],
        # Every case is a CHD case, so there is no split to give.
        "n_chd_case_carriers": "",
        "comparator": "control_cohort",
        "n_control_carriers": record["n_het_controls"],
        "n_controls": record["n_total_controls"],
        "expected_count": "",
        "effect": "" if unbounded else shortest(odds_ratio),
        "effect_measure": "odds_ratio",
        "effect_bound": "unbounded_above" if unbounded else "",
        "ci_low": shortest(record["fet.ci_95_lower"]),
        "ci_high": (
            "" if record["fet.ci_95_upper"] == "Infinity" else shortest(record["fet.ci_95_upper"])
        ),
        "pvalue": shortest(record["fet.p_value"]),
        "pvalue_test": "fisher_exact",
        # Supplementary Data 3 publishes `fet.p_value` and no corrected column at
        # all -- which is why the gene page names the study's 138,609 comparisons.
        "pvalue_adjusted": "",
        "pvalue_adjustment": "",
        "case_cohorts": CASE_COHORTS,
        "control_cohorts": CONTROL_COHORTS,
        "method_note": "",
        "source": SOURCE,
    }


def _sort_key(burden: dict[str, str]) -> tuple[tuple[int, object], ...]:
    """The table's sort key with ``maf_max`` compared as a number, empty smallest.

    ``validate_sort_order`` compares the values polars parsed, so ``maf_max`` is a
    float there; as text, "0.001" would sort after "0.0001".
    """
    out: list[tuple[int, object]] = []
    for column in _SORT_KEY:
        value = burden[column]
        if column == "maf_max":
            out.append((1, float(value)) if value else (0, 0.0))
        else:
            out.append((1, value) if value else (0, ""))
    return tuple(out)


def build(
    workbook: Path,
    *,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Write the ``audain2026_sd3`` package into the empty directory ``out`` and return it.

    ``revision`` has no default: a published version is permanent (F13), so the
    operator states whether this is the first build of the input (1) or a recipe
    change or a rebuild whose ``manifest.json`` bytes differ (2, 3, ...).
    ``built_by["recipe_commit"]`` is set by the build from ``commit_of``.
    """
    built_by, _ = stamp(built_by, commit_of)
    sha256 = hashlib.sha256(workbook.read_bytes()).hexdigest()
    if sha256 != PINNED_SHA256:
        raise PackageRefused(
            f"{workbook.name} has sha256 {sha256}, not the pinned {PINNED_SHA256}; "
            "this recipe was written against the pinned workbook only"
        )
    source = load_source(HERE / "source.yaml")

    rows: list[dict[str, str]] = []
    read = 0
    seen: set[str] = set()
    corrupted: list[str] = []
    for sheet, (stratum, carriers, denominator) in SHEETS.items():
        for record in read_records(workbook, sheet):
            read += 1
            symbol = record["SYMBOL"]
            seen.add(symbol)
            if _DATE_CORRUPTED.fullmatch(symbol):
                corrupted.append(symbol)
                continue
            rows.append(row(record, stratum, carriers, denominator))

    dropped: list[dict[str, Any]] = []
    if corrupted:
        dropped.append(
            {
                "reason": "Excel date corruption",
                "count": len(corrupted),
                "examples": sorted(set(corrupted))[:_EXAMPLES],
            }
        )

    ordered = sorted(rows, key=_sort_key)
    cells: list[Sequence[str | None]] = [[r[column] for column in ORDER] for r in ordered]
    return write_package(
        out,
        package=SOURCE,
        version=package_version(upstream_release=None, input_sha256=sha256, revision=revision),
        ids="symbol",
        coverage="every_gene",
        tables=[TableOut("burden", 1, ORDER, cells, "gene")],
        inputs=[
            {
                "name": WORKBOOK_NAME,
                "url": source["url"],
                "sha256": sha256,
                "retrieved_on": RETRIEVED_ON.isoformat(),
                "upstream_release": None,
            }
        ],
        input_counts={"rows": read, "genes": len(seen)},
        dropped=dropped,
        built_by=built_by,
        source={
            "name": source["name"],
            "version": f"sha256:{sha256}",
            "retrieved_on": RETRIEVED_ON.isoformat(),
            "url": source["url"],
            "licence": source["licence"],
            "redistribution": source["redistribution"],
        },
        licence_evidence=source["licence_evidence"],
    )
