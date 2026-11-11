from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from feeds.package import TableOut

BURDEN_COLUMNS = (
    "study gene cohort_stratum lesion_group variant_class consequence_class origin maf_max "
    "count_unit n_case_carriers n_cases n_chd_case_carriers comparator n_control_carriers "
    "n_controls expected_count effect effect_measure effect_bound ci_low ci_high pvalue "
    "pvalue_test pvalue_adjusted pvalue_adjustment case_cohorts control_cohorts method_note source"
).split()


def burden_row(gene: str, **over: str | None) -> list[str | None]:
    base: dict[str, str | None] = {c: None for c in BURDEN_COLUMNS}
    base |= {
        "study": "PMID:1",
        "gene": gene,
        "cohort_stratum": "all",
        "variant_class": "snv_indel",
        "consequence_class": "lof",
        "origin": "any",
        "count_unit": "individuals",
        "n_case_carriers": "3",
        "n_cases": "100",
        "comparator": "none",
        "case_cohorts": "pcgc",
        "source": "demo",
    }
    base |= over
    return [base[c] for c in BURDEN_COLUMNS]


def burden_table(rows: list[list[str | None]] | None = None) -> TableOut:
    rows = rows if rows is not None else [burden_row("TBX5"), burden_row("GATA4")]
    return TableOut("burden", 1, tuple(BURDEN_COLUMNS), rows, "gene")


def package_kwargs(tables: list[TableOut] | None = None, ids: str = "symbol") -> dict[str, Any]:
    return {
        "package": "demo",
        "version": "2026-11-03",
        "ids": ids,
        "coverage": "every_gene",
        "tables": tables if tables is not None else [burden_table()],
        "inputs": [
            {
                "name": "in.tsv",
                "url": "https://example.org/in.tsv",
                "sha256": "a" * 64,
                "retrieved_on": "2026-11-03",
                "upstream_release": None,
            }
        ],
        "input_counts": {"rows": 2, "genes": 2},
        "dropped": [],
        "built_by": {"recipe_commit": "b" * 40, "host": "test"},
        "source": {
            "name": "Demo",
            "version": "2026-11-03",
            "retrieved_on": "2026-11-03",
            "url": "https://example.org",
            "licence": "CC-BY-4.0",
            "redistribution": "permitted_with_attribution",
        },
        "licence_evidence": "stated on the landing page",
    }


@pytest.fixture
def out(tmp_path: Path) -> Path:
    return tmp_path / "pkg"
