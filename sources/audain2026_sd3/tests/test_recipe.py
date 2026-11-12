"""The Audain 2026 SD3 recipe on a synthetic three-sheet workbook: no network, no mirror.

The recipe refuses any workbook whose sha256 is not the pinned one, so the build
tests patch ``recipe.PINNED_SHA256`` to the synthetic file's hash;
``test_a_workbook_that_is_not_the_pinned_one_is_refused`` runs unpatched.

Reading the real workbook is out of scope here (no test reads ``local/``); its
counts are in the README, measured by a real build.
"""

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from feeds.package import PackageRefused
from feeds.version import package_version
from sources.audain2026_sd3 import recipe
from tests.xlsx_fixture import write_workbook_sheets

BUILT_BY = {"host": "login01"}
RECIPE_COMMIT = "a" * 40
SHEETS = ("aCHD_vs_controls", "sCHD_vs_controls", "nsCHD_vs_controls")


def _record(symbol: str, csq: str = "hcLOF", **over: str) -> dict[str, str]:
    base = {
        "SYMBOL": symbol,
        "analysis": "all_cases",
        "csq_group": csq,
        "maf": "1E-3",
        "n_het_cases": "30",
        "n_total_cases": "3876",
        "n_het_syndromic": "11",
        "n_total_syndromic": "1471",
        "n_het_nonsyndromic": "19",
        "n_total_nonsyndromic": "2405",
        "n_het_controls": "40",
        "n_total_controls": "45082",
        "fet.odds_ratio": "6.1099999999999999",
        "fet.ci_95_lower": "2.5",
        "fet.ci_95_upper": "9.75",
        "fet.p_value": "0.00012",
    }
    return base | over


def _table(records: list[dict[str, str]]) -> list[list[str]]:
    header = list(_record("X"))
    return [header, *[[r.get(c, "") for c in header] for r in records]]


ANALYSIS = {
    "aCHD_vs_controls": "all_cases",
    "sCHD_vs_controls": "syndromic",
    "nsCHD_vs_controls": "nonsyndromic",
}
# Each sheet's own statistics, so a build that reads one sheet under another's
# name publishes a different p-value and odds ratio.
STATS = {
    "aCHD_vs_controls": {"fet.p_value": "0.01", "fet.odds_ratio": "2.5"},
    "sCHD_vs_controls": {"fet.p_value": "0.002", "fet.odds_ratio": "3.5"},
    "nsCHD_vs_controls": {"fet.p_value": "0.3", "fet.odds_ratio": "1.25"},
}


def _workbook(path: Path, records: list[dict[str, str]]) -> Path:
    sheets = {"Info": [["Supplementary Data 3"]]}
    sheets |= {
        name: _table([r | {"analysis": ANALYSIS[name]} | STATS[name] for r in records])
        for name in SHEETS
    }
    return write_workbook_sheets(path, sheets)


@pytest.fixture
def pinned(monkeypatch: pytest.MonkeyPatch):  # type: ignore[no-untyped-def]
    def make(path: Path, records: list[dict[str, str]]) -> Path:
        _workbook(path, records)
        monkeypatch.setattr(recipe, "PINNED_SHA256", hashlib.sha256(path.read_bytes()).hexdigest())
        return path

    return make


def _build(workbook: Path, tmp_path: Path, **over: Any) -> Path:
    kwargs: dict[str, Any] = {
        "revision": 1,
        "built_by": BUILT_BY,
        "out": tmp_path / "pkg",
        "commit_of": lambda: RECIPE_COMMIT,
    }
    return recipe.build(workbook, **(kwargs | over))


def _read(package: Path, table: str) -> list[dict[str, str]]:
    lines = gzip.decompress((package / f"{table}.tsv.gz").read_bytes()).decode().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _manifest(package: Path) -> dict[str, Any]:
    loaded: dict[str, Any] = json.loads((package / "manifest.json").read_text())
    return loaded


def test_each_sheet_is_its_stratum_and_carries_its_own_carrier_columns() -> None:
    """Driven through the row builder: the stratum -> carrier/denominator pair is the contract."""
    record = _record("TBX5")
    got = {
        sheet: recipe.row(record, stratum, carriers, denominator)
        for sheet, (stratum, carriers, denominator) in recipe.SHEETS.items()
    }
    assert {
        s: (r["cohort_stratum"], r["n_case_carriers"], r["n_cases"]) for s, r in got.items()
    } == {
        "aCHD_vs_controls": ("all", "30", "3876"),
        "sCHD_vs_controls": ("syndromic", "11", "1471"),
        "nsCHD_vs_controls": ("nonsyndromic", "19", "2405"),
    }
    row = got["sCHD_vs_controls"]
    assert (row["n_control_carriers"], row["n_controls"]) == ("40", "45082")
    assert (row["study"], row["source"], row["gene"]) == ("PMID:42230622", "audain2026_sd3", "TBX5")
    assert (row["comparator"], row["count_unit"], row["origin"]) == (
        "control_cohort",
        "individuals",
        "any",
    )
    assert (row["case_cohorts"], row["control_cohorts"]) == ("cnchd;ddd;nottingham", "ukbb")
    assert (row["maf_max"], row["effect"], row["pvalue"]) == ("0.001", "6.11", "0.00012")
    assert (row["ci_low"], row["ci_high"], row["pvalue_test"]) == ("2.5", "9.75", "fisher_exact")
    assert tuple(row) == recipe.ORDER


def test_each_stratum_carries_its_own_sheets_statistics(tmp_path: Path, pinned: Any) -> None:
    """The sheet read for a stratum is the sheet it is named for (I1 of the final review)."""
    out = _build(pinned(tmp_path / "w.xlsx", [_record("TBX5")]), tmp_path)
    got = {r["cohort_stratum"]: (r["pvalue"], r["effect"]) for r in _read(out, "burden")}
    assert got == {
        "all": ("0.01", "2.5"),
        "syndromic": ("0.002", "3.5"),
        "nonsyndromic": ("0.3", "1.25"),
    }


def test_a_sheet_whose_analysis_column_names_another_comparison_is_refused(tmp_path: Path) -> None:
    sheets = {name: _table([_record("TBX5") | {"analysis": ANALYSIS[name]}]) for name in SHEETS}
    sheets["sCHD_vs_controls"] = _table(
        [_record("TBX5", analysis="syndromic"), _record("ZIC3", analysis="all_cases")]
    )
    path = write_workbook_sheets(tmp_path / "w.xlsx", sheets)
    assert len(recipe.read_records(path, "aCHD_vs_controls")) == 1
    with pytest.raises(PackageRefused, match="sCHD_vs_controls.*row 3.*'all_cases'.*'syndromic'"):
        recipe.read_records(path, "sCHD_vs_controls")


def test_an_infinite_odds_ratio_is_unbounded_with_no_number() -> None:
    record = _record("TAB2", **{"fet.odds_ratio": "Infinity", "fet.ci_95_upper": "Infinity"})
    row = recipe.row(record, "all", "n_het_cases", "n_total_cases")
    assert (row["effect"], row["ci_high"]) == ("", "")
    assert row["effect_bound"] == "unbounded_above"
    assert row["effect_measure"] == "odds_ratio"
    assert row["ci_low"] == "2.5"
    finite = recipe.row(_record("TAB2"), "all", "n_het_cases", "n_total_cases")
    assert finite["effect_bound"] == ""


def test_an_excel_date_symbol_is_dropped_and_named(tmp_path: Path, pinned: Any) -> None:
    """Dropped rows are counted in ``rows`` and listed in ``dropped``; nothing is guessed."""
    records = [
        _record("TBX5"),
        _record("2001-03-01"),
        _record("2001-03-01", "syn"),
        _record("2002-09-02"),
        _record("MARCH1"),
    ]
    out = _build(pinned(tmp_path / "w.xlsx", records), tmp_path)
    assert {r["gene"] for r in _read(out, "burden")} == {"TBX5", "MARCH1"}
    manifest = _manifest(out)
    # Three sheets x three date rows.
    assert manifest["dropped"] == [
        {
            "reason": "Excel date corruption",
            "count": 9,
            "examples": ["2001-03-01", "2002-09-02"],
        }
    ]
    assert manifest["input_counts"] == {"rows": 15, "genes": 4}
    assert [(t["name"], t["rows"], t["genes"]) for t in manifest["tables"]] == [("burden", 6, 2)]


def test_every_symbol_is_carried_as_the_workbook_spells_it(tmp_path: Path, pinned: Any) -> None:
    """TDGF1, BVES and C1orf127 stay as written; the atlas resolves them at ingest (F22)."""
    records = [_record(s) for s in ("TDGF1", "BVES", "C1orf127", "ACTC1")]
    out = _build(pinned(tmp_path / "w.xlsx", records), tmp_path)
    rows = _read(out, "burden")
    assert [r["gene"] for r in rows] == ["ACTC1"] * 3 + ["BVES"] * 3 + ["C1orf127"] * 3 + [
        "TDGF1"
    ] * 3
    manifest = _manifest(out)
    assert (manifest["ids"], manifest["coverage"], manifest["package"]) == (
        "symbol",
        "every_gene",
        "audain2026_sd3",
    )
    assert manifest["dropped"] == []
    (entry,) = manifest["inputs"]
    assert entry["name"] == "41525_2026_582_MOESM4_ESM.xlsx"
    assert (entry["retrieved_on"], entry["upstream_release"]) == ("2026-08-04", None)
    assert entry["url"].endswith("MediaObjects/41525_2026_582_MOESM4_ESM.xlsx")
    assert manifest["source"]["version"] == f"sha256:{recipe.PINNED_SHA256}"
    assert manifest["source"]["licence"] == "CC-BY-4.0"
    assert manifest["version"] == package_version(
        upstream_release=None, input_sha256=recipe.PINNED_SHA256
    )


def test_rows_are_sorted_with_maf_compared_as_a_number(tmp_path: Path, pinned: Any) -> None:
    """Literals, not a re-sort: gene, stratum, consequence, then maf numerically."""
    records = [
        _record("ZIC3", "syn"),
        _record("ZIC3", "hcLOF", maf="1E-4"),
        _record("ZIC3", "hcLOF", maf="1E-3"),
        _record("ACTC1", "missC"),
    ]
    out = _build(pinned(tmp_path / "w.xlsx", records), tmp_path)
    got = [
        (r["gene"], r["cohort_stratum"], r["consequence_class"], r["maf_max"])
        for r in _read(out, "burden")
        if r["cohort_stratum"] == "all"
    ]
    assert got == [
        ("ACTC1", "all", "missense_damaging", "0.001"),
        ("ZIC3", "all", "lof", "0.0001"),
        ("ZIC3", "all", "lof", "0.001"),
        ("ZIC3", "all", "synonymous", "0.001"),
    ]


def test_a_missing_sheet_or_column_or_vocabulary_is_refused(tmp_path: Path, pinned: Any) -> None:
    bad = _record("TBX5", "unheard_of")
    with pytest.raises(PackageRefused, match="csq_group.*unheard_of"):
        _build(pinned(tmp_path / "a.xlsx", [bad]), tmp_path)

    sheets = {name: _table([_record("TBX5")]) for name in SHEETS[:2]}
    path = write_workbook_sheets(tmp_path / "b.xlsx", sheets)
    with pytest.raises(PackageRefused, match="nsCHD_vs_controls"):
        recipe.read_records(path, "nsCHD_vs_controls")

    short = [row[:-1] for row in _table([_record("TBX5")])]
    path = write_workbook_sheets(tmp_path / "c.xlsx", {n: short for n in SHEETS})
    with pytest.raises(PackageRefused, match="fet.p_value"):
        recipe.read_records(path, "aCHD_vs_controls")


def test_a_workbook_that_is_not_the_pinned_one_is_refused(tmp_path: Path) -> None:
    workbook = _workbook(tmp_path / "w.xlsx", [_record("TBX5")])
    with pytest.raises(PackageRefused, match="not the pinned") as raised:
        _build(workbook, tmp_path)
    assert recipe.PINNED_SHA256 in str(raised.value)
    assert not (tmp_path / "pkg").exists()


def test_the_pinned_workbook_gives_version_8ebb6fdc() -> None:
    assert recipe.PINNED_SHA256 == (
        "8ebb6fdceb38b0099d6229a1ebf4d77db9f9492a62cfafb569d2865f9645ee2b"
    )
    assert package_version(upstream_release=None, input_sha256=recipe.PINNED_SHA256) == "8ebb6fdc"
