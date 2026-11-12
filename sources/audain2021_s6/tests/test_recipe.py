"""The Audain 2021 recipe on synthetic S6 and S8 workbooks: no network, no mirror.

Ported from the atlas's ``tests/unit/test_convert_audain2021.py`` (at 4b1f8f6),
minus the one of its four tests that reads the committed mirror. The converter's tests drove
``convert`` with the sheet readers monkeypatched; these write real workbooks and
call ``recipe.build``, so the reading path is exercised too. The recipe refuses
any workbook whose sha256 is not the pinned one, so the build tests patch the
two pins to the synthetic files' hashes; ``test_a_workbook_that_is_not_the_pinned_one_is_refused``
runs unpatched.

No gene or number here is a fact about a real gene. Reading the real workbooks
is out of scope (no test reads ``local/``); their counts are in the README.
"""

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from feeds.package import PackageRefused
from feeds.version import package_version
from sources.audain2021_s6 import recipe
from tests.xlsx_fixture import write_workbook, write_workbook_sheets

BUILT_BY = {"host": "login01"}
RECIPE_COMMIT = "a" * 40

# Written out rather than imported, so a change to `METHOD_NOTE` has to be made
# here too, deliberately. It names no case count of its own: the recipe's
# comment on `METHOD_NOTE` records why.
_METHOD_NOTE = (
    "Cases combine congenital heart disease with a series of 777 sporadic "
    "thoracic aortic aneurysm probands among the study's 7,958 cases; the study "
    "does not report how many of the 777 are in this row's case set. The study "
    "reports the CHD-only contribution separately for significant genes. A "
    "carrier is an individual with a rare deletion (MAF below 0.01) overlapping "
    "the gene."
)

# S8 Table's two header rows as the workbook has them, written out rather than
# built from the recipe's constants, then one gene per line: symbol, chr,
# start, end, and each band's case, control, p, corrected p. Note the trailing
# space in "significant CHD ", as the sheet spells it.
_S8_GROUPS = ["", "", "", "", "CHD+TAA", "", "", "", "only CHD", "", "", "", "only TAA"]
_S8_HEADER = [
    "symbol", "chr", "start", "end",
    "case", "control", "p_value", "p_value_adj",
    "case", "control", "p_value", "p_value_adj",
    "case", "control", "p_value", "p_value_adj",
    "significant CHD ", "significant TAA", "LOEUF", "DNM", "*Expression RNA-Seq", "Info",
]  # fmt: skip
# Synthetic genes, synthetic numbers. SPLIT's combined `case` is its S6 `case`
# (120), split 117 CHD-only and 3 TAA-only, and no two of its cells are alike.
_SPLIT_COMBINED = ["120", "21", "1E-4", "2E-4"]
_SPLIT_CHD_ONLY = ["117", "17", "3E-4", "4E-4"]
_SPLIT_TAA_ONLY = ["3", "13", "5E-4", "6E-4"]
_S8_SPLIT = ["SPLIT", "22", "11", "12", *_SPLIT_COMBINED, *_SPLIT_CHD_ONLY, *_SPLIT_TAA_ONLY]
# S8 leaves a band blank where it reports nothing for it.
_S8_BLANK = ["BLANK", "8", "31", "32", "29", "25", "7E-4", "8E-4", "", "", "", "", "29"]

_S6_HEADER = ["symbol", "chr", "case", "control", "p_value", "p_value_adj"]


def _s6_record(symbol: str, case: str = "120") -> list[str]:
    return [symbol, "1", case, "21", "1E-4", "2E-4"]


# SPLIT and BLANK, which S8 lists; UNLISTED, which it does not.
_S6_ROWS = [
    _s6_record("SPLIT"),
    _s6_record("BLANK", "29"),
    _s6_record("UNLISTED", "30"),
    _s6_record("UNREGISTERED"),
]


def _read(package: Path, table: str = "burden") -> list[dict[str, str]]:
    lines = gzip.decompress((package / f"{table}.tsv.gz").read_bytes()).decode().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _manifest(package: Path) -> dict[str, Any]:
    loaded: dict[str, Any] = json.loads((package / "manifest.json").read_text())
    return loaded


@pytest.fixture
def make(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):  # type: ignore[no-untyped-def]
    """Writes both workbooks, pins their hashes, builds, and returns the package path."""

    def run(s8: list[list[str]], s6_rows: list[list[str]] | None = None, **over: Any) -> Path:
        s6 = write_workbook(
            tmp_path / "s011.xlsx",
            [_S6_HEADER, *(_S6_ROWS if s6_rows is None else s6_rows)],
            "Sheet A",
        )
        s8_path = write_workbook(tmp_path / "s013.xlsx", s8, "S8 Table")
        monkeypatch.setattr(recipe, "PINNED_S6_SHA256", hashlib.sha256(s6.read_bytes()).hexdigest())
        monkeypatch.setattr(
            recipe, "PINNED_S8_SHA256", hashlib.sha256(s8_path.read_bytes()).hexdigest()
        )
        kwargs: dict[str, Any] = {
            "revision": 1,
            "built_by": BUILT_BY,
            "out": tmp_path / "pkg",
            "commit_of": lambda: RECIPE_COMMIT,
        }
        return recipe.build(s6, s8_path, **(kwargs | over))

    return run


def test_the_row_the_recipe_builds_carries_the_current_note_and_denominators(make: Any) -> None:
    """Driven through ``build`` on synthetic S6 and S8 workbooks.

    ``n_chd_case_carriers`` is S8's "only CHD" ``case`` (117 for SPLIT, not the
    combined 120 or the TAA-only 3), empty where S8 leaves that cell blank
    (BLANK) and where S8 does not list the gene (UNLISTED). Every row carries
    exactly the schema's columns, in the schema's order, and every S6 symbol
    gets a row, the unregistered one included.
    """
    package = make([_S8_GROUPS, _S8_HEADER, _S8_SPLIT, _S8_BLANK])
    rows = _read(package)

    assert [(r["gene"], r["n_chd_case_carriers"]) for r in rows] == [
        ("BLANK", ""),
        ("SPLIT", "117"),
        ("UNLISTED", ""),
        ("UNREGISTERED", ""),
    ]
    assert list(rows[0]) == list(recipe.ORDER)
    row = rows[1]
    assert row["method_note"] == _METHOD_NOTE
    # Every case source registered for this study, as a literal: the four
    # registered with the study, then the seven CHD series S2 identifies by
    # PMID with 50 or more cases (owner decision R12, 2026-10-26).
    assert row["case_cohorts"].split(";") == [
        "decipher",
        "isca",
        "pmid18535492_22349727_cases",
        "pmid22318994_cases",
        "pmid22912587_cases",
        "pmid22929023_cases",
        "pmid22939634_cases",
        "pmid22969434_cases",
        "pmid26655555_cases",
        "signature",
        "taa_cases",
    ]
    assert row["control_cohorts"] == "gain_controls;ottawa_controls;wtccc_controls"
    # Literals, not `N_CASES`/`N_CONTROLS`: a test reading the constants
    # compares the recipe to itself.
    assert (row["n_cases"], row["n_controls"]) == ("6746", "14024")
    assert (row["study"], row["source"]) == ("PMID:34324492", "audain2021_s6")
    assert (row["n_case_carriers"], row["n_control_carriers"]) == ("120", "21")
    assert (row["pvalue"], row["pvalue_adjusted"]) == ("0.0001", "0.0002")
    assert (row["pvalue_test"], row["pvalue_adjustment"]) == (
        "permutation",
        "familywise_permutation",
    )
    assert (row["variant_class"], row["consequence_class"], row["maf_max"]) == (
        "cnv_deletion",
        "lof",
        "0.01",
    )
    assert (row["cohort_stratum"], row["count_unit"], row["comparator"], row["origin"]) == (
        "all",
        "individuals",
        "control_cohort",
        "any",
    )
    assert (row["effect"], row["effect_measure"], row["ci_low"]) == ("", "", "")


def test_the_package_records_both_inputs_every_gene_and_the_symbol_ids(make: Any) -> None:
    package = make([_S8_GROUPS, _S8_HEADER, _S8_SPLIT, _S8_BLANK])
    manifest = _manifest(package)
    assert (manifest["package"], manifest["ids"], manifest["coverage"]) == (
        "audain2021_s6",
        "symbol",
        "every_gene",
    )
    assert manifest["input_counts"] == {"rows": 4, "genes": 4}
    assert manifest["dropped"] == []
    s6, s8 = manifest["inputs"]
    assert (s6["name"], s6["sha256"]) == ("s011.xlsx", recipe.PINNED_S6_SHA256)
    assert (s8["name"], s8["sha256"]) == ("s013.xlsx", recipe.PINNED_S8_SHA256)
    assert s6["url"].endswith("pgen.1009679.s011&type=supplementary")
    assert s8["url"].endswith("pgen.1009679.s013&type=supplementary")
    assert (s6["retrieved_on"], s8["retrieved_on"]) == ("2026-08-05", "2026-08-05")
    assert manifest["source"]["version"] == f"sha256:{recipe.PINNED_S6_SHA256}"
    assert (manifest["source"]["licence"], manifest["source"]["redistribution"]) == (
        "CC0-1.0",
        "permitted",
    )
    assert manifest["version"] == package_version(
        upstream_release=None, input_sha256=recipe.PINNED_S6_SHA256
    )


def test_the_burden_columns_are_the_vendored_schemas() -> None:
    schema = json.loads(
        (Path(__file__).resolve().parents[3] / "schemas/atlas/tables/burden.v1.json").read_text()
    )
    assert [field["name"] for field in schema["fields"]] == list(recipe.ORDER)


@pytest.mark.parametrize(
    ("groups", "header", "split"),
    [
        # The TAA-only band before the CHD-only one: an index fixed at the
        # workbook's column 8 would publish the TAA-only 3.
        pytest.param(
            [*_S8_GROUPS[:8], "only TAA", "", "", "", "only CHD"],
            _S8_HEADER,
            ["SPLIT", "22", "11", "12", *_SPLIT_COMBINED, *_SPLIT_TAA_ONLY, *_SPLIT_CHD_ONLY],
            id="bands-reordered",
        ),
        # `case` second in every band: a column taken from the band's label
        # alone, without reading the second header row, is each band's
        # `control`.
        pytest.param(
            _S8_GROUPS,
            ["symbol", "chr", "start", "end", *("control", "case", "p_value", "p_value_adj") * 3],
            [
                "SPLIT", "22", "11", "12",
                "21", "120", "1E-4", "2E-4",
                "17", "117", "3E-4", "4E-4",
                "13", "3", "5E-4", "6E-4",
            ],
            id="case-second-in-each-band",
        ),
    ],
)  # fmt: skip
def test_the_chd_only_count_is_read_from_the_column_its_headers_name(
    make: Any, groups: list[str], header: list[str], split: list[str]
) -> None:
    """Located from S8's header rows, so a moved column is still the right one.

    BUR003 cannot catch a misread band: every band's ``case`` is at most the
    row's case carriers (see the recipe's CHD-ONLY CARRIERS), so the wrong one
    validates clean. The combined band is located the same way and still
    checked against S6, which these layouts pass only if it is found too.
    """
    package = make([groups, header, split], s6_rows=[_s6_record("SPLIT")])
    assert [r["n_chd_case_carriers"] for r in _read(package)] == ["117"]


@pytest.mark.parametrize(
    ("s8", "match"),
    [
        pytest.param(
            [[*_S8_GROUPS[:8], "CHD only", "", "", "", "only TAA"], _S8_HEADER, _S8_SPLIT],
            "first header row reads 'only CHD' 0 times",
            id="band-unlabelled",
        ),
        pytest.param(
            [[*_S8_GROUPS[:8], "only CHD", "", "", "", "only CHD"], _S8_HEADER, _S8_SPLIT],
            "first header row reads 'only CHD' 2 times",
            id="band-labelled-twice",
        ),
        # The band's columns renamed: nothing under it reads `case`.
        pytest.param(
            [_S8_GROUPS, [*_S8_HEADER[:8], "cases", *_S8_HEADER[9:]], _S8_SPLIT],
            "second header row under 'only CHD' reads 'case' 0 times",
            id="case-unnamed",
        ),
        # With the next band unlabelled, the CHD-only band runs to the end of
        # the row and holds two `case` columns, and neither can be chosen.
        pytest.param(
            [_S8_GROUPS[:12], _S8_HEADER, _S8_SPLIT],
            "second header row under 'only CHD' reads 'case' 2 times",
            id="next-band-unlabelled",
        ),
        # S8's combined count is not S6's: the two tables are not one analysis.
        pytest.param(
            [_S8_GROUPS, _S8_HEADER, [*_S8_SPLIT[:4], "121", *_S8_SPLIT[5:]]],
            "differs from S6",
            id="combined-is-not-s6",
        ),
        # A gene S6 does not list cannot be checked against it.
        pytest.param(
            [_S8_GROUPS, _S8_HEADER, ["ELSEWHERE", *_S8_SPLIT[1:]]],
            "S6 lists no such gene, for 1 gene\\(s\\): ELSEWHERE",
            id="gene-not-in-s6",
        ),
        pytest.param(
            [_S8_GROUPS, _S8_HEADER, _S8_SPLIT, _S8_SPLIT],
            "S8 Table lists SPLIT twice",
            id="gene-listed-twice",
        ),
    ],
)
def test_the_recipe_stops_on_an_s8_table_it_cannot_read_as_s6s_split(
    make: Any, s8: list[list[str]], match: str
) -> None:
    """Refusals rather than reports: each would otherwise write a wrong count, or none."""
    with pytest.raises(PackageRefused, match=match):
        make(s8, s6_rows=[_s6_record("SPLIT")])


def test_a_gene_s6_lists_twice_is_refused(make: Any) -> None:
    with pytest.raises(PackageRefused, match="S6 Table lists SPLIT more than once"):
        make([_S8_GROUPS, _S8_HEADER], s6_rows=[_s6_record("SPLIT"), _s6_record("SPLIT")])


def test_an_s6_sheet_without_a_needed_column_is_refused(make: Any, tmp_path: Path) -> None:
    header = [c for c in _S6_HEADER if c != "p_value_adj"]
    path = write_workbook(
        tmp_path / "x.xlsx", [header, ["SPLIT", "1", "120", "21", "1E-4"]], "Sheet A"
    )
    with pytest.raises(PackageRefused, match="no column p_value_adj"):
        recipe.read_s6(path)


@pytest.mark.parametrize("wrong", ["s011.xlsx", "s013.xlsx"])
def test_a_workbook_that_is_not_the_pinned_one_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, wrong: str
) -> None:
    """Each input is checked: the other is pinned to its synthetic hash, this one is not."""
    s6 = write_workbook_sheets(tmp_path / "s011.xlsx", {"Sheet A": [_S6_HEADER, *_S6_ROWS]})
    s8 = write_workbook_sheets(tmp_path / "s013.xlsx", {"S8 Table": [_S8_GROUPS, _S8_HEADER]})
    if wrong != "s011.xlsx":
        monkeypatch.setattr(recipe, "PINNED_S6_SHA256", hashlib.sha256(s6.read_bytes()).hexdigest())
    if wrong != "s013.xlsx":
        monkeypatch.setattr(recipe, "PINNED_S8_SHA256", hashlib.sha256(s8.read_bytes()).hexdigest())
    with pytest.raises(PackageRefused, match="not the pinned") as raised:
        recipe.build(
            s6,
            s8,
            revision=1,
            built_by=BUILT_BY,
            out=tmp_path / "pkg",
            commit_of=lambda: RECIPE_COMMIT,
        )
    assert str(raised.value).startswith(wrong)
    assert not (tmp_path / "pkg").exists()


def test_the_pins_are_the_atlas_sources_and_give_version_e8bd50c0() -> None:
    assert recipe.PINNED_S6_SHA256 == (
        "e8bd50c0a454d3f5563a6c1ee52c62206d98840cc696c3c8db4c0e0eccb63e62"
    )
    assert recipe.PINNED_S8_SHA256 == (
        "d070834994949311b3cff7a8a8470944c6ddd14c4644351062b0243c57db0ebc"
    )
    assert (
        package_version(upstream_release=None, input_sha256=recipe.PINNED_S6_SHA256) == "e8bd50c0"
    )
