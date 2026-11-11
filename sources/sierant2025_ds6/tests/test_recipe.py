"""The Sierant 2025 recipe on a synthetic workbook: no network, no committed mirror.

The workbook writer is the one the atlas's ``tests/unit/test_convert_sierant2025.py``
uses, ported. The recipe refuses any workbook whose sha256 is not the pinned one,
so the build tests patch ``recipe.PINNED_SHA256`` to the synthetic file's hash;
``test_a_workbook_that_is_not_the_pinned_one_is_refused`` runs unpatched.

The guards each have a test that fails without them (the mutations were applied
by hand from an in-memory copy; the results are in the commit message):
``_check_layout``, ``_check_bands``, the Excel-date drop, the pinned-hash refusal
in ``build`` and in ``fetch``.
"""

import gzip
import hashlib
import html
import json
import zipfile
from pathlib import Path
from typing import Any

import pytest

from feeds.package import PackageRefused
from feeds.version import package_version
from sources.sierant2025_ds6 import fetch, recipe

BUILT_BY = {"recipe_commit": "a" * 40, "host": "login01"}
STUDY = "PMID:40127276"

# Row 7 of the real workbook, written out so a change to `_EXPECTED_HEADERS` has to
# be made twice deliberately instead of the guard agreeing with whatever it became.
REAL_HEADER = {
    2: "pLI",
    3: "MisZ",
    5: "LOF",
    6: "MisD",
    7: "Obs",
    8: "Exp",
    9: "Enrich",
    10: "P-Val",
    11: "FDR",
    13: "LOF",
    14: "MisD",
    15: "Obs",
    16: "Total",
    17: "AF",
    18: "LOF",
    19: "MisD",
    20: "Obs",
    21: "Total",
    22: "AF",
    23: "Enrich",
    24: "P-Val",
    25: "FDR",
    27: "P-Val",
    28: "FDR",
    30: "Obs",
    31: "Exp",
    32: "Enrich",
    33: "P-Val",
    34: "FDR",
    36: "Obs",
    37: "Total",
    38: "AF",
    39: "Obs",
    40: "Total",
    41: "AF",
    42: "Enrich",
    43: "P-Val",
    44: "FDR",
    46: "P-Val",
    47: "FDR",
    49: "Obs",
    50: "Exp",
    51: "Enrich",
    52: "P-Val",
    53: "FDR",
    55: "Obs",
    56: "Total",
    57: "AF",
    58: "Obs",
    59: "Total",
    60: "AF",
    61: "Enrich",
    62: "P-Val",
    63: "FDR",
    65: "P-Val",
    66: "FDR",
}
REAL_BANDS = {
    23: "Fisher's Exact",
    27: "JL-FDR",
    42: "Fisher's Exact",
    46: "JL-FDR",
    61: "Fisher's Exact",
    65: "JL-FDR",
}

# The workbook's strict OOXML namespace, which a reader assuming the transitional
# one reads as zero rows.
_STRICT = "http://purl.oclc.org/ooxml/spreadsheetml/main"


def _column_letters(index: int) -> str:
    letters = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(ord("A") + remainder) + letters
    return letters


def _workbook(path: Path, rows: list[list[str]]) -> Path:
    """The smallest xlsx ``read_sheet_rows`` accepts: ``rows`` as inline strings on SHEET."""
    data = "".join(
        f'<row r="{number}">'
        + "".join(
            f'<c r="{_column_letters(index)}{number}" t="inlineStr">'
            f"<is><t>{html.escape(value)}</t></is></c>"
            for index, value in enumerate(row)
            if value
        )
        + "</row>"
        for number, row in enumerate(rows, start=1)
    )
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(
            "xl/workbook.xml",
            f'<workbook xmlns="{_STRICT}" '
            'xmlns:r="http://purl.oclc.org/ooxml/officeDocument/relationships">'
            f'<sheets><sheet name="{html.escape(recipe.SHEET)}" sheetId="1" r:id="rId1"/></sheets>'
            "</workbook>",
        )
        archive.writestr(
            "xl/_rels/workbook.xml.rels",
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>',
        )
        archive.writestr(
            "xl/worksheets/sheet1.xml",
            f'<worksheet xmlns="{_STRICT}"><sheetData>{data}</sheetData></worksheet>',
        )
    return path


def _band_row() -> list[str]:
    row = [""] * 70
    for index, label in REAL_BANDS.items():
        row[index] = label
    return row


def _header_row() -> list[str]:
    row = [""] * 70
    for index, label in REAL_HEADER.items():
        row[index] = label
    return row


def _record(symbol: str, seed: int) -> list[str]:
    """A data row whose statistic cells are distinct per (gene, column), so a misread shows."""
    row = [f"0.{seed:02d}{index:03d}" for index in range(70)]
    row[recipe.GENE_COLUMN] = symbol
    for block in recipe._BLOCKS.values():
        # The count columns must be whole numbers the table schema accepts.
        row[block["dnm_obs"]] = "2"
        row[block["cc_obs"]], row[block["cc_total"]] = "3", "1000"
        row[block["cc_ctrl"]], row[block["cc_ctrl_total"]] = "1", "5000"
    return row


def _sheet(symbols: list[str], bands: list[str] | None = None, header: list[str] | None = None):  # type: ignore[no-untyped-def]
    records = [_record(symbol, n) for n, symbol in enumerate(symbols, start=1)]
    return [[] for _ in range(recipe._BAND_ROW)] + [
        bands or _band_row(),
        header or _header_row(),
        *records,
    ]


@pytest.fixture
def pinned(monkeypatch: pytest.MonkeyPatch):  # type: ignore[no-untyped-def]
    """Make a synthetic workbook the 'pinned' one: returns a writer that patches the hash."""

    def make(path: Path, rows: list[list[str]]) -> Path:
        _workbook(path, rows)
        monkeypatch.setattr(recipe, "PINNED_SHA256", hashlib.sha256(path.read_bytes()).hexdigest())
        return path

    return make


def _read(package: Path, table: str) -> list[dict[str, str]]:
    lines = gzip.decompress((package / f"{table}.tsv.gz").read_bytes()).decode().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _manifest(package: Path) -> dict[str, Any]:
    loaded: dict[str, Any] = json.loads((package / "manifest.json").read_text())
    return loaded


def test_every_panel_gene_is_carried_under_the_papers_symbol(tmp_path: Path, pinned: Any) -> None:
    """Three genes given in reverse order: six arm rows and three combined rows each, sorted."""
    workbook = pinned(tmp_path / "w.xlsx", _sheet(["ZIC3", "NOT-A-HGNC-SYMBOL", "ACTC1"]))
    out = recipe.build(workbook, built_by=BUILT_BY, out=tmp_path / "pkg")

    arms = _read(out, "burden")
    combined = _read(out, "burden_combined")
    # Total sort over the schema's key, literals: gene, then consequence, then origin.
    assert [(r["gene"], r["consequence_class"], r["origin"]) for r in arms[:6]] == [
        ("ACTC1", "damaging", "de_novo"),
        ("ACTC1", "damaging", "transmitted_or_unphased"),
        ("ACTC1", "lof", "de_novo"),
        ("ACTC1", "lof", "transmitted_or_unphased"),
        ("ACTC1", "missense_damaging", "de_novo"),
        ("ACTC1", "missense_damaging", "transmitted_or_unphased"),
    ]
    assert [r["gene"] for r in arms] == sorted(
        [g for g in ("ACTC1",) * 6 + ("NOT-A-HGNC-SYMBOL",) * 6 + ("ZIC3",) * 6]
    )
    assert [(r["gene"], r["consequence_class"]) for r in combined] == [
        (g, c)
        for g in ("ACTC1", "NOT-A-HGNC-SYMBOL", "ZIC3")
        for c in ("damaging", "lof", "missense_damaging")
    ]
    # The symbol is verbatim, never looked up: no HGNC id appears anywhere.
    assert not any(r["gene"].startswith("HGNC:") for r in arms + combined)
    assert {r["source"] for r in arms + combined} == {"sierant2025_ds6"}

    manifest = _manifest(out)
    assert (manifest["ids"], manifest["coverage"]) == ("symbol", "every_gene")
    assert manifest["package"] == "sierant2025_ds6"
    assert [(t["name"], t["rows"], t["genes"]) for t in manifest["tables"]] == [
        ("burden", 18, 3),
        ("burden_combined", 9, 3),
    ]
    assert manifest["dropped"] == []
    assert manifest["input_counts"] == {"rows": 3, "genes": 3}
    (entry,) = manifest["inputs"]
    assert entry["name"] == fetch.WORKBOOK_NAME
    assert entry["url"] == fetch.BUNDLE_URL
    assert entry["sha256"] == recipe.PINNED_SHA256
    assert (entry["retrieved_on"], entry["upstream_release"]) == ("2026-08-05", None)
    assert manifest["version"] == package_version(
        upstream_release=None, input_sha256=recipe.PINNED_SHA256
    )
    assert manifest["source"]["version"] == f"sha256:{recipe.PINNED_SHA256}"
    assert (manifest["source"]["retrieved_on"], manifest["source"]["licence"]) == (
        "2026-08-05",
        "CC-BY-4.0",
    )


def test_the_pinned_workbook_gives_version_d8e9f777() -> None:
    assert package_version(upstream_release=None, input_sha256=fetch.PINNED_SHA256) == "d8e9f777"
    assert fetch.PINNED_SHA256 == (
        "d8e9f7777b329c810adde8b9f4a0f91a4d1a6211b760daeb9329ae6f36606daa"
    )


def test_an_excel_date_symbol_is_dropped_and_named(tmp_path: Path, pinned: Any) -> None:
    """Both date shapes are dropped with their count and examples; nothing is guessed."""
    symbols = ["TBX5", "2001-03-01", "1-Mar", "12-Sep", "GATA4", "2002-09-02"]
    workbook = pinned(tmp_path / "w.xlsx", _sheet(symbols))
    out = recipe.build(workbook, built_by=BUILT_BY, out=tmp_path / "pkg")

    assert {r["gene"] for r in _read(out, "burden")} == {"GATA4", "TBX5"}
    assert {r["gene"] for r in _read(out, "burden_combined")} == {"GATA4", "TBX5"}
    manifest = _manifest(out)
    assert manifest["dropped"] == [
        {
            "reason": "Excel date corruption",
            "count": 4,
            "examples": ["1-Mar", "12-Sep", "2001-03-01", "2002-09-02"],
        }
    ]
    # `rows` counts the panel rows read, dropped ones included.
    assert manifest["input_counts"] == {"rows": 6, "genes": 6}
    assert [t["genes"] for t in manifest["tables"]] == [2, 2]


def test_a_gene_that_merely_looks_like_a_date_is_kept(tmp_path: Path, pinned: Any) -> None:
    """``MARC1``-style names and symbols with a digit-dash tail are real genes."""
    workbook = pinned(tmp_path / "w.xlsx", _sheet(["MARCH1", "HLA-B", "NKX2-5", "SEPT7"]))
    out = recipe.build(workbook, built_by=BUILT_BY, out=tmp_path / "pkg")
    assert {r["gene"] for r in _read(out, "burden")} == {"MARCH1", "HLA-B", "NKX2-5", "SEPT7"}
    assert _manifest(out)["dropped"] == []


def _synthetic(**cells: str) -> list[str]:
    record = [""] * 70
    record[recipe.GENE_COLUMN] = "KMT2D"
    for key, value in cells.items():
        record[recipe._BLOCKS["lof"][key]] = value
    return record


def test_the_arm_rows_keep_their_comparators_and_count_units() -> None:
    """Driven through the row builders (the labelling contract is decided there).

    The de novo arm counts mutations against 3,887 trios; the case-control arm
    counts alleles against a per-gene allele total.
    """
    record = _synthetic(
        dnm_obs="21",
        dnm_exp="0.16558",
        dnm_effect="126.82",
        dnm_p="6.65e-37",
        dnm_fdr="1.6491e-34",
        cc_obs="12",
        cc_total="21768",
        cc_ctrl="35",
        cc_ctrl_total="227965",
        cc_effect="3.5919",
        cc_p="5.4935e-4",
        cc_fdr="1.2385e-2",
    )
    ix = recipe._BLOCKS["lof"]

    de_novo = recipe._de_novo_row("KMT2D", "lof", record, ix)
    assert de_novo["count_unit"] == "de_novo_mutations"
    assert de_novo["origin"] == "de_novo"
    assert de_novo["comparator"] == "mutation_model"
    assert (de_novo["n_case_carriers"], de_novo["n_cases"]) == ("21", "3887")
    assert de_novo["expected_count"] == "0.16558"
    assert de_novo["effect_measure"] == "enrichment_ratio"
    assert de_novo["pvalue_test"] == "poisson"
    assert de_novo["n_control_carriers"] == de_novo["n_controls"] == ""
    assert de_novo["control_cohorts"] == ""
    assert de_novo["maf_max"] == ""

    case_control = recipe._case_control_row("KMT2D", "lof", record, ix)
    assert case_control["count_unit"] == "alleles"
    assert case_control["origin"] == "transmitted_or_unphased"
    assert case_control["comparator"] == "control_cohort"
    assert (case_control["n_case_carriers"], case_control["n_cases"]) == ("12", "21768")
    assert (case_control["n_control_carriers"], case_control["n_controls"]) == ("35", "227965")
    assert case_control["effect_measure"] == "odds_ratio"
    assert case_control["pvalue_test"] == "fisher_exact"
    assert case_control["expected_count"] == ""
    assert case_control["maf_max"] == "1e-05"

    for row in (de_novo, case_control):
        assert row["pvalue_adjustment"] == "benjamini_hochberg"
        assert row["case_cohorts"] == "pcgc"
        assert row["study"] == STUDY
        assert tuple(row) == recipe.ORDER


def test_an_infinite_odds_ratio_publishes_no_number_and_still_names_its_measure() -> None:
    """The workbook stores a literal ``Inf`` where no control allele carries."""
    record = _synthetic(
        cc_obs="1",
        cc_total="21901",
        cc_ctrl="0",
        cc_ctrl_total="270942",
        cc_effect="Inf",
        cc_p="0.0747875141287311",
        cc_fdr="0.3",
    )
    ix = recipe._BLOCKS["lof"]
    row = recipe._case_control_row("RPL5", "lof", record, ix)
    assert row["effect"] == ""
    assert row["effect_bound"] == "unbounded_above"
    assert row["effect_measure"] == "odds_ratio"
    assert row["ci_low"] == row["ci_high"] == ""

    record[ix["cc_effect"]] = "3.5919514463315401"
    finite = recipe._case_control_row("KMT2D", "lof", record, ix)
    assert finite["effect"] == "3.59195144633154"
    assert finite["effect_bound"] == ""


def test_the_combination_is_read_from_its_own_columns_never_from_fishers() -> None:
    record = _synthetic(
        cc_p="5.4935e-4", cc_fdr="1.2385e-2", meta_p="1.6582e-38", meta_fdr="2.6973e-34"
    )
    row = recipe._combined_row("KMT2D", "lof", record, recipe._BLOCKS["lof"])
    assert (row["pvalue"], row["local_fdr"]) == ("1.6582e-38", "2.6973e-34")
    assert row["combines"] == "mutation_model;control_cohort"
    assert row["method"] == "joint_local_fdr"
    assert (row["study"], row["gene"], row["consequence_class"]) == (STUDY, "KMT2D", "lof")
    assert tuple(row) == recipe.COMBINED_ORDER
    assert row["method_note"].endswith("The study's genome-wide threshold is p < 2.6e-6.")


def test_each_consequence_block_reads_its_own_band_of_the_sheet(
    tmp_path: Path, pinned: Any
) -> None:
    """Every cell in the synthetic record is distinct, so a misplaced index shows."""
    workbook = pinned(tmp_path / "w.xlsx", _sheet(["TBX5"]))
    out = recipe.build(workbook, built_by=BUILT_BY, out=tmp_path / "pkg")
    got = [
        (r["consequence_class"], r["pvalue"], r["local_fdr"]) for r in _read(out, "burden_combined")
    ]
    assert got == [
        ("damaging", "0.01027", "0.01028"),
        ("lof", "0.01046", "0.01047"),
        ("missense_damaging", "0.01065", "0.01066"),
    ]
    used = [index for block in recipe._BLOCKS.values() for index in block.values()]
    assert len(used) == len(set(used)) == 42, "two blocks read the same column"


def test_a_shifted_header_is_refused_rather_than_transcribed(tmp_path: Path, pinned: Any) -> None:
    header = _header_row()
    shifted = header[:7] + ["Obs"] + header[7:]
    workbook = pinned(tmp_path / "w.xlsx", _sheet(["TBX5"], header=shifted))
    with pytest.raises(PackageRefused, match="does not match the layout") as raised:
        recipe.build(workbook, built_by=BUILT_BY, out=tmp_path / "pkg")
    assert "column 8" in str(raised.value)
    assert not (tmp_path / "pkg").exists()


def test_a_workbook_whose_fisher_and_jl_fdr_bands_are_swapped_is_refused(
    tmp_path: Path, pinned: Any
) -> None:
    bands = _band_row()
    bands[23], bands[27] = bands[27], bands[23]
    workbook = pinned(tmp_path / "w.xlsx", _sheet(["TBX5"], bands=bands))
    with pytest.raises(PackageRefused, match="band row") as raised:
        recipe.build(workbook, built_by=BUILT_BY, out=tmp_path / "pkg")
    assert "column 23" in str(raised.value)


def test_a_header_row_too_short_to_hold_the_layout_is_refused_by_name() -> None:
    with pytest.raises(PackageRefused, match="past the end of the row") as raised:
        recipe._check_layout(["", ""])
    assert "column 7" in str(raised.value)


def test_a_workbook_that_is_not_the_pinned_one_is_refused(tmp_path: Path) -> None:
    """Unpatched: a perfectly good synthetic workbook has the wrong hash."""
    workbook = _workbook(tmp_path / "w.xlsx", _sheet(["TBX5"]))
    with pytest.raises(PackageRefused, match="not the pinned") as raised:
        recipe.build(workbook, built_by=BUILT_BY, out=tmp_path / "pkg")
    assert recipe.PINNED_SHA256 in str(raised.value)
    assert not (tmp_path / "pkg").exists()


def test_a_symbol_on_two_panel_rows_is_refused(tmp_path: Path, pinned: Any) -> None:
    workbook = pinned(tmp_path / "w.xlsx", _sheet(["TBX5", "GATA4", "TBX5"]))
    with pytest.raises(PackageRefused, match="more than one row.*TBX5"):
        recipe.build(workbook, built_by=BUILT_BY, out=tmp_path / "pkg")


# ---- fetch: a local file:// bundle stands in for Europe PMC ----------------------------


def _bundle(path: Path, members: dict[str, bytes]) -> str:
    with zipfile.ZipFile(path, "w") as archive:
        for name, data in members.items():
            archive.writestr(name, data)
    return path.as_uri()


def test_fetch_extracts_only_the_pinned_workbook(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"workbook bytes"
    monkeypatch.setattr(fetch, "PINNED_SHA256", hashlib.sha256(content).hexdigest())
    uri = _bundle(
        tmp_path / "b.zip",
        {f"supp/{fetch.WORKBOOK_NAME}": content, "supp/other.xlsx": b"x", "readme.txt": b"y"},
    )
    monkeypatch.setattr(fetch, "BUNDLE_URL", uri)
    dest = tmp_path / "dest"
    path = fetch.fetch(dest)
    assert path == dest / fetch.WORKBOOK_NAME
    assert path.read_bytes() == content
    assert sorted(p.name for p in dest.iterdir()) == [fetch.WORKBOOK_NAME]


def test_fetch_refuses_another_sha256_and_leaves_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    uri = _bundle(tmp_path / "b.zip", {fetch.WORKBOOK_NAME: b"not the workbook"})
    monkeypatch.setattr(fetch, "BUNDLE_URL", uri)
    dest = tmp_path / "dest"
    with pytest.raises(fetch.FetchRefused, match="not the pinned"):
        fetch.fetch(dest)
    assert list(dest.iterdir()) == []


@pytest.mark.parametrize("member", [f"../{fetch.WORKBOOK_NAME}", f"/etc/{fetch.WORKBOOK_NAME}"])
def test_fetch_refuses_a_member_path_that_escapes_the_destination(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, member: str
) -> None:
    content = b"workbook bytes"
    monkeypatch.setattr(fetch, "PINNED_SHA256", hashlib.sha256(content).hexdigest())
    monkeypatch.setattr(fetch, "BUNDLE_URL", _bundle(tmp_path / "b.zip", {member: content}))
    dest = tmp_path / "dest"
    with pytest.raises(fetch.FetchRefused, match="escapes"):
        fetch.fetch(dest)
    assert list(dest.iterdir()) == []
    assert not (tmp_path / fetch.WORKBOOK_NAME).exists()


def test_fetch_refuses_a_bundle_without_exactly_one_workbook(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(fetch, "BUNDLE_URL", _bundle(tmp_path / "b.zip", {"other.xlsx": b"x"}))
    with pytest.raises(fetch.FetchRefused, match="exactly one"):
        fetch.fetch(tmp_path / "dest")


def test_fetch_refuses_a_backslash_in_a_member_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``a\\b/<name>`` has a clean POSIX path; only the backslash clause refuses it."""
    content = b"workbook bytes"
    monkeypatch.setattr(fetch, "PINNED_SHA256", hashlib.sha256(content).hexdigest())
    uri = _bundle(tmp_path / "b.zip", {f"a\\..\\b/{fetch.WORKBOOK_NAME}": content})
    monkeypatch.setattr(fetch, "BUNDLE_URL", uri)
    with pytest.raises(fetch.FetchRefused, match="escapes"):
        fetch.fetch(tmp_path / "dest")


def test_fetch_refuses_an_existing_target(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    content = b"workbook bytes"
    monkeypatch.setattr(fetch, "PINNED_SHA256", hashlib.sha256(content).hexdigest())
    monkeypatch.setattr(
        fetch, "BUNDLE_URL", _bundle(tmp_path / "b.zip", {fetch.WORKBOOK_NAME: content})
    )
    dest = tmp_path / "dest"
    dest.mkdir()
    (dest / fetch.WORKBOOK_NAME).write_bytes(b"earlier")
    with pytest.raises(fetch.FetchRefused, match="exists"):
        fetch.fetch(dest)
    assert (dest / fetch.WORKBOOK_NAME).read_bytes() == b"earlier"


def test_fetch_refuses_a_response_that_is_not_a_zip(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    not_zip = tmp_path / "page.html"
    not_zip.write_bytes(b"<html>500</html>")
    monkeypatch.setattr(fetch, "BUNDLE_URL", not_zip.as_uri())
    with pytest.raises(fetch.FetchRefused, match="not return a zip"):
        fetch.fetch(tmp_path / "dest")


def test_fetch_refuses_a_member_over_the_size_cap_before_reading_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"x" * 100
    monkeypatch.setattr(fetch, "PINNED_SHA256", hashlib.sha256(content).hexdigest())
    monkeypatch.setattr(
        fetch, "BUNDLE_URL", _bundle(tmp_path / "b.zip", {fetch.WORKBOOK_NAME: content})
    )
    monkeypatch.setattr(fetch, "MAX_MEMBER_BYTES", 99)
    with pytest.raises(fetch.FetchRefused, match="cap"):
        fetch.fetch(tmp_path / "dest")
    monkeypatch.setattr(fetch, "MAX_MEMBER_BYTES", 100)  # the boundary itself is allowed
    assert fetch.fetch(tmp_path / "dest2").read_bytes() == content


def test_the_real_workbook_fits_under_the_size_cap() -> None:
    assert fetch.MAX_MEMBER_BYTES > 233_828 * 10


def test_dropped_examples_are_capped_at_five_and_the_count_is_not(
    tmp_path: Path, pinned: Any
) -> None:
    dates = [f"{day}-Mar" for day in range(1, 9)]  # eight corrupted symbols
    workbook = pinned(tmp_path / "w.xlsx", _sheet(["TBX5", *dates]))
    out = recipe.build(workbook, built_by=BUILT_BY, out=tmp_path / "pkg")
    (entry,) = _manifest(out)["dropped"]
    assert entry["count"] == 8
    assert entry["examples"] == sorted(dates)[:5]
