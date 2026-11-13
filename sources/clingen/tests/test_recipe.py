"""The ClinGen recipe on synthetic inputs: no network, no committed mirror.

The banner and rows are the atlas's ``tests/unit/test_convert_clingen.py`` fixture,
ported. The ClinGen CSV preamble is a trap and the trap is silent: lines 1-3 are a
title block, line 4 and line 6 are ``+++`` separator rows, and line 5 is the real
header. A reader with default arguments takes line 1 as the header and produces ten
columns named after the title, with no exception and no empty table.
"""

import gzip
import hashlib
import json
import re
from datetime import date
from pathlib import Path

import pytest

from feeds.package import PackageRefused
from sources.clingen import recipe

BUILT_BY = {"host": "login01"}
RECIPE_COMMIT = "a" * 40

COLUMNS = [
    "gene",
    "gene_symbol",
    "disease",
    "disease_label",
    "moi",
    "sop",
    "classification",
    "classification_date",
    "gcep",
    "report_url",
]

BANNER = (
    '"CLINGEN GENE DISEASE VALIDITY CURATIONS","","","","","","","","",""\n'
    '"FILE CREATED: 2026-11-11","","","","","","","","",""\n'
    '"WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity","","","","","","","","",""\n'
    '"+++++++++++","++++++++++++++","+++++++++++++","++++++++++++++++++","+++++++++",'
    '"+++++++++","++++++++++++++","+++++++++++++","+++++++++++++++++++","+++++++++++++++++++"\n'
    '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP",'
    '"CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n'
    '"+++++++++++","++++++++++++++","+++++++++++++","++++++++++++++++++","+++++++++",'
    '"+++++++++","++++++++++++++","+++++++++++++","+++++++++++++++++++","+++++++++++++++++++"\n'
)

# Out of order on purpose: TBX5 is HGNC:11604 and GATA4 is HGNC:4173, so symbol order
# and id order disagree, and string comparison puts "HGNC:11604" before "HGNC:4173".
ROWS = (
    '"TBX5","HGNC:11604","Holt-Oram syndrome","MONDO:0007732","AD","SOP11","Definitive",'
    '"https://example.org/a","2025-03-25T16:00:00.000Z","Syndromic Disorders GCEP"\n'
    '"GATA4","HGNC:4173","congenital heart disease","MONDO:0005453","AD","SOP9","Limited",'
    '"https://example.org/b","2020-01-01T16:00:00.000Z","Congenital Heart Disease GCEP"\n'
    '"GATA4","HGNC:4173","atrial septal defect","MONDO:0000001","AR","SOP9","Disputed",'
    '"https://example.org/d","2022-01-01T16:00:00.000Z","Congenital Heart Disease GCEP"\n'
    '"GATA4","HGNC:4173","atrial septal defect","MONDO:0000001","AD","SOP9","Moderate",'
    '"https://example.org/c","2021-01-01T16:00:00.000Z","Congenital Heart Disease GCEP"\n'
)


def _csv(tmp_path: Path, text: str = BANNER + ROWS) -> Path:
    path = tmp_path / "clingen.csv"
    path.write_text(text, encoding="utf-8")
    return path


def _read_package(root: Path) -> tuple[dict[str, object], list[list[str]]]:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    raw = gzip.decompress((root / "clingen_validity.tsv.gz").read_bytes()).decode("utf-8")
    return manifest, [line.split("\t") for line in raw.splitlines()]


def test_the_banner_is_the_version_and_separator_rows_are_not_data(tmp_path: Path) -> None:
    banner, rows = recipe.read_rows(_csv(tmp_path))
    assert banner == "FILE CREATED: 2026-11-11"
    assert len(rows) == 4, "four data rows, no separator rows"
    assert not any("+++" in cell for row in rows for cell in row)

    out = recipe.build(
        _csv(tmp_path),
        retrieved_on=date(2026, 8, 4),
        revision=1,
        built_by=BUILT_BY,
        out=tmp_path / "pkg",
        commit_of=lambda: RECIPE_COMMIT,
    )
    manifest, table = _read_package(out)
    assert manifest["version"] == "2026-11-11"
    assert manifest["ids"] == "hgnc_id"
    assert manifest["coverage"] == "every_gene"
    assert manifest["source"]["version"] == "FILE CREATED: 2026-11-11"  # type: ignore[index]
    assert manifest["source"]["retrieved_on"] == "2026-08-04"  # type: ignore[index]
    assert manifest["built_by"] == {**BUILT_BY, "recipe_commit": RECIPE_COMMIT}
    assert manifest["input_counts"] == {"rows": 4, "genes": 2}
    (entry,) = manifest["inputs"]  # type: ignore[misc]
    assert entry["upstream_release"] == "FILE CREATED: 2026-11-11"
    assert entry["retrieved_on"] == "2026-08-04"
    assert entry["sha256"] == hashlib.sha256(_csv(tmp_path).read_bytes()).hexdigest()
    assert table[0] == COLUMNS
    assert len(table) == 5


def test_rows_are_written_in_table_order_and_sorted_by_the_key(tmp_path: Path) -> None:
    _, rows = recipe.read_rows(_csv(tmp_path))
    assert [(r[0], r[2], r[4]) for r in rows] == [
        ("HGNC:11604", "MONDO:0007732", "AD"),
        ("HGNC:4173", "MONDO:0000001", "AD"),
        ("HGNC:4173", "MONDO:0000001", "AR"),  # input order is AR first: moi decides
        ("HGNC:4173", "MONDO:0005453", "AD"),
    ]
    # Table column order, not the CSV's: ONLINE REPORT is last, GENE SYMBOL second.
    assert rows[0] == (
        "HGNC:11604",
        "TBX5",
        "MONDO:0007732",
        "Holt-Oram syndrome",
        "AD",
        "SOP11",
        "Definitive",
        "2025-03-25T16:00:00.000Z",
        "Syndromic Disorders GCEP",
        "https://example.org/a",
    )


def test_a_renamed_upstream_column_stops_the_build(tmp_path: Path) -> None:
    """Loud, not partial: a renamed column must not yield a table missing a field."""
    path = _csv(tmp_path, BANNER.replace('"MOI"', '"MODE OF INHERITANCE"') + ROWS)
    with pytest.raises(PackageRefused, match="MOI"):
        recipe.read_rows(path)
    with pytest.raises(PackageRefused, match="MOI"):
        recipe.build(
            path,
            retrieved_on=date(2026, 8, 4),
            revision=1,
            built_by=BUILT_BY,
            out=tmp_path / "pkg",
            commit_of=lambda: RECIPE_COMMIT,
        )
    assert not (tmp_path / "pkg").exists()


def test_a_banner_without_a_creation_date_is_refused(tmp_path: Path) -> None:
    path = _csv(tmp_path, BANNER.replace("FILE CREATED: 2026-11-11", "CREATED SOMETIME") + ROWS)
    with pytest.raises(PackageRefused, match="FILE CREATED"):
        recipe.read_rows(path)


# FOXE3 / familial TAAD / AD, curated by two GCEPs: the real pair in ClinGen's 2026-11-13
# file (F28). The Heritable (HTAAD) row comes first here, the order that sorts wrong:
# "Hereditary" < "Heritable", so a sort that ignores the panel keeps input order.
_FOXE3_HTAAD = (
    '"FOXE3","HGNC:3808","familial thoracic aortic aneurysm and aortic dissection",'
    '"MONDO:0019625","AD","SOP4","Moderate",'
    '"https://search.clinicalgenome.org/kb/gene-validity/CGGCIEX:assertion_8261",'
    '"2016-06-15T00:00:00",'
    '"Heritable Thoracic Aortic Aneurysm and Dissection Gene Curation Expert Panel"\n'
)
_FOXE3_HCVD = (
    '"FOXE3","HGNC:3808","familial thoracic aortic aneurysm and aortic dissection",'
    '"MONDO:0019625","AD","SOP12","Limited",'
    '"https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_982f5f02-596f-4dfa-'
    'a019-147934f6509b-2026-06-19T160000.000Z","2026-06-19T16:00:00.000Z",'
    '"Hereditary Cardiovascular Disease Gene Curation Expert Panel"\n'
)


def test_two_panels_on_one_triple_are_both_written_in_panel_order(tmp_path: Path) -> None:
    out = recipe.build(
        _csv(tmp_path, BANNER + _FOXE3_HTAAD + _FOXE3_HCVD),
        retrieved_on=date(2026, 11, 13),
        revision=1,
        built_by=BUILT_BY,
        out=tmp_path / "pkg",
        commit_of=lambda: RECIPE_COMMIT,
    )
    manifest, table = _read_package(out)
    assert [row[COLUMNS.index("gcep")] for row in table[1:]] == [
        "Hereditary Cardiovascular Disease Gene Curation Expert Panel",
        "Heritable Thoracic Aortic Aneurysm and Dissection Gene Curation Expert Panel",
    ]
    assert [row[COLUMNS.index("classification")] for row in table[1:]] == ["Limited", "Moderate"]
    (entry,) = manifest["tables"]  # type: ignore[misc]
    assert entry["format_version"] == 2
    assert manifest["input_counts"] == {"rows": 2, "genes": 1}


_FOXE3_HTAAD_RELABELLED = _FOXE3_HTAAD.replace(
    "Heritable Thoracic Aortic Aneurysm and Dissection Gene Curation Expert Panel",
    "Hereditary Cardiovascular Disease Gene Curation Expert Panel",
)


@pytest.mark.parametrize(
    "first",
    [_FOXE3_HTAAD_RELABELLED, _FOXE3_HCVD],
    ids=["another-record-same-panel", "verbatim-repeat"],
)
def test_one_panel_twice_on_a_triple_is_refused(tmp_path: Path, first: str) -> None:
    """The same GCEP twice on one triple is a duplicate the owner rules on, never deduplicated.

    Both shapes: a different record under the same panel, and the same line repeated
    verbatim, which a set-based tidy would silently collapse.
    """
    assert _FOXE3_HTAAD_RELABELLED != _FOXE3_HTAAD
    # The v2 key, panel included: under v1 the same file is refused on the triple alone.
    key = "('HGNC:3808', 'MONDO:0019625', 'AD', 'Hereditary Cardiovascular Disease"
    with pytest.raises(PackageRefused, match=re.escape(f"primary key {key}")):
        recipe.build(
            _csv(tmp_path, BANNER + first + _FOXE3_HCVD),
            retrieved_on=date(2026, 11, 13),
            revision=1,
            built_by=BUILT_BY,
            out=tmp_path / "pkg",
            commit_of=lambda: RECIPE_COMMIT,
        )
    assert not (tmp_path / "pkg").exists()
