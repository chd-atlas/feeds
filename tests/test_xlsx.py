import zipfile
from pathlib import Path

import pytest

from feeds import xlsx
from tests.xlsx_fixture import write_workbook

ROWS = [["gene", "p"], ["TBX5", "0.5"], ["GATA4", "1e-3"]]


@pytest.mark.parametrize("strict", [True, False])
def test_a_workbook_reads_its_rows_in_either_namespace(tmp_path: Path, strict: bool) -> None:
    path = write_workbook(tmp_path / "w.xlsx", ROWS, "S", strict=strict)
    assert xlsx.read_sheet_rows(path, "S") == ROWS
    assert xlsx.read_sheet(path, "S") == [
        {"gene": "TBX5", "p": "0.5"},
        {"gene": "GATA4", "p": "1e-3"},
    ]


def test_a_sparse_row_keeps_column_positions(tmp_path: Path) -> None:
    path = write_workbook(tmp_path / "w.xlsx", [["a", "", "", "d"], ["", "b"]], "S")
    assert xlsx.read_sheet_rows(path, "S") == [["a", "", "", "d"], ["", "b"]]


def test_shared_strings_are_resolved(tmp_path: Path) -> None:
    path = tmp_path / "w.xlsx"
    ns = "http://purl.oclc.org/ooxml/spreadsheetml/main"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr(
            "xl/workbook.xml",
            f'<workbook xmlns="{ns}" xmlns:r="http://purl.oclc.org/ooxml/officeDocument/relationships">'
            '<sheets><sheet name="S" sheetId="1" r:id="rId1"/></sheets></workbook>',
        )
        z.writestr(
            "xl/_rels/workbook.xml.rels",
            '<Relationships xmlns="x">'
            '<Relationship Id="rId1" Target="/xl/worksheets/sheet1.xml"/></Relationships>',
        )
        z.writestr(
            "xl/sharedStrings.xml",
            f'<sst xmlns="{ns}"><si><t>one</t></si><si><t>t</t><t>wo</t></si></sst>',
        )
        z.writestr(
            "xl/worksheets/sheet1.xml",
            f'<worksheet xmlns="{ns}"><sheetData><row><c r="A1" t="s"><v>0</v></c>'
            '<c r="C1" t="s"><v>1</v></c><c r="D1"><v>7</v></c></row></sheetData></worksheet>',
        )
    assert xlsx.read_sheet_rows(path, "S") == [["one", "", "two", "7"]]


def test_a_missing_sheet_is_refused(tmp_path: Path) -> None:
    path = write_workbook(tmp_path / "w.xlsx", ROWS, "S")
    with pytest.raises(ValueError, match="no sheet named 'T'"):
        xlsx.read_sheet_rows(path, "T")
