"""A synthetic-workbook writer: the smallest xlsx `feeds.xlsx` reads."""

from __future__ import annotations

import html
import zipfile
from pathlib import Path

# The strict OOXML namespace, which a reader assuming the transitional one reads
# as zero rows.
STRICT = "http://purl.oclc.org/ooxml/spreadsheetml/main"
STRICT_REL = "http://purl.oclc.org/ooxml/officeDocument/relationships"


def column_letters(index: int) -> str:
    letters = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(ord("A") + remainder) + letters
    return letters


def write_workbook(path: Path, rows: list[list[str]], sheet: str, *, strict: bool = True) -> Path:
    """``rows`` as inline strings on ``sheet``; empty cells are omitted (a sparse row)."""
    return write_workbook_sheets(path, {sheet: rows}, strict=strict)


def write_workbook_sheets(
    path: Path, sheets: dict[str, list[list[str]]], *, strict: bool = True
) -> Path:
    """Several sheets, in the dict's order, each written as ``write_workbook`` writes one."""
    ns = STRICT if strict else "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel = (
        STRICT_REL
        if strict
        else "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    )

    def data(rows: list[list[str]]) -> str:
        return "".join(
            f'<row r="{number}">'
            + "".join(
                f'<c r="{column_letters(index)}{number}" t="inlineStr">'
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
            f'<workbook xmlns="{ns}" xmlns:r="{rel}"><sheets>'
            + "".join(
                f'<sheet name="{html.escape(name)}" sheetId="{n}" r:id="rId{n}"/>'
                for n, name in enumerate(sheets, start=1)
            )
            + "</sheets></workbook>",
        )
        archive.writestr(
            "xl/_rels/workbook.xml.rels",
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + "".join(
                f'<Relationship Id="rId{n}" Target="worksheets/sheet{n}.xml"/>'
                for n in range(1, len(sheets) + 1)
            )
            + "</Relationships>",
        )
        for n, rows in enumerate(sheets.values(), start=1):
            archive.writestr(
                f"xl/worksheets/sheet{n}.xml",
                f'<worksheet xmlns="{ns}"><sheetData>{data(rows)}</sheetData></worksheet>',
            )
    return path
