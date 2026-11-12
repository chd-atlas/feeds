"""A minimal read-only xlsx reader (openpyxl is not a dependency).

Shared by the recipes that read a journal supplement workbook. The body is the
one `sources/sierant2025_ds6` carried, which in turn matched the atlas's
`scripts/convert_burden.py` (at 4b1f8f6) line for line; both were compared on
real workbooks before this module replaced them (see the commit message).

The namespace is derived from the root tag rather than hardcoded: PNAS
supplements use the *strict* OOXML namespace
(`http://purl.oclc.org/ooxml/spreadsheetml/main`), and a reader assuming the
transitional `schemas.openxmlformats.org` one silently finds zero rows.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from typing import Final
from xml.etree import ElementTree as ET

_CELL_REF: Final = re.compile(r"([A-Z]+)\d+")


def _namespace(root: ET.Element) -> str:
    return root.tag.split("}")[0][1:] if root.tag.startswith("{") else ""


def _tag(namespace: str, name: str) -> str:
    return f"{{{namespace}}}{name}" if namespace else name


def _column_index(reference: str) -> int:
    """Zero-based column index of a cell reference such as "AB12"."""
    match = _CELL_REF.match(reference)
    if match is None:
        raise ValueError(f"unparseable cell reference {reference!r}")
    index = 0
    for character in match.group(1):
        index = index * 26 + (ord(character) - 64)
    return index - 1


def read_sheet(path: Path, name: str) -> list[dict[str, str]]:
    """One worksheet as a list of dicts keyed on its header row.

    Only correct for a sheet whose first row *is* the header. PNAS supplements
    put a title, a caption and a three-level header above the data, which keys
    every column on `""`; those callers use `read_sheet_rows` and index by position.
    """
    rows = read_sheet_rows(path, name)
    header = rows[0]
    named = [h for h in header if h]
    repeated = sorted({h for h in named if named.count(h) > 1})
    if repeated:
        # A dict keyed on the header would keep only the last of two columns
        # sharing a name, under a package that builds.
        raise ValueError(
            f"sheet {name!r} repeats the header {', '.join(map(repr, repeated))}; "
            "a record keyed on it would silently drop all but one of those columns"
        )
    return [
        {name: (row[i] if i < len(row) else "") for i, name in enumerate(header)}
        for row in rows[1:]
    ]


def read_sheet_rows(path: Path, name: str) -> list[list[str]]:
    """One worksheet as raw rows of cell text, blanks preserved by position.

    Position is the point: a sheet with a merged multi-level header has no row
    that can key a dict, and the columns still have to be addressable. Cells are
    placed by their `r` reference rather than by order of appearance, because a
    sparse row omits its empty cells entirely -- so a run of blanks in the middle
    of a row would otherwise shift every column after it left.
    """
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        workbook_ns = _namespace(workbook)
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {
            element.get("Id"): element.get("Target")
            for element in relationships.findall(_tag(_namespace(relationships), "Relationship"))
        }

        sheets = workbook.find(_tag(workbook_ns, "sheets"))
        if sheets is None:
            raise ValueError(f"{path} declares no sheets")
        target: str | None = None
        for sheet in sheets.findall(_tag(workbook_ns, "sheet")):
            if sheet.get("name") != name:
                continue
            # The r:id attribute is namespaced, and the namespace differs
            # between the strict and transitional OOXML flavours, so it is found
            # by suffix rather than spelled out.
            key = next(attribute for attribute in sheet.attrib if attribute.endswith("}id"))
            target = targets[sheet.get(key)]
        if target is None:
            raise ValueError(f"{path} has no sheet named {name!r}")

        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            table = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            table_ns = _namespace(table)
            shared = [
                "".join(node.text or "" for node in item.iter(_tag(table_ns, "t")))
                for item in table.findall(_tag(table_ns, "si"))
            ]

        worksheet = ET.fromstring(archive.read("xl/" + target.lstrip("/").removeprefix("xl/")))
        sheet_ns = _namespace(worksheet)
        data = worksheet.find(_tag(sheet_ns, "sheetData"))
        if data is None:
            raise ValueError(f"sheet {name!r} of {path} has no sheetData")

        rows: list[list[str]] = []
        for row in data.findall(_tag(sheet_ns, "row")):
            cells: dict[int, str] = {}
            for cell in row.findall(_tag(sheet_ns, "c")):
                reference = cell.get("r")
                index = _column_index(reference) if reference else len(cells)
                kind = cell.get("t")
                if kind == "inlineStr":
                    inline = cell.find(_tag(sheet_ns, "is"))
                    value = (
                        "".join(node.text or "" for node in inline.iter(_tag(sheet_ns, "t")))
                        if inline is not None
                        else ""
                    )
                else:
                    node = cell.find(_tag(sheet_ns, "v"))
                    value = (node.text or "") if node is not None else ""
                    if kind == "s" and value:
                        value = shared[int(value)]
                cells[index] = value
            rows.append([cells.get(i, "") for i in range(max(cells) + 1 if cells else 0)])

    return rows
