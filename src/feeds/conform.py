"""Check a table against the atlas's exported table schema (``schemas/atlas/tables``)."""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from typing import Any

_BOOLEANS = frozenset({"true", "false"})
# Strict forms: no underscores, no surrounding whitespace (polars' i64/f64 parse, atlas TBL000).
_INTEGER = re.compile(r"[+-]?[0-9]+")
_I64 = (-(2**63), 2**63 - 1)  # polars' i64 parse; beyond it the atlas's TBL000 -> ING010
_NUMBER = re.compile(r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
_NON_FINITE = frozenset({"nan", "inf", "+inf", "-inf", "infinity", "+infinity", "-infinity"})


def _problem(kind: str, cell: str) -> str | None:
    if kind == "integer":
        if not _INTEGER.fullmatch(cell):
            return "is not an integer"
        try:
            inside = _I64[0] <= int(cell) <= _I64[1]
        except ValueError:  # Python refuses to parse more than 4,300 digits
            inside = False
        if not inside:
            return "is outside the 64-bit integer range"
    elif kind == "number":
        if cell.lower() in _NON_FINITE:
            return "is not finite"
        if not _NUMBER.fullmatch(cell):
            return "is not a number"
        if not math.isfinite(float(cell)):
            return "is not finite"
    elif kind == "boolean" and cell not in _BOOLEANS:
        return "is not 'true' or 'false'"
    return None


def check_table(
    columns: Sequence[str],
    rows: Sequence[Sequence[str | None]],
    table_schema: dict[str, Any],
    *,
    ids: str,
) -> list[str]:
    """Human-readable errors (row = 1-based data row); an empty list means it conforms."""
    fields: list[dict[str, Any]] = table_schema["fields"]
    names = [f["name"] for f in fields]
    if list(columns) != names:
        return [f"header {list(columns)} differs from the schema's {names}"]
    missing = set(table_schema.get("missingValues", [""]))  # only the literal "", never " "
    gene_column = table_schema.get("geneColumn")
    errors: list[str] = []
    for number, row in enumerate(rows, start=1):
        if len(row) != len(names):
            errors.append(f"row {number}: {len(row)} cells, expected {len(names)}")
            continue
        for field, cell in zip(fields, row, strict=False):
            name, rules = field["name"], field.get("constraints", {})
            if cell is None or cell in missing:
                if rules.get("required"):
                    errors.append(f"row {number}, column '{name}': required but empty")
                continue
            if "enum" in rules and cell not in rules["enum"]:
                errors.append(f"row {number}, column '{name}': '{cell}' is not in {rules['enum']}")
            skip = name == gene_column and ids != "hgnc_id"
            if "pattern" in rules and not skip and not re.fullmatch(rules["pattern"], cell):
                errors.append(f"row {number}, column '{name}': '{cell}' !~ {rules['pattern']}")
            if (why := _problem(field["type"], cell)) is not None:
                errors.append(f"row {number}, column '{name}': '{cell}' {why}")
            number_value = None
            if field["type"] in ("integer", "number") and _problem(field["type"], cell) is None:
                number_value = float(cell)
            if number_value is not None:
                if "minimum" in rules and number_value < rules["minimum"]:
                    errors.append(f"row {number}, column '{name}': {cell} < {rules['minimum']}")
                if "maximum" in rules and number_value > rules["maximum"]:
                    errors.append(f"row {number}, column '{name}': {cell} > {rules['maximum']}")
    key = table_schema.get("primaryKey")
    if key:
        index = [names.index(k) for k in key]
        seen: dict[tuple[str, ...], int] = {}
        for number, row in enumerate(rows, start=1):
            if len(row) != len(names):
                continue
            tup = tuple((row[i] or "").strip() for i in index)
            if tup in seen:
                errors.append(f"row {number}: primary key {tup} repeats row {seen[tup]}")
            else:
                seen[tup] = number
    return errors
