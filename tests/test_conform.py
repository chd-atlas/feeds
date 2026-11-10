from __future__ import annotations

from typing import Any

from feeds.conform import check_table

SCHEMA: dict[str, Any] = {
    "geneColumn": "g",
    "missingValues": [""],
    "primaryKey": ["g", "k"],
    "fields": [
        {
            "name": "g",
            "type": "string",
            "constraints": {"required": True, "pattern": "^HGNC:\\d+$"},
        },
        {"name": "k", "type": "string", "constraints": {"enum": ["a", "b"], "required": False}},
        {"name": "n", "type": "integer", "constraints": {"minimum": 0, "required": False}},
        {"name": "x", "type": "number", "constraints": {"maximum": 1, "required": False}},
        {"name": "f", "type": "boolean", "constraints": {"required": False}},
    ],
}
COLS = ["g", "k", "n", "x", "f"]


def test_a_good_table_has_no_errors() -> None:
    rows = [["HGNC:1", "a", "3", "0.5", "true"], ["HGNC:1", None, None, "", "false"]]
    assert check_table(COLS, rows, SCHEMA, ids="hgnc_id") == []


def test_each_violation_is_reported_with_its_row_and_column() -> None:
    rows = [
        ["", "z", "x", "2", "yes"],
        ["HGNC:1", "a", "-1", "abc", None],
        ["HGNC:1", "a", None, None, None],
    ]
    text = "\n".join(check_table(COLS, rows, SCHEMA, ids="hgnc_id"))
    for needle in [
        "row 1, column 'g': required",
        "row 1, column 'k'",
        "row 1, column 'n': 'x' is not an integer",
        "row 1, column 'x': 2 > 1",
        "row 1, column 'f'",
        "row 2, column 'n': -1 < 0",
        "row 2, column 'x': 'abc' is not a number",
        "row 3: primary key",
    ]:
        assert needle in text, needle


def test_the_header_must_match_in_order() -> None:
    assert check_table(["k", "g", "n", "x", "f"], [], SCHEMA, ids="hgnc_id")


def test_the_gene_pattern_is_skipped_for_every_id_kind_but_hgnc() -> None:
    rows = [["TBX5", None, None, None, None]]
    assert check_table(COLS, rows, SCHEMA, ids="symbol") == []
    assert check_table(COLS, rows, SCHEMA, ids="ensembl_gene") == []
    assert check_table(COLS, rows, SCHEMA, ids="hgnc_id")
    # required still applies to a non-HGNC gene column
    assert check_table(COLS, [["", None, None, None, None]], SCHEMA, ids="symbol")
