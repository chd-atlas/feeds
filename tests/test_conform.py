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


def test_numbers_are_parsed_strictly_and_must_be_finite() -> None:
    def errs(col: int, cell: str) -> list[str]:
        row: list[str | None] = ["HGNC:1", None, None, None, None]
        row[col] = cell
        return check_table(COLS, [row], SCHEMA, ids="hgnc_id")

    for cell in ["1_000", " 100 ", "1.5", "", "+"]:
        if cell:
            assert errs(2, cell), cell
    for cell in ["nan", "inf", "-inf", "NaN", "Infinity", "1e999", "1_0", " 0.5", "0.5 "]:
        assert errs(3, cell), cell
    assert not errs(2, "+3") and not errs(3, "1e-3") and not errs(3, ".5") and not errs(3, "1.")
    assert "not finite" in errs(3, "nan")[0]


def test_only_the_empty_string_is_missing() -> None:
    sch: dict[str, Any] = {
        "missingValues": [""],
        "fields": [
            {"name": "k", "type": "string", "constraints": {"enum": ["a"], "required": True}}
        ],
    }
    assert "not in" in check_table(["k"], [[" "]], sch, ids="symbol")[0]
    assert check_table(["k"], [[""]], sch, ids="symbol")
    assert check_table(["k"], [["a"]], sch, ids="symbol") == []


def test_a_ragged_row_is_reported_not_raised() -> None:
    errors = check_table(COLS, [["HGNC:1", "a"]], SCHEMA, ids="hgnc_id")
    assert errors == ["row 1: 2 cells, expected 5"]


def test_a_number_that_overflows_to_infinity_is_not_finite() -> None:
    sch: dict[str, Any] = {
        "fields": [{"name": "x", "type": "number", "constraints": {"required": True}}]
    }
    assert "not finite" in check_table(["x"], [["1e999"]], sch, ids="symbol")[0]
    assert check_table(["x"], [["1e300"]], sch, ids="symbol") == []


def test_an_integer_must_fit_the_64_bits_the_atlas_parses() -> None:
    schema = {"fields": [{"name": "n", "type": "integer", "constraints": {}}]}
    fits = ["9223372036854775807", "-9223372036854775808", "+0"]
    over = ["9223372036854775808", "-9223372036854775809", "9" * 5000]
    assert check_table(["n"], [[c] for c in fits], schema, ids="symbol") == []
    errors = check_table(["n"], [[c] for c in over], schema, ids="symbol")
    assert len(errors) == 3 and all("64-bit" in e for e in errors)
