"""The HGNC recipe on synthetic input: five genes, no network.

The expected rows are literals, written to equal what the atlas's
``scripts/convert_hgnc.py`` writes for the same raw rows (null stays an empty
cell, the RefSeq half is chosen by prefix), except the unpinned list, which the
converter refuses and the recipe keeps verbatim (the atlas's decision 3).
"""

import gzip
import hashlib
import json
import re
from datetime import date
from pathlib import Path
from typing import Any

import pytest

from feeds.package import PackageRefused
from sources.hgnc import recipe

BUILT_BY = {"recipe_commit": "a" * 40, "host": "login01", "hvantk_release": "4c90a241"}

RAW_HEADER = "hgnc_id\tsymbol\tname\tstatus\tuniprot_ids\tmane_select"
# The complete set quotes its multi-value cells; the hvantk export does not.
RAW_ROWS = [
    "HGNC:11604\tTBX5\tT-box transcription factor 5\tApproved\tQ99593"
    '\t"ENST00000405440.7|NM_181486.4"',
    'HGNC:9896\tRBM10\tRNA binding motif protein 10\tApproved\t"P0DW28|P98175"'
    '\t"ENST00000377604.9|NM_005676.5"',
    'HGNC:2488\tNKX2-5\tNK2 homeobox 5\tApproved\t"P52952|Q00000"\t"ENST00000329198.5|NM_004387.4"',
    "HGNC:50001\tLINC00001\tlong intergenic non-protein coding RNA 1\tApproved\t\t",
    "HGNC:50002\tTEST1\ttest gene 1\tApproved\tP12345\tENST00000999999.1",
    # File order is the reverse of sorted order, and HGNC:52 appears twice (rows != ids).
    "HGNC:51\tOLD1\told symbol\tSymbol Withdrawn\t\t",
    "HGNC:52\tOLD2\told entry\tEntry Withdrawn\t\t",
    "HGNC:52\tOLD2B\told entry again\tEntry Withdrawn\t\t",
    "HGNC:53\tOLD3\tanother old entry\tEntry Withdrawn\t\t",
]
HVANTK_HEADER = "\t".join(recipe.HVANTK_FIELDS)
HVANTK_ROWS = [
    # hgnc_id gene_symbol gene_name alias_symbols prev_symbols ensembl entrez location
    # uniprot_ids mane_select
    [
        "HGNC:50002",
        "TEST1",
        "test gene 1",
        "",
        "",
        "ENSG00000000001",
        "5",
        "1p1",
        "P12345",
        "ENST00000999999.1",
    ],
    [
        "HGNC:11604",
        "TBX5",
        "T-box transcription factor 5",
        "",
        "HOS",
        "ENSG00000089225",
        "6910",
        "12q24.21",
        "Q99593",
        "ENST00000405440.7|NM_181486.4",
    ],
    [
        "HGNC:9896",
        "RBM10",
        "RNA binding motif protein 10",
        "ZRANB5|DXS8237E",
        "",
        "ENSG00000182872",
        "8241",
        "Xp11.3",
        "P0DW28|P98175",
        "ENST00000377604.9|NM_005676.5",
    ],
    [
        "HGNC:2488",
        "NKX2-5",
        "NK2 homeobox 5",
        "CSX",
        "",
        "ENSG00000183072",
        "1482",
        "5q35.1",
        "P52952|Q00000",
        "ENST00000329198.5|NM_004387.4",
    ],
    [
        "HGNC:50001",
        "LINC00001",
        "long intergenic non-protein coding RNA 1",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
]


def _tsv(path: Path, header: str, rows: list[Any]) -> Path:
    lines = [header, *("\t".join(r) if isinstance(r, list) else r for r in rows)]
    data = ("\n".join(lines) + "\n").encode()
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(data, mtime=0))
    else:
        path.write_bytes(data)
    return path


def _inputs(
    tmp_path: Path, *, raw_rows: list[str] | None = None, hvantk_rows: list[Any] | None = None
):
    hvantk = _tsv(tmp_path / "hgnc_lookup.tsv.gz", HVANTK_HEADER, hvantk_rows or HVANTK_ROWS)
    raw = _tsv(tmp_path / "hgnc_complete_set.txt", RAW_HEADER, raw_rows or RAW_ROWS)
    return hvantk, raw


def _build(tmp_path: Path, hvantk: Path, raw: Path, **over: Any) -> Path:
    kwargs: dict[str, Any] = {
        "raw_url": "https://example.org/hgnc_complete_set.txt",
        "retrieved_on": date(2026, 11, 11),
        "release": "2026-10-16",
        "revision": 1,
        "built_by": BUILT_BY,
        "out": tmp_path / "pkg",
    }
    return recipe.build(hvantk, raw, **(kwargs | over))


def _rows(pkg: Path) -> list[list[str]]:
    text = gzip.decompress((pkg / "genes.tsv.gz").read_bytes()).decode()
    return [line.split("\t") for line in text.splitlines()]


def test_the_package_rows_are_the_rules_applied(tmp_path: Path) -> None:
    hvantk, raw = _inputs(tmp_path)
    pkg = _build(tmp_path, hvantk, raw)

    assert _rows(pkg) == [
        [
            "hgnc_id",
            "symbol",
            "name",
            "aliases",
            "prev_symbols",
            "ensembl_gene",
            "ncbi_gene",
            "locus",
            "uniprot",
            "mane_select",
        ],
        # sorted on the string: HGNC:11604 < HGNC:2488 < HGNC:50001 < HGNC:50002 < HGNC:9896
        [
            "HGNC:11604",
            "TBX5",
            "T-box transcription factor 5",
            "",
            "HOS",
            "ENSG00000089225",
            "6910",
            "12q24.21",
            "Q99593",
            "NM_181486.4",
        ],
        [
            "HGNC:2488",
            "NKX2-5",
            "NK2 homeobox 5",
            "CSX",
            "",
            "ENSG00000183072",
            "1482",
            "5q35.1",
            "P52952|Q00000",
            "NM_004387.4",
        ],
        [
            "HGNC:50001",
            "LINC00001",
            "long intergenic non-protein coding RNA 1",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
        ],
        ["HGNC:50002", "TEST1", "test gene 1", "", "", "ENSG00000000001", "5", "1p1", "P12345", ""],
        [
            "HGNC:9896",
            "RBM10",
            "RNA binding motif protein 10",
            "ZRANB5|DXS8237E",
            "",
            "ENSG00000182872",
            "8241",
            "Xp11.3",
            "P98175",
            "NM_005676.5",
        ],
    ]

    manifest = json.loads((pkg / "manifest.json").read_text())
    raw_sha = hashlib.sha256(raw.read_bytes()).hexdigest()
    assert manifest["package"] == "hgnc"
    assert manifest["version"] == "2026-10-16"
    assert manifest["ids"] == "hgnc_id"
    assert manifest["coverage"] == "every_gene"
    assert manifest["inputs"] == [
        {
            "name": "hgnc_complete_set.txt",
            "url": "https://example.org/hgnc_complete_set.txt",
            "sha256": raw_sha,
            "retrieved_on": "2026-11-11",
            "upstream_release": "2026-10-16",
        }
    ]
    assert manifest["input_counts"] == {"rows": 9, "genes": 8}
    assert manifest["dropped"] == [
        {
            "reason": "status is Entry Withdrawn, not Approved; hvantk keeps Approved genes only",
            "count": 3,
            "examples": ["HGNC:52", "HGNC:53"],
        },
        {
            "reason": "status is Symbol Withdrawn, not Approved; hvantk keeps Approved genes only",
            "count": 1,
            "examples": ["HGNC:51"],
        },
    ]
    assert manifest["source"] == {
        "name": "HUGO Gene Nomenclature Committee",
        "version": "2026-10-16",
        "retrieved_on": "2026-11-11",
        "url": "https://www.genenames.org/",
        "licence": "CC0-1.0",
        "redistribution": "permitted",
    }
    assert "CC0" in manifest["licence_evidence"]
    assert manifest["built_by"] == BUILT_BY
    (table,) = manifest["tables"]
    assert (table["name"], table["format_version"], table["rows"], table["genes"]) == (
        "genes",
        2,
        5,
        5,
    )


def test_a_revision_is_part_of_the_version(tmp_path: Path) -> None:
    hvantk, raw = _inputs(tmp_path)
    pkg = _build(tmp_path, hvantk, raw, revision=2)
    assert json.loads((pkg / "manifest.json").read_text())["version"] == "2026-10-16-r2"


@pytest.mark.parametrize(
    ("cell", "message"),
    [
        pytest.param(
            "P0DW28", r"pinned UniProt accession P98175 is not among.*\(P0DW28\)", id="stale"
        ),
        pytest.param("", r"pinned UniProt accession P98175 is not among.*\(None\)", id="empty"),
    ],
)
def test_a_pin_hgnc_no_longer_lists_refuses_the_build(
    tmp_path: Path, cell: str, message: str
) -> None:
    rows = [list(r) for r in HVANTK_ROWS]
    next(r for r in rows if r[0] == "HGNC:9896")[8] = cell
    hvantk, raw = _inputs(tmp_path, hvantk_rows=rows)
    with pytest.raises(PackageRefused, match=message):
        _build(tmp_path, hvantk, raw)
    assert not (tmp_path / "pkg").exists()


def test_a_pinned_gene_missing_from_the_table_refuses_the_build(tmp_path: Path) -> None:
    rows = [r for r in HVANTK_ROWS if r[0] != "HGNC:9896"]
    raw_rows = [r for r in RAW_ROWS if not r.startswith("HGNC:9896")]
    hvantk, raw = _inputs(tmp_path, raw_rows=raw_rows, hvantk_rows=rows)
    with pytest.raises(PackageRefused, match="pinned gene.*HGNC:9896"):
        _build(tmp_path, hvantk, raw)


def test_rows_hvantk_dropped_must_be_accounted_for(tmp_path: Path) -> None:
    # An Approved gene hvantk lost: the complete set has it, the table does not.
    rows = [r for r in HVANTK_ROWS if r[0] != "HGNC:50002"]
    hvantk, raw = _inputs(tmp_path, hvantk_rows=rows)
    with pytest.raises(PackageRefused, match=r"1 missing \(e\.g\. \['HGNC:50002'\]\)"):
        _build(tmp_path, hvantk, raw)

    # A gene the table holds although the complete set withdrew it.
    extra = [*HVANTK_ROWS, ["HGNC:51", "OLD1", "old symbol", "", "", "", "", "", "", ""]]
    hvantk, raw = _inputs(tmp_path, hvantk_rows=extra)
    with pytest.raises(PackageRefused, match=r"1 not Approved there \(e\.g\. \['HGNC:51'\]\)"):
        _build(tmp_path, hvantk, raw)

    # A repeated row hides a missing one behind a matching count.
    dup = [*HVANTK_ROWS[:-1], HVANTK_ROWS[0]]
    hvantk, raw = _inputs(tmp_path, hvantk_rows=dup)
    with pytest.raises(PackageRefused, match="repeats an hgnc_id"):
        _build(tmp_path, hvantk, raw)

    # The complete set repeating an Approved row would count as a dropped one.
    hvantk, raw = _inputs(tmp_path, raw_rows=[*RAW_ROWS, RAW_ROWS[0]])
    with pytest.raises(PackageRefused, match=r"5 complete-set rows are absent.*but 4 are not"):
        _build(tmp_path, hvantk, raw)
    assert not (tmp_path / "pkg").exists()


def test_a_column_the_recipe_reads_that_hvantk_did_not_export_is_named(tmp_path: Path) -> None:
    header = "\t".join(f for f in recipe.HVANTK_FIELDS if f != "mane_select")
    rows = [r[:-1] for r in HVANTK_ROWS]
    hvantk = _tsv(tmp_path / "hgnc_lookup.tsv.gz", header, rows)
    raw = _tsv(tmp_path / "hgnc_complete_set.txt", RAW_HEADER, RAW_ROWS)
    with pytest.raises(PackageRefused, match="lacks column.*mane_select"):
        _build(tmp_path, hvantk, raw)


def test_the_sbatch_exports_exactly_the_fields_the_recipe_reads() -> None:
    sbatch = (Path(recipe.__file__).parent / "hpc" / "hgnc_build.sbatch").read_text()
    match = re.search(r'^FIELDS="([^"]+)"$', sbatch, re.M)
    assert match is not None
    assert tuple(match.group(1).split(",")) == recipe.HVANTK_FIELDS
    assert "--mail-user" not in re.sub(r"^#.*$", "", sbatch, flags=re.M)
    assert not re.search(r"^#SBATCH --mail-user", sbatch, re.M)


def test_a_quoted_hvantk_cell_is_refused_not_unquoted(tmp_path: Path) -> None:
    rows = [list(r) for r in HVANTK_ROWS]
    rows[0][2] = '"quoted" name'
    hvantk, raw = _inputs(tmp_path, hvantk_rows=rows)
    with pytest.raises(PackageRefused, match="opens with a double quote"):
        _build(tmp_path, hvantk, raw)


@pytest.mark.parametrize(
    ("raw_name", "release", "message"),
    [
        ("hgnc_complete_set_2026-10-06.txt", "2026-10-16", "2026-10-06 archive.*2026-10-16"),
        ("hgnc_complete_set.txt", "October 2026", "not an ISO date"),
        ("hgnc_complete_set.txt", "2026-10-1", "not an ISO date"),
    ],
)
def test_the_release_is_an_iso_date_and_matches_an_archive_name(
    tmp_path: Path, raw_name: str, release: str, message: str
) -> None:
    hvantk, raw = _inputs(tmp_path)
    renamed = raw.rename(tmp_path / raw_name)
    with pytest.raises(PackageRefused, match=message):
        _build(tmp_path, hvantk, renamed, release=release)
    ok = tmp_path / "hgnc_complete_set_2026-10-16.txt"
    renamed.rename(ok)
    _build(tmp_path, hvantk, ok, release="2026-10-16")
