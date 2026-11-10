from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pytest

from feeds.package import PackageRefused, TableOut, write_package
from tests.conftest import BURDEN_COLUMNS, burden_row, burden_table, package_kwargs


def _bytes(directory: Path) -> dict[str, bytes]:
    return {p.name: p.read_bytes() for p in sorted(directory.iterdir())}


def test_a_package_is_byte_identical_across_two_runs(tmp_path: Path) -> None:
    write_package(tmp_path / "a", **package_kwargs())
    write_package(tmp_path / "b", **package_kwargs())
    assert _bytes(tmp_path / "a") == _bytes(tmp_path / "b")
    blob = (tmp_path / "a" / "burden.tsv.gz").read_bytes()
    assert blob[4:8] == b"\0\0\0\0"  # gzip mtime field
    assert blob[3] == 0  # no file name in the header
    text = gzip.decompress(blob).decode()
    assert text.endswith("\n") and "\r" not in text and '"' not in text
    manifest = json.loads((tmp_path / "a" / "manifest.json").read_text())
    entry = manifest["tables"][0]
    assert entry["sha256"] == hashlib.sha256(blob).hexdigest()
    assert (entry["rows"], entry["genes"], entry["file"]) == (2, 2, "burden.tsv.gz")
    assert (tmp_path / "a" / "manifest.json").read_text().endswith("}\n")


@pytest.mark.parametrize("cell", ["a\tb", "a\nb", "a\rb", '"open'])
def test_the_writer_refuses_a_cell_the_atlas_would_misread(out: Path, cell: str) -> None:
    rows = [burden_row("TBX5", method_note=cell), burden_row("GATA4")]
    with pytest.raises(PackageRefused):
        write_package(out, **package_kwargs([burden_table(rows)]))
    assert not out.exists()


def test_a_package_conforms_to_the_vendored_schemas(tmp_path: Path) -> None:
    write_package(tmp_path / "ok", **package_kwargs(ids="symbol"))
    with pytest.raises(PackageRefused, match="HGNC"):
        write_package(tmp_path / "bad", **package_kwargs(ids="hgnc_id"))
    assert not (tmp_path / "bad").exists()
    good = [burden_row("HGNC:1100"), burden_row("HGNC:11634")]
    write_package(tmp_path / "hgnc", **package_kwargs([burden_table(good)], ids="hgnc_id"))


def test_genes_counts_distinct_stripped_non_empty_values(out: Path) -> None:
    rows = [burden_row("TBX5"), burden_row(" TBX5 ", study="PMID:2"), burden_row("GATA4")]
    write_package(out, **package_kwargs([burden_table(rows)]))
    assert json.loads((out / "manifest.json").read_text())["tables"][0]["genes"] == 2


def _quantiles(shard: str | None) -> TableOut:
    cols = ("dataset", "tissue", "stage", "percentile", "value", "unit", "n_genes")
    return TableOut(
        "profile_quantiles",
        1,
        cols,
        [["d", "heart", "4wpc", "50", "1.5", "tpm", "10"]],
        None,
        shard,
    )


def test_a_shard_names_its_file_and_layout_mismatches_are_refused(tmp_path: Path) -> None:
    write_package(tmp_path / "ok", **package_kwargs([_quantiles("s1")]))
    assert (tmp_path / "ok" / "profile_quantiles.s1.tsv.gz").is_file()
    with pytest.raises(PackageRefused, match="sharded and gives no shard"):
        write_package(tmp_path / "a", **package_kwargs([_quantiles(None)]))
    flat = TableOut("burden", 1, tuple(BURDEN_COLUMNS), [burden_row("TBX5")], "gene", "s1")
    with pytest.raises(PackageRefused, match="flat and gives a shard"):
        write_package(tmp_path / "b", **package_kwargs([flat]))
    with pytest.raises(PackageRefused):
        write_package(tmp_path / "c", **package_kwargs([_quantiles("..")]))
    with pytest.raises(PackageRefused):
        write_package(tmp_path / "d", **{**package_kwargs(), "version": ".."})
    assert [p.name for p in tmp_path.iterdir()] == ["ok"]


def test_every_row_of_a_source_column_must_name_the_package_exactly(tmp_path: Path) -> None:
    for bad in ["other", "demo ", None]:
        rows = [burden_row("TBX5"), burden_row("GATA4", source=bad)]
        with pytest.raises(PackageRefused, match="source"):
            write_package(tmp_path / "x", **package_kwargs([burden_table(rows)]))
    assert not (tmp_path / "x").exists()


def test_the_manifest_is_held_to_what_the_atlas_model_refuses(tmp_path: Path) -> None:
    kw = package_kwargs()
    dated = [{**kw["inputs"][0], "retrieved_on": d} for d in ("not-a-date", "2026-02-30")]
    for inputs in ([dated[0]], [dated[1]]):
        with pytest.raises(PackageRefused, match="manifest"):
            write_package(tmp_path / "d", **{**kw, "inputs": inputs})
    with pytest.raises(PackageRefused, match="repeated"):
        write_package(tmp_path / "e", **package_kwargs([burden_table(), burden_table()]))
    for blank in ["", "  \n"]:
        with pytest.raises(PackageRefused, match="blank"):
            write_package(tmp_path / "f", **{**kw, "licence_evidence": blank})
    assert not any(tmp_path.iterdir())


def test_an_unknown_table_and_unencodable_text_are_refused(tmp_path: Path) -> None:
    ghost = TableOut("nosuch", 1, ("a",), [["x"]], None)
    with pytest.raises(PackageRefused, match="no vendored schema"):
        write_package(tmp_path / "a", **package_kwargs([ghost]))
    rows = [burden_row("TBX5", method_note="bad\ud800")]
    with pytest.raises(PackageRefused, match="UTF-8"):
        write_package(tmp_path / "b", **package_kwargs([burden_table(rows)]))
    with pytest.raises(PackageRefused, match="UTF-8"):
        write_package(tmp_path / "c", **{**package_kwargs(), "licence_evidence": "x\ud800"})
    assert not any(tmp_path.iterdir())


def test_the_write_is_atomic_and_the_manifest_order_ignores_the_callers(tmp_path: Path) -> None:
    tables = [_quantiles("b"), burden_table()]
    write_package(tmp_path / "x", **package_kwargs(tables))
    write_package(tmp_path / "y", **package_kwargs(tables[::-1]))
    assert _bytes(tmp_path / "x") == _bytes(tmp_path / "y")
    listed = [
        (t["name"], t.get("shard"))
        for t in json.loads((tmp_path / "x/manifest.json").read_text())["tables"]
    ]
    assert listed == [("burden", None), ("profile_quantiles", "b")]
    assert sorted(p.name for p in tmp_path.iterdir()) == ["x", "y"]  # no staging directory left
    (tmp_path / "pre").mkdir()
    write_package(tmp_path / "pre", **package_kwargs())  # an existing empty directory is fine
    assert (tmp_path / "pre" / "manifest.json").is_file()


def test_the_version_pattern_equals_the_vendored_manifest_schema() -> None:
    from feeds.package import VENDORED
    from feeds.version import SEGMENT

    schema = json.loads((VENDORED / "package_manifest.schema.json").read_text())
    assert schema["properties"]["version"]["pattern"] == SEGMENT


def test_the_writer_refuses_a_wrong_gene_column_a_bad_manifest_and_a_used_directory(
    tmp_path: Path,
) -> None:
    wrong = TableOut("burden", 1, tuple(BURDEN_COLUMNS), [burden_row("TBX5")], None)
    with pytest.raises(PackageRefused, match="gene_column"):
        write_package(tmp_path / "a", **package_kwargs([wrong]))
    with pytest.raises(PackageRefused, match="manifest"):
        write_package(tmp_path / "b", **{**package_kwargs(), "ids": "bogus"})
    used = tmp_path / "c"
    used.mkdir()
    (used / "stale").write_text("x")
    with pytest.raises(PackageRefused, match="not empty"):
        write_package(used, **package_kwargs())
    assert not (tmp_path / "a").exists() and not (tmp_path / "b").exists()
    assert [p.name for p in used.iterdir()] == ["stale"]
