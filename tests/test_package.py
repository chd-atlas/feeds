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


def test_a_shard_names_its_file_and_a_bad_shard_or_version_is_refused(tmp_path: Path) -> None:
    kw = package_kwargs(
        [TableOut("burden", 1, tuple(BURDEN_COLUMNS), [burden_row("TBX5")], "gene", "s1")]
    )
    write_package(tmp_path / "ok", **kw)
    assert (tmp_path / "ok" / "burden.s1.tsv.gz").is_file()
    bad = TableOut("burden", 1, tuple(BURDEN_COLUMNS), [burden_row("TBX5")], "gene", "..")
    with pytest.raises(PackageRefused):
        write_package(tmp_path / "x", **package_kwargs([bad]))
    with pytest.raises(PackageRefused):
        write_package(tmp_path / "y", **{**package_kwargs(), "version": ".."})
    assert not (tmp_path / "x").exists() and not (tmp_path / "y").exists()


def test_the_writer_refuses_a_wrong_gene_column_a_bad_manifest_and_a_used_directory(
    tmp_path: Path,
) -> None:
    wrong = TableOut("burden", 1, tuple(BURDEN_COLUMNS), [burden_row("TBX5")], None)
    with pytest.raises(PackageRefused, match="gene_column"):
        write_package(tmp_path / "a", **package_kwargs([wrong]))
    with pytest.raises(PackageRefused, match="manifest"):
        write_package(tmp_path / "b", **{**package_kwargs(), "licence_evidence": ""})
    used = tmp_path / "c"
    used.mkdir()
    (used / "stale").write_text("x")
    with pytest.raises(PackageRefused, match="not empty"):
        write_package(used, **package_kwargs())
    assert not (tmp_path / "a").exists() and not (tmp_path / "b").exists()
    assert [p.name for p in used.iterdir()] == ["stale"]
