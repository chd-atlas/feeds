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
from datetime import date
from pathlib import Path

import pytest

from feeds.package import PackageRefused
from sources.clingen import recipe

BUILT_BY = {"host": "login01"}
RECIPE_COMMIT = "a" * 40
ATLAS_COMMIT = "4b1f8f68cbc30beaf932fbbc58e3617ebe401b6d"

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
    # Not the bootstrap constant's date, so a build that took MIRROR_BANNER would differ.
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


def _mirror(tmp_path: Path) -> Path:
    _, rows = recipe.read_rows(_csv(tmp_path))
    text = "\n".join("\t".join(r) for r in [tuple(COLUMNS), *rows]) + "\n"
    path = tmp_path / "clingen_gene_validity.tsv"
    path.write_text(text, encoding="utf-8")
    return path


def test_the_bootstrap_refuses_any_other_bytes_and_states_its_origin(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mirror = _mirror(tmp_path)
    kwargs = {"revision": 1, "built_by": BUILT_BY, "commit_of": lambda: RECIPE_COMMIT}

    # Another file is refused, unpatched: the real pinned hash is not this file's.
    with pytest.raises(PackageRefused, match="sha256"):
        recipe.bootstrap(mirror, atlas_commit=ATLAS_COMMIT, out=tmp_path / "a", **kwargs)  # type: ignore[arg-type]
    assert not (tmp_path / "a").exists()

    monkeypatch.setattr(recipe, "MIRROR_SHA256", hashlib.sha256(mirror.read_bytes()).hexdigest())
    with pytest.raises(PackageRefused, match="pinned"):
        recipe.bootstrap(mirror, atlas_commit="4b1f8f6", out=tmp_path / "b", **kwargs)  # type: ignore[arg-type]
    # A well-formed 40-hex commit that is not the pinned one is refused too.
    with pytest.raises(PackageRefused, match="pinned"):
        recipe.bootstrap(mirror, atlas_commit="b" * 40, out=tmp_path / "b2", **kwargs)  # type: ignore[arg-type]

    out = recipe.bootstrap(mirror, atlas_commit=ATLAS_COMMIT, out=tmp_path / "c", **kwargs)  # type: ignore[arg-type]
    manifest, table = _read_package(out)
    assert manifest["version"] == "2026-08-03"
    assert manifest["source"]["version"] == "FILE CREATED: 2026-08-03"  # type: ignore[index]
    assert manifest["source"]["retrieved_on"] == "2026-08-03"  # type: ignore[index]
    (entry,) = manifest["inputs"]  # type: ignore[misc]
    assert entry == {
        "name": "clingen_gene_validity.tsv",
        "url": "https://raw.githubusercontent.com/chd-atlas/chd-atlas.github.io/"
        f"{ATLAS_COMMIT}/mirrors/clingen_gene_validity.tsv",
        "sha256": hashlib.sha256(mirror.read_bytes()).hexdigest(),
        "retrieved_on": "2026-08-03",
        "upstream_release": "FILE CREATED: 2026-08-03",
    }
    assert manifest["input_counts"] == {"rows": 4, "genes": 2}
    assert table[0] == COLUMNS and len(table) == 5

    # A mirror whose header is not the table's column order is refused.
    bad = tmp_path / "bad.tsv"
    bad.write_text(mirror.read_text().replace("report_url", "url"), encoding="utf-8")
    monkeypatch.setattr(recipe, "MIRROR_SHA256", hashlib.sha256(bad.read_bytes()).hexdigest())
    with pytest.raises(PackageRefused, match="header"):
        recipe.bootstrap(bad, atlas_commit=ATLAS_COMMIT, out=tmp_path / "d", **kwargs)  # type: ignore[arg-type]
