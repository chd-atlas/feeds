"""The GenCC recipe on synthetic inputs: no network, no committed mirror.

The real export is quoted TSV: measured 2026-11-11, 46,968 physical lines hold 30,376
records because uncarried free-text cells contain quoted newlines. The fixture
reproduces that with one cell.
"""

import gzip
import hashlib
import json
from datetime import date
from pathlib import Path

import pytest

from feeds.package import PackageRefused
from sources.gencc import recipe

BUILT_BY = {"host": "login01"}
RECIPE_COMMIT = "a" * 40
ATLAS_COMMIT = "b" * 40

COLUMNS = [
    "sgc_id",
    "gene",
    "gene_symbol",
    "disease",
    "disease_label",
    "moi",
    "classification",
    "submitter",
    "submitted_on",
    "report_url",
]

# Upstream names, in an order that is not the table's, plus two uncarried columns.
HEADER = [
    "uuid",
    "gene_curie",
    "gene_symbol",
    "disease_curie",
    "disease_title",
    "moi_title",
    "classification_title",
    "submitter_title",
    "submitted_as_date",
    "submitted_as_public_report_url",
    "sgc_id",
    "assertion_criteria_url",
]


def _line(sgc: str, gene: str, disease: str, submitter: str, note: str = "x", **over: str) -> str:
    cells = {
        "uuid": "u-" + sgc,
        "gene_curie": gene,
        "gene_symbol": "SYM",
        "disease_curie": disease,
        "disease_title": "a disease",
        "moi_title": "Autosomal dominant",
        "classification_title": "Strong",
        "submitter_title": submitter,
        "submitted_as_date": "2018-05-31 20:04:32",
        "submitted_as_public_report_url": "",
        "sgc_id": sgc,
        "assertion_criteria_url": note,
    }
    cells.update(over)
    return "\t".join(cells[h] for h in HEADER)


def _tsv(tmp_path: Path, *lines: str, name: str = "gencc.tsv") -> Path:
    path = tmp_path / name
    path.write_text("\t".join(HEADER) + "\n" + "\n".join(lines) + "\n", encoding="utf-8")
    return path


GOOD = (
    _line("SGC-3", "HGNC:4173", "MONDO:0005453", "B lab"),
    _line("SGC-2", "HGNC:11604", "MONDO:0007732", "B lab"),
    _line("SGC-1", "HGNC:11604", "MONDO:0007732", "A lab", classification_title="Limited"),
)


def _read_package(root: Path) -> tuple[dict[str, object], list[list[str]]]:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    raw = gzip.decompress((root / "gencc_submissions.tsv.gz").read_bytes()).decode("utf-8")
    return manifest, [line.split("\t") for line in raw.splitlines()]


def test_a_quoted_free_text_cell_spanning_lines_is_one_record(tmp_path: Path) -> None:
    path = _tsv(
        tmp_path,
        _line("SGC-1", "HGNC:1", "MONDO:0000001", "A lab", note='"a\nb"'),
        _line("SGC-2", "HGNC:2", "MONDO:0000002", "A lab"),
    )
    assert len(path.read_text().splitlines()) == 4, "header + 3 physical lines"
    rows, dropped = recipe.read_rows(path)
    assert [r[0] for r in rows] == ["SGC-1", "SGC-2"]
    assert dropped == []


def test_rows_follow_table_order_and_the_key_and_the_package_states_its_origin(
    tmp_path: Path,
) -> None:
    rows, dropped = recipe.read_rows(_tsv(tmp_path, *GOOD))
    # (gene, disease, moi, submitter, sgc_id): "HGNC:11604" sorts before "HGNC:4173" as text.
    assert [(r[1], r[7], r[0]) for r in rows] == [
        ("HGNC:11604", "A lab", "SGC-1"),
        ("HGNC:11604", "B lab", "SGC-2"),
        ("HGNC:4173", "B lab", "SGC-3"),
    ]
    assert rows[0] == (
        "SGC-1",
        "HGNC:11604",
        "SYM",
        "MONDO:0007732",
        "a disease",
        "Autosomal dominant",
        "Limited",
        "A lab",
        "2018-05-31 20:04:32",
        "",
    )
    assert dropped == []

    path = _tsv(tmp_path, *GOOD)
    out = recipe.build(
        path,
        retrieved_on=date(2026, 11, 11),
        revision=1,
        built_by=BUILT_BY,
        out=tmp_path / "pkg",
        commit_of=lambda: RECIPE_COMMIT,
    )
    manifest, table = _read_package(out)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    assert manifest["version"] == sha[:8]
    assert manifest["ids"] == "hgnc_id"
    assert manifest["source"]["version"] == "retrieved 2026-11-11"  # type: ignore[index]
    assert manifest["input_counts"] == {"rows": 3, "genes": 2}
    assert manifest["dropped"] == []
    (entry,) = manifest["inputs"]  # type: ignore[misc]
    assert entry["sha256"] == sha and entry["upstream_release"] is None
    assert table[0] == COLUMNS and len(table) == 4


def test_a_row_outside_the_curie_forms_is_dropped_and_counted(tmp_path: Path) -> None:
    path = _tsv(
        tmp_path,
        *GOOD,
        _line("SGC-8", "ENSG1", "MONDO:0000001", "A lab"),
        _line("SGC-9", "HGNC:5", "OMIM:123456", "A lab"),
    )
    rows, dropped = recipe.read_rows(path)
    assert len(rows) == 3
    assert [(d["count"], d["examples"]) for d in dropped] == [(2, ["SGC-8", "SGC-9"])]
    out = recipe.build(
        path,
        retrieved_on=date(2026, 11, 11),
        revision=1,
        built_by=BUILT_BY,
        out=tmp_path / "pkg",
        commit_of=lambda: RECIPE_COMMIT,
    )
    manifest, _ = _read_package(out)
    assert manifest["dropped"][0]["count"] == 2  # type: ignore[index]
    assert manifest["input_counts"] == {"rows": 3, "genes": 2}


def test_a_carried_cell_with_a_quote_refuses_naming_the_record(tmp_path: Path) -> None:
    path = _tsv(
        tmp_path,
        *GOOD,
        _line("SGC-77", "HGNC:5", "MONDO:0000005", "A lab", disease_title='"""quoted"" disease"'),
    )
    with pytest.raises(PackageRefused, match="SGC-77"):
        recipe.read_rows(path)


def test_a_renamed_upstream_column_stops_the_build(tmp_path: Path) -> None:
    path = _tsv(tmp_path, *GOOD)
    path.write_text(path.read_text().replace("moi_title", "moi"), encoding="utf-8")
    with pytest.raises(PackageRefused, match="moi_title"):
        recipe.read_rows(path)


def _mirror(tmp_path: Path) -> Path:
    rows, _ = recipe.read_rows(_tsv(tmp_path, *GOOD))
    text = "\n".join("\t".join(r) for r in [tuple(COLUMNS), *rows]) + "\n"
    path = tmp_path / "gencc_submissions.tsv"
    path.write_text(text, encoding="utf-8")
    return path


def test_the_bootstrap_refuses_any_other_bytes_and_states_its_origin(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mirror = _mirror(tmp_path)
    kwargs = {"revision": 1, "built_by": BUILT_BY, "commit_of": lambda: RECIPE_COMMIT}

    with pytest.raises(PackageRefused, match="sha256"):
        recipe.bootstrap(mirror, atlas_commit=ATLAS_COMMIT, out=tmp_path / "a", **kwargs)  # type: ignore[arg-type]
    assert not (tmp_path / "a").exists()

    sha = hashlib.sha256(mirror.read_bytes()).hexdigest()
    monkeypatch.setattr(recipe, "MIRROR_SHA256", sha)
    with pytest.raises(PackageRefused, match="40"):
        recipe.bootstrap(mirror, atlas_commit="4b1f8f6", out=tmp_path / "b", **kwargs)  # type: ignore[arg-type]

    out = recipe.bootstrap(mirror, atlas_commit=ATLAS_COMMIT, out=tmp_path / "c", **kwargs)  # type: ignore[arg-type]
    manifest, table = _read_package(out)
    assert manifest["version"] == sha[:8]
    assert manifest["source"]["version"] == "retrieved 2026-08-03"  # type: ignore[index]
    assert manifest["source"]["retrieved_on"] == "2026-08-03"  # type: ignore[index]
    (entry,) = manifest["inputs"]  # type: ignore[misc]
    assert entry == {
        "name": "gencc_submissions.tsv",
        "url": "https://raw.githubusercontent.com/chd-atlas/chd-atlas.github.io/"
        f"{ATLAS_COMMIT}/mirrors/gencc_submissions.tsv",
        "sha256": sha,
        "retrieved_on": "2026-08-03",
        "upstream_release": None,
    }
    assert manifest["input_counts"] == {"rows": 3, "genes": 2}
    assert table[0] == COLUMNS and len(table) == 4

    bad = tmp_path / "bad.tsv"
    bad.write_text(mirror.read_text().replace("report_url", "url"), encoding="utf-8")
    monkeypatch.setattr(recipe, "MIRROR_SHA256", hashlib.sha256(bad.read_bytes()).hexdigest())
    with pytest.raises(PackageRefused, match="header"):
        recipe.bootstrap(bad, atlas_commit=ATLAS_COMMIT, out=tmp_path / "d", **kwargs)  # type: ignore[arg-type]
