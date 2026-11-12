"""The Mondo SSSOM recipe on synthetic files: no network, no committed file."""

import gzip
import hashlib
import json
from pathlib import Path

import pytest

from feeds.package import PackageRefused
from sources.mondo_sssom import recipe

BUILT_BY = {"host": "login01"}
RECIPE_COMMIT = "a" * 40
HEADER = "subject_id\tpredicate_id\tobject_id\tmapping_justification"
COMMENTS = "# curie_map:\n#   OMIM: https://omim.org/entry/\n# mapping_set_id: mondo\n"


def _sssom(tmp_path: Path, *triples: tuple[str, str, str], comments: str = COMMENTS) -> Path:
    path = tmp_path / "mondo.sssom.tsv"
    body = "".join(f"{s}\t{p}\t{o}\tsemapv:ManualMappingCuration\n" for s, p, o in triples)
    path.write_text(comments + HEADER + "\n" + body, encoding="utf-8")
    return path


def test_an_xref_is_an_exact_match_to_orphanet_or_omim_for_any_mondo_term(tmp_path: Path) -> None:
    """MONDO spells Orphanet `Orphanet:`; the table spells it `ORPHA:`. No term filter."""
    path = _sssom(
        tmp_path,
        ("MONDO:0018997", "skos:exactMatch", "Orphanet:648"),
        ("MONDO:0018997", "skos:exactMatch", "OMIM:163950"),
        ("MONDO:0018997", "skos:closeMatch", "Orphanet:98733"),  # not exact
        ("MONDO:0018997", "skos:exactMatch", "DOID:3490"),  # not a target the rule follows
        ("MONDO:0018997", "skos:exactMatch", "OMIMPS:163950"),  # a phenotypic series
        ("DOID:1", "skos:exactMatch", "Orphanet:2"),  # not a MONDO subject
        ("MONDO:0999999", "skos:exactMatch", "Orphanet:1"),  # any term, mirrored or not
    )
    rows, count = recipe.read_rows(path)
    assert rows == [
        ("MONDO:0018997", "OMIM:163950", "mondo_sssom"),
        ("MONDO:0018997", "ORPHA:648", "mondo_sssom"),
        ("MONDO:0999999", "ORPHA:1", "mondo_sssom"),
    ]
    assert count == 3


def test_comment_lines_are_skipped_and_orphanet_is_normalised(tmp_path: Path) -> None:
    path = _sssom(
        tmp_path,
        ("MONDO:0000002", "skos:exactMatch", "Orphanet:20"),
        ("MONDO:0000001", "skos:exactMatch", "OMIM:100"),
        ("MONDO:0000001", "skos:exactMatch", "OMIM:100"),  # a repeated pair is one row
        ("MONDO:0000002", "skos:exactMatch", "Orphanet:3"),
    )
    rows, _ = recipe.read_rows(path)
    assert rows == [
        ("MONDO:0000001", "OMIM:100", "mondo_sssom"),
        ("MONDO:0000002", "ORPHA:20", "mondo_sssom"),
        ("MONDO:0000002", "ORPHA:3", "mondo_sssom"),
    ]


def test_a_carried_cell_the_table_cannot_hold_is_refused(tmp_path: Path) -> None:
    path = _sssom(tmp_path, ("MONDO:0000001", "skos:exactMatch", 'OMIM:1"2'))
    with pytest.raises(PackageRefused, match="MONDO:0000001"):
        recipe.read_rows(path)


def test_the_package_states_the_pinned_release(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _sssom(
        tmp_path,
        ("MONDO:0000001", "skos:exactMatch", "Orphanet:648"),
        ("MONDO:0000001", "skos:broadMatch", "Orphanet:9"),
    )
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    monkeypatch.setattr(recipe, "PINNED_SHA256", sha)
    out = recipe.build(
        path, revision=1, built_by=BUILT_BY, out=tmp_path / "pkg", commit_of=lambda: RECIPE_COMMIT
    )
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    table = gzip.decompress((out / "mondo_xrefs.tsv.gz").read_bytes()).decode().splitlines()
    assert manifest["version"] == "2026-10-06"
    assert manifest["ids"] == "hgnc_id" and manifest["coverage"] == "every_gene"
    assert manifest["source"]["version"] == f"sha256:{sha}"
    assert manifest["source"]["retrieved_on"] == "2026-10-25"
    assert manifest["source"]["licence"] == "CC-BY-4.0"
    assert manifest["source"]["redistribution"] == "permitted_with_attribution"
    assert manifest["source"]["url"] == (
        "https://raw.githubusercontent.com/monarch-initiative/mondo/v2026-10-06"
        "/src/ontology/mappings/mondo.sssom.tsv"
    )
    (entry,) = manifest["inputs"]
    assert entry["name"] == "mondo.sssom.tsv" and entry["sha256"] == sha
    assert entry["url"] == manifest["source"]["url"]
    assert entry["upstream_release"] == "v2026-10-06"
    assert manifest["input_counts"] == {"rows": 1} and manifest["dropped"] == []
    assert table == ["mondo\ttarget\tsource", "MONDO:0000001\tORPHA:648\tmondo_sssom"]


def test_bytes_other_than_the_pinned_download_are_refused(tmp_path: Path) -> None:
    path = _sssom(tmp_path, ("MONDO:0000001", "skos:exactMatch", "Orphanet:648"))
    with pytest.raises(PackageRefused, match="not the pinned"):
        recipe.build(
            path,
            revision=1,
            built_by=BUILT_BY,
            out=tmp_path / "pkg",
            commit_of=lambda: RECIPE_COMMIT,
        )
    assert not (tmp_path / "pkg").exists()


def test_a_local_download_is_the_pinned_one() -> None:
    """Where the download is present -- a curator's machine, never CI -- it must agree."""
    path = Path(__file__).resolve().parents[3] / "local" / "mondo" / "mondo.sssom.tsv"
    if not path.is_file():
        pytest.skip(f"{path} is not present (gitignored; a curator's copy)")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == recipe.PINNED_SHA256
