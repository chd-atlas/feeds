"""The Orphanet product 1 recipe on synthetic XML: no network, no committed file."""

import gzip
import hashlib
import json
from pathlib import Path

import pytest

from feeds.orphadata import jdbor_release
from feeds.package import PackageRefused
from sources.orphanet_product1 import recipe

BUILT_BY = {"host": "login01"}
RECIPE_COMMIT = "a" * 40

_HEAD = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<JDBOR date="2026-06-23 07:53:50" version="1.3.42 / 4.1.8 [2025-03-03]" copyright="x">\n'
    '<DisorderList count="2">\n'
)
_TAIL = "</DisorderList>\n</JDBOR>\n"

NOONAN = """
<Disorder id="17601">
  <OrphaCode>648</OrphaCode>
  <Name lang="en">Noonan syndrome</Name>
  <ExternalReferenceList count="3">
    <ExternalReference id="1">
      <Source>ICD-10</Source>
      <Reference>Q87.1</Reference>
      <DisorderMappingRelation id="2">
        <Name lang="en">NTBT (narrower)</Name>
      </DisorderMappingRelation>
      <DisorderMappingValidationStatus id="3">
        <Name lang="en">Validated</Name>
      </DisorderMappingValidationStatus>
    </ExternalReference>
    <ExternalReference id="1">
      <Source>ORPHA-ICD11</Source>
      <Reference>Q87.1</Reference>
      <DisorderMappingRelation id="2">
        <Name lang="en">NTBT (narrower)</Name>
      </DisorderMappingRelation>
      <DisorderMappingValidationStatus id="3">
        <Name lang="en">Validated</Name>
      </DisorderMappingValidationStatus>
    </ExternalReference>
    <ExternalReference id="4">
      <Source>OMIM</Source>
      <Reference>610733</Reference>
      <DisorderMappingRelation id="5">
        <Name lang="en">BTNT (ORPHAcode is broader than the targeted code ...)</Name>
      </DisorderMappingRelation>
      <DisorderMappingValidationStatus id="6">
        <Name lang="en">Validated</Name>
      </DisorderMappingValidationStatus>
    </ExternalReference>
  </ExternalReferenceList>
</Disorder>
"""


def _disorder(orpha: str, *omims: tuple[str, str, str]) -> str:
    refs = "".join(
        f"<ExternalReference id='1'><Source>OMIM</Source><Reference>{omim}</Reference>"
        f"<DisorderMappingRelation id='1'><Name lang='en'>{rel}</Name></DisorderMappingRelation>"
        f"<DisorderMappingValidationStatus id='1'><Name lang='en'>{status}</Name>"
        "</DisorderMappingValidationStatus></ExternalReference>"
        for omim, rel, status in omims
    )
    return (
        f'<Disorder id="1"><OrphaCode>{orpha}</OrphaCode><Name lang="en">d</Name>'
        f'<ExternalReferenceList count="1">{refs}</ExternalReferenceList></Disorder>\n'
    )


def _xml(tmp_path: Path, *disorders: str) -> Path:
    path = tmp_path / "en_product1.xml"
    path.write_text(_HEAD + "".join(disorders) + _TAIL, encoding="utf-8")
    return path


def _build(path: Path, out: Path) -> Path:
    return recipe.build(
        path, revision=1, built_by=BUILT_BY, out=out, commit_of=lambda: RECIPE_COMMIT
    )


def test_an_omim_row_takes_the_relation_code_and_skips_other_sources(tmp_path: Path) -> None:
    rows = recipe.read_rows(_xml(tmp_path, NOONAN))
    assert rows == [("ORPHA:648", "OMIM:610733", "BTNT", "Validated", "orphanet_product1")]


def test_every_relation_and_status_is_kept_and_rows_sort_by_orpha_then_omim(
    tmp_path: Path,
) -> None:
    xml = _xml(
        tmp_path,
        _disorder("20", ("100", "E (exact)", "Validated")),
        _disorder(
            "100",
            ("300", "NTBT (narrower)", "Not yet validated"),
            ("200", "BTNT (broader)", "Validated"),
        ),
    )
    rows = recipe.read_rows(xml)
    assert [(r[0], r[1], r[2], r[3]) for r in rows] == [
        ("ORPHA:100", "OMIM:200", "BTNT", "Validated"),
        ("ORPHA:100", "OMIM:300", "NTBT", "Not yet validated"),
        ("ORPHA:20", "OMIM:100", "E", "Validated"),
    ]


def test_the_source_block_carries_the_files_own_header(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    xml = _xml(tmp_path, NOONAN)
    monkeypatch.setattr(recipe, "PINNED_SHA256", hashlib.sha256(xml.read_bytes()).hexdigest())
    out = _build(xml, tmp_path / "pkg")
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    table = gzip.decompress((out / "orphanet_omim.tsv.gz").read_bytes()).decode().splitlines()
    sha = hashlib.sha256(xml.read_bytes()).hexdigest()
    assert manifest["version"] == sha[:8]
    assert manifest["ids"] == "hgnc_id" and manifest["coverage"] == "every_gene"
    assert manifest["source"] == {
        "name": manifest["source"]["name"],
        "version": f"sha256:{sha}",
        "file_version": "1.3.42 / 4.1.8 [2025-03-03]",
        "file_date": "2026-06-23 07:53:50",
        "retrieved_on": "2026-10-25",
        "url": "https://www.orphadata.com/data/xml/en_product1.xml",
        "licence": "CC-BY-4.0",
        "redistribution": "permitted_with_attribution",
    }
    (entry,) = manifest["inputs"]
    assert entry["sha256"] == sha and entry["upstream_release"] == "1.3.42 / 4.1.8 [2025-03-03]"
    assert manifest["input_counts"] == {"rows": 1} and manifest["dropped"] == []
    assert table == [
        "orpha\tomim\trelation\tvalidation\tsource",
        "ORPHA:648\tOMIM:610733\tBTNT\tValidated\torphanet_product1",
    ]


def test_bytes_other_than_the_pinned_download_are_refused(tmp_path: Path) -> None:
    with pytest.raises(PackageRefused, match="not the pinned"):
        _build(_xml(tmp_path, NOONAN), tmp_path / "pkg")
    assert not (tmp_path / "pkg").exists()


def test_a_local_download_states_the_recorded_release() -> None:
    """Where the download is present -- a curator's machine, never CI -- it must agree."""
    path = Path(__file__).resolve().parents[3] / "local" / "orphanet" / "en_product1.xml"
    if not path.is_file():
        pytest.skip(f"{path} is not present (gitignored; a curator's copy)")
    assert jdbor_release(path) == (recipe.FILE_VERSION, recipe.FILE_DATE)
