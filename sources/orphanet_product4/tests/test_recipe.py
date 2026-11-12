"""The Orphanet product 4 recipe on synthetic XML: no network, no committed file.

The package carries every HPO annotation (owner decision F23); the atlas filters at
ingest, so Macrocephaly below is a row.
"""

import gzip
import hashlib
import json
from pathlib import Path

import pytest

from feeds.orphadata import jdbor_release
from feeds.package import PackageRefused
from sources.orphanet_product4 import recipe

BUILT_BY = {"host": "login01"}
RECIPE_COMMIT = "a" * 40

_HEAD = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<JDBOR date="2030-01-02 03:04:05" version="9.9.9 / test [2030-01-01]" copyright="x">\n'
    '<HPODisorderSetStatusList count="2">\n'
)
_TAIL = "</HPODisorderSetStatusList>\n</JDBOR>\n"

NOONAN = """
<HPODisorderSetStatus id="1">
  <Disorder id="2">
    <OrphaCode>648</OrphaCode>
    <Name lang="en">Noonan syndrome</Name>
    <HPODisorderAssociationList count="3">
      <HPODisorderAssociation id="10">
        <HPO id="4"><HPOId>HP:0001642</HPOId><HPOTerm>Pulmonic stenosis</HPOTerm></HPO>
        <HPOFrequency id="5"><Name lang="en">Frequent (79-30%)</Name></HPOFrequency>
      </HPODisorderAssociation>
      <HPODisorderAssociation id="11">
        <HPO id="6"><HPOId>HP:0000256</HPOId><HPOTerm>Macrocephaly</HPOTerm></HPO>
        <HPOFrequency id="7"><Name lang="en">Very frequent (99-80%)</Name></HPOFrequency>
      </HPODisorderAssociation>
      <HPODisorderAssociation id="12">
        <HPO id="8"><HPOId>HP:0001631</HPOId><HPOTerm>Atrial septal defect</HPOTerm></HPO>
        <HPOFrequency id="9"/>
      </HPODisorderAssociation>
    </HPODisorderAssociationList>
  </Disorder>
</HPODisorderSetStatus>
"""


def _assoc(hpo: str, term: str, band: str) -> str:
    freq = f'<HPOFrequency id="1"><Name lang="en">{band}</Name></HPOFrequency>' if band else ""
    return (
        f'<HPODisorderAssociation id="1"><HPO id="1"><HPOId>{hpo}</HPOId>'
        f"<HPOTerm>{term}</HPOTerm></HPO>{freq}</HPODisorderAssociation>"
    )


def _disorder(orpha: str, name: str, *assocs: str) -> str:
    return (
        '<HPODisorderSetStatus id="1"><Disorder id="1">'
        f'<OrphaCode>{orpha}</OrphaCode><Name lang="en">{name}</Name>'
        f'<HPODisorderAssociationList count="{len(assocs)}">{"".join(assocs)}'
        "</HPODisorderAssociationList></Disorder></HPODisorderSetStatus>\n"
    )


def _xml(tmp_path: Path, *disorders: str, name: str = "en_product4.xml") -> Path:
    path = tmp_path / name
    path.write_text(_HEAD + "".join(disorders) + _TAIL, encoding="utf-8")
    return path


@pytest.fixture
def pinned(monkeypatch: pytest.MonkeyPatch) -> object:
    """Pin whatever bytes a test writes, so ``build`` accepts a synthetic file."""

    def _pin(path: Path) -> None:
        monkeypatch.setattr(recipe, "PINNED_SHA256", hashlib.sha256(path.read_bytes()).hexdigest())

    return _pin


def _build(path: Path, out: Path) -> Path:
    return recipe.build(
        path, revision=1, built_by=BUILT_BY, out=out, commit_of=lambda: RECIPE_COMMIT
    )


def _read(root: Path) -> tuple[dict[str, object], list[list[str]]]:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    raw = gzip.decompress((root / "orphanet_phenotypes.tsv.gz").read_bytes()).decode("utf-8")
    return manifest, [line.split("\t") for line in raw.splitlines()]


def test_a_phenotype_row_keeps_orphanets_band_verbatim_for_every_annotation(
    tmp_path: Path,
) -> None:
    """Macrocephaly is kept (the atlas filters); an empty band stays ''; rows are keyed."""
    rows, dropped = recipe.read_rows(_xml(tmp_path, NOONAN))
    assert dropped == []
    assert rows == [
        (
            "ORPHA:648",
            "Noonan syndrome",
            "HP:0000256",
            "Macrocephaly",
            "Very frequent (99-80%)",
            "orphanet_product4",
        ),
        (
            "ORPHA:648",
            "Noonan syndrome",
            "HP:0001631",
            "Atrial septal defect",
            "",
            "orphanet_product4",
        ),
        (
            "ORPHA:648",
            "Noonan syndrome",
            "HP:0001642",
            "Pulmonic stenosis",
            "Frequent (79-30%)",
            "orphanet_product4",
        ),
    ]


def test_the_sort_key_is_orpha_hpo_then_frequency(tmp_path: Path) -> None:
    xml = _xml(
        tmp_path,
        _disorder("20", "B", _assoc("HP:0000002", "t", "Occasional (29-5%)")),
        _disorder(
            "100",
            "A",
            _assoc("HP:0000009", "t", "Obligate (100%)"),
            _assoc("HP:0000001", "t", "Very rare (&lt;4-1%)"),
            _assoc("HP:0000001", "t", "Frequent (79-30%)"),
        ),
    )
    rows, _ = recipe.read_rows(xml)
    # ORPHA:100 < ORPHA:20 as text; the same (orpha, hpo) is ordered by its band.
    assert [(r[0], r[2], r[4]) for r in rows] == [
        ("ORPHA:100", "HP:0000001", "Frequent (79-30%)"),
        ("ORPHA:100", "HP:0000001", "Very rare (<4-1%)"),
        ("ORPHA:100", "HP:0000009", "Obligate (100%)"),
        ("ORPHA:20", "HP:0000002", "Occasional (29-5%)"),
    ]


def test_exact_duplicate_annotations_collapse_into_one_dropped_entry(
    tmp_path: Path, pinned: object
) -> None:
    """Identical rows collapse; a differing band keeps both; the count is one entry."""
    dup = _assoc("HP:0001642", "Pulmonic stenosis", "Frequent (79-30%)")
    xml = _xml(
        tmp_path,
        _disorder(
            "648",
            "Noonan syndrome",
            dup,
            dup,
            dup,
            _assoc("HP:0001642", "Pulmonic stenosis", "Occasional (29-5%)"),
        ),
        _disorder("7", "Other", _assoc("HP:0000001", "t", ""), _assoc("HP:0000001", "t", "")),
        _disorder("9", "Unique", _assoc("HP:0000005", "t", "")),
    )
    rows, dropped = recipe.read_rows(xml)
    assert len(rows) == 4, "one Frequent, one Occasional, one empty band, one unique pair"
    assert len(dropped) == 1
    assert dropped[0]["count"] == 3
    assert dropped[0]["examples"] == ["ORPHA:648 HP:0001642", "ORPHA:7 HP:0000001"]
    assert dropped[0]["reason"] == (
        "annotation identical to another in every table column in en_product4.xml"
    )

    pinned(xml)  # type: ignore[operator]
    manifest, table = _read(_build(xml, tmp_path / "pkg"))
    # Decision 4: every annotation read, kept or collapsed.
    assert manifest["input_counts"] == {"rows": 7}
    assert manifest["dropped"] == dropped
    assert len(table) == 1 + 4


def test_at_most_five_distinct_examples_are_named(tmp_path: Path) -> None:
    a = _assoc("HP:0000001", "t", "")
    xml = _xml(tmp_path, *(_disorder(str(n), "d", a, a) for n in range(1, 9)))
    _, dropped = recipe.read_rows(xml)
    assert dropped[0]["count"] == 8
    assert dropped[0]["examples"] == [f"ORPHA:{n} HP:0000001" for n in range(1, 6)]


def test_the_source_block_carries_the_files_own_header(tmp_path: Path, pinned: object) -> None:
    xml = _xml(tmp_path, NOONAN)
    pinned(xml)  # type: ignore[operator]
    out = _build(xml, tmp_path / "pkg")
    manifest, table = _read(out)
    sha = hashlib.sha256(xml.read_bytes()).hexdigest()
    # The synthetic header differs from the recorded constants, so the manifest can only
    # have taken it from the file.
    assert jdbor_release(xml) == ("9.9.9 / test [2030-01-01]", "2030-01-02 03:04:05")
    assert manifest["version"] == sha[:8]
    assert manifest["ids"] == "hgnc_id" and manifest["coverage"] == "every_gene"
    src = manifest["source"]
    assert src["version"] == f"sha256:{sha}"  # type: ignore[index]
    assert src["file_version"] == "9.9.9 / test [2030-01-01]"  # type: ignore[index]
    assert src["file_date"] == "2030-01-02 03:04:05"  # type: ignore[index]
    assert src["retrieved_on"] == "2026-10-25"  # type: ignore[index]
    assert src["url"] == "https://www.orphadata.com/data/xml/en_product4.xml"  # type: ignore[index]
    assert src["licence"] == "CC-BY-4.0"  # type: ignore[index]
    assert src["redistribution"] == "permitted_with_attribution"  # type: ignore[index]
    (entry,) = manifest["inputs"]  # type: ignore[misc]
    assert entry["sha256"] == sha and entry["retrieved_on"] == "2026-10-25"
    assert entry["upstream_release"] == "9.9.9 / test [2030-01-01]"
    assert manifest["input_counts"] == {"rows": 3}
    assert table[0] == ["orpha", "disorder_name", "hpo", "hpo_label", "frequency", "source"]


def test_bytes_other_than_the_pinned_download_are_refused(tmp_path: Path) -> None:
    xml = _xml(tmp_path, NOONAN)
    with pytest.raises(PackageRefused, match="not the pinned"):
        _build(xml, tmp_path / "pkg")
    assert not (tmp_path / "pkg").exists()


def test_a_local_download_states_the_recorded_release() -> None:
    """Where the download is present -- a curator's machine, never CI -- it must agree."""
    path = Path(__file__).resolve().parents[3] / "local" / "orphanet" / "en_product4.xml"
    if not path.is_file():
        pytest.skip(f"{path} is not present (gitignored; a curator's copy)")
    assert jdbor_release(path) == (recipe.FILE_VERSION, recipe.FILE_DATE)
