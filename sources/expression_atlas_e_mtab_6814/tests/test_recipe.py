"""The E-MTAB-6814 recipe on synthetic records: no network, no committed mirror.

Ported from the atlas's ``tests/unit/test_convert_cardoso_moreira.py`` (all 15 cases,
which call the builders directly), plus the cases for what changed when the recipe
moved: every matrix gene is carried under its Ensembl id (no registry join), an
unlabelled assay group refuses instead of being skipped, and the percentile grid is
computed over every matrix gene.

The recipe refuses any input whose sha256 is not the pinned one, so the build tests
patch the pins to the synthetic files' hashes; ``test_other_bytes_are_refused`` runs
unpatched.
"""

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from feeds.package import PackageRefused
from feeds.version import package_version
from sources.expression_atlas_e_mtab_6814 import recipe
from sources.expression_atlas_e_mtab_6814.recipe import (
    UNIT,
    AssayGroup,
    parse_cell,
    profile_row,
    quantile_linear,
    quantile_rows,
)

BUILT_BY = {"host": "workstation"}
RECIPE_COMMIT = "a" * 40


@pytest.mark.parametrize(
    ("cell", "expected"),
    [
        ("8,15,18,20,21", (15.0, 18.0, 20.0)),
        ("13,13,13,13,13", (13.0, 13.0, 13.0)),
        ("7", (7.0, 7.0, 7.0)),
    ],
)
def test_a_quintuple_yields_q25_median_q75(cell: str, expected: tuple[float, float, float]) -> None:
    """The Atlas hands over the quartiles; nothing here recomputes them.

    A single value is the degenerate quintuple the Atlas emits for a one-assay
    group, and is accepted as all five being equal.
    """
    assert parse_cell(cell) == expected


@pytest.mark.parametrize("cell", ["", "   ", "NA", "nan", "1,2,3"])
def test_an_unquantified_cell_is_none_and_never_zero(cell: str) -> None:
    """ "Not measured" and "measured as zero" are different facts.

    Returning 0.0 would put a fabricated floor-level value into the quantile grid
    and give the gene a ``median_abundance`` the source never reported. A malformed
    cell (``1,2,3``: neither one value nor five) is refused for the same reason.
    """
    assert parse_cell(cell) is None


def test_quartiles_are_dropped_below_three_replicates_but_the_median_survives() -> None:
    """Quartiles of two points are not quartiles; the median still publishes."""
    three = AssayGroup("g1", "7 week post conception", "heart", 3)
    two = AssayGroup("g2", "4 week post conception", "heart", 2)

    rich = profile_row("E-MTAB-6814", "ENSG00000000001", three, 15.0, 18.0, 20.0)
    assert (rich["q25"], rich["q75"], rich["n_samples"]) == (15.0, 20.0, 3)

    thin = profile_row("E-MTAB-6814", "ENSG00000000001", two, 15.0, 18.0, 20.0)
    assert (thin["q25"], thin["q75"]) == (None, None)
    assert thin["median_abundance"] == 18.0
    assert thin["n_samples"] == 2
    assert thin["unit"] == UNIT == "tpm"


@pytest.mark.parametrize(
    ("percentile", "expected"),
    [(0, 0.0), (25, 24.75), (50, 49.5), (75, 74.25), (100, 99.0)],
)
def test_the_quantile_estimator_is_r_type_7(percentile: int, expected: float) -> None:
    """Pinned against worked literals, not against whatever a library returns."""
    assert quantile_linear(list(range(100)), percentile) == pytest.approx(expected)


def test_a_grid_is_101_rows_non_decreasing_and_states_its_denominator() -> None:
    """Ties at the floor are the normal case; ``n_genes`` names the population ranked over."""
    group = AssayGroup("g1", "7 week post conception", "heart", 3)
    medians = [0.0] * 46 + [float(i) for i in range(54)]
    rows = list(quantile_rows("E-MTAB-6814", group, medians))

    assert len(rows) == 101
    assert [row["percentile"] for row in rows] == list(range(101))
    values = [row["value"] for row in rows]
    assert all(a <= b for a, b in zip(values, values[1:], strict=False))
    assert values[0] == 0.0
    assert {row["n_genes"] for row in rows} == {100}
    assert {row["unit"] for row in rows} == {"tpm"}


def _config(groups: dict[str, tuple[str, int]]) -> str:
    body = ""
    for group_id, (label, n_assays) in groups.items():
        assays = "".join(f"<assay>ERR{group_id}{n}</assay>" for n in range(n_assays))
        body += f'<assay_group id="{group_id}" label="{label}">{assays}</assay_group>\n'
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n<configuration><analytics><assay_groups>'
        f"{body}</assay_groups></analytics></configuration>\n"
    )


GROUPS = {
    "g1": ("10 week post conception; heart", 3),
    "g2": ("9 week post conception; heart", 2),
    "g3": ("10 week post conception; brain", 4),
}


def _matrix(rows: list[tuple[str, list[str]]]) -> str:
    lines = ["GeneID\tGene Name\tg1\tg2\tg3"]
    lines += [f"{gene}\tSYM\t" + "\t".join(cells) for gene, cells in rows]
    return "\n".join(lines) + "\n"


# Gene order is deliberately not sorted, and ENSG...9 sorts after ENSG...10 only
# as a number: the package sorts as strings, like the atlas's mirrors.
ROWS = [
    ("ENSG00000000009", ["1,2,3,4,5", "NA", "10"]),
    ("ENSG00000000002", ["0,0,0,0,0", "6,7,8,9,10", "20,21,22,23,24"]),
    ("ENSG00000000005", ["", "", ""]),
]


@pytest.fixture
def pinned(monkeypatch: pytest.MonkeyPatch):  # type: ignore[no-untyped-def]
    def make(tmp_path: Path, config: str, matrix: str) -> tuple[Path, Path]:
        tpms = tmp_path / "m.tsv"
        xml = tmp_path / "c.xml"
        tpms.write_text(matrix)
        xml.write_text(config)
        monkeypatch.setattr(recipe, "MATRIX_SHA256", hashlib.sha256(tpms.read_bytes()).hexdigest())
        monkeypatch.setattr(
            recipe, "CONFIGURATION_SHA256", hashlib.sha256(xml.read_bytes()).hexdigest()
        )
        return tpms, xml

    return make


def _build(matrix: Path, config: Path, tmp_path: Path, **over: Any) -> Path:
    kwargs: dict[str, Any] = {
        "revision": 1,
        "built_by": BUILT_BY,
        "out": tmp_path / "pkg",
        "commit_of": lambda: RECIPE_COMMIT,
    }
    return recipe.build(matrix, config, **(kwargs | over))


def _read(package: Path, name: str) -> list[dict[str, str]]:
    blob = gzip.decompress((package / f"{name}.E-MTAB-6814.tsv.gz").read_bytes()).decode()
    lines = blob.splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _manifest(package: Path) -> dict[str, Any]:
    loaded: dict[str, Any] = json.loads((package / "manifest.json").read_text())
    return loaded


def test_every_matrix_gene_is_carried_under_its_ensembl_id(tmp_path: Path, pinned: Any) -> None:
    """No registry join: three matrix genes in, the two with a measurement out, by Ensembl id."""
    matrix, config = pinned(tmp_path, _config(GROUPS), _matrix(ROWS))
    out = _build(matrix, config, tmp_path)

    profiles = _read(out, "profiles")
    # Sorted (gene, tissue, stage) as strings; literals, not a re-sort of the output.
    assert [(r["gene"], r["tissue"], r["stage"]) for r in profiles] == [
        ("ENSG00000000002", "brain", "10 week post conception"),
        ("ENSG00000000002", "heart", "10 week post conception"),
        ("ENSG00000000002", "heart", "9 week post conception"),
        ("ENSG00000000009", "brain", "10 week post conception"),
        ("ENSG00000000009", "heart", "10 week post conception"),
    ]
    by_key = {(r["gene"], r["tissue"], r["stage"]): r for r in profiles}
    heart10 = by_key[("ENSG00000000009", "heart", "10 week post conception")]
    assert heart10 == {
        "dataset": "E-MTAB-6814",
        "gene": "ENSG00000000009",
        "tissue": "heart",
        "stage": "10 week post conception",
        "median_abundance": "3.0",
        "unit": "tpm",
        "q25": "2.0",
        "q75": "4.0",
        "n_samples": "3",
    }
    # A one-value cell in a four-assay group: quartiles present (n_samples >= 3).
    brain = by_key[("ENSG00000000009", "brain", "10 week post conception")]
    assert (brain["median_abundance"], brain["q25"], brain["n_samples"]) == ("10.0", "10.0", "4")
    # A two-assay group: median published, quartiles empty cells.
    thin = by_key[("ENSG00000000002", "heart", "9 week post conception")]
    assert (thin["median_abundance"], thin["q25"], thin["q75"], thin["n_samples"]) == (
        "8.0",
        "",
        "",
        "2",
    )

    manifest = _manifest(out)
    assert (manifest["ids"], manifest["coverage"]) == ("ensembl_gene_id", "every_gene")
    assert manifest["package"] == "expression_atlas_e_mtab_6814"
    assert manifest["version"] == package_version(
        upstream_release=None, input_sha256=recipe.MATRIX_SHA256
    )
    assert manifest["input_counts"] == {"rows": 3, "genes": 3}
    assert manifest["built_by"]["host"] == "workstation"
    assert [(t["name"], t["shard"], t["rows"], t["genes"]) for t in manifest["tables"]] == [
        ("profile_quantiles", "E-MTAB-6814", 303, None),
        ("profiles", "E-MTAB-6814", 5, 2),
    ]
    assert [i["name"] for i in manifest["inputs"]] == [
        "E-MTAB-6814-tpms.tsv",
        "E-MTAB-6814-configuration.xml",
    ]
    assert manifest["source"]["retrieved_on"] == "2026-08-19"


def test_an_unlabelled_group_refuses(tmp_path: Path, pinned: Any) -> None:
    """A label with no semicolon is not stage-and-organ; skipping it would drop a column."""
    groups = GROUPS | {"g4": ("a group with no split", 2)}
    matrix, config = pinned(tmp_path, _config(groups), _matrix(ROWS))
    with pytest.raises(PackageRefused, match="g4"):
        _build(matrix, config, tmp_path)
    assert not (tmp_path / "pkg").exists()


def test_quantiles_are_computed_over_every_gene_not_a_registry(tmp_path: Path, pinned: Any) -> None:
    """The heart/10-week grid ranks the three heart-measured genes, wherever they are."""
    rows = [(f"ENSG{n:011d}", [f"{n},{n},{n},{n},{n}", "", ""]) for n in range(1, 11)]
    rows.append(("ENSG00000000099", ["", "5", "5"]))
    matrix, config = pinned(tmp_path, _config(GROUPS), _matrix(rows))
    out = _build(matrix, config, tmp_path)

    grid = _read(out, "profile_quantiles")
    heart10 = [r for r in grid if (r["tissue"], r["stage"]) == ("heart", "10 week post conception")]
    assert [int(r["percentile"]) for r in heart10] == list(range(101))
    assert {r["n_genes"] for r in heart10} == {"10"}
    assert (heart10[0]["value"], heart10[50]["value"], heart10[100]["value"]) == (
        "1.0",
        "5.5",
        "10.0",
    )
    # The 9-week heart group is measured for one gene only.
    heart9 = [r for r in grid if r["stage"] == "9 week post conception"]
    assert {r["n_genes"] for r in heart9} == {"1"}
    # Order is (tissue, stage, percentile) with the percentile numeric: 2 before 10.
    assert [int(r["percentile"]) for r in heart9][:3] == [0, 1, 2]
    assert [int(r["percentile"]) for r in heart9][-3:] == [98, 99, 100]


def test_a_duplicate_matrix_gene_is_refused(tmp_path: Path, pinned: Any) -> None:
    matrix, config = pinned(tmp_path, _config(GROUPS), _matrix(ROWS + [ROWS[0]]))
    with pytest.raises(PackageRefused, match="ENSG00000000009"):
        _build(matrix, config, tmp_path)


def test_other_bytes_are_refused(tmp_path: Path) -> None:
    tpms = tmp_path / "m.tsv"
    xml = tmp_path / "c.xml"
    tpms.write_text(_matrix(ROWS))
    xml.write_text(_config(GROUPS))
    with pytest.raises(PackageRefused, match="sha256"):
        _build(tpms, xml, tmp_path)


def test_the_pins_are_the_measured_files() -> None:
    assert (
        recipe.MATRIX_SHA256 == "92484a979fb8e262fe1f535b2b902ba8088310bbbec76969ba873b51d593f4be"
    )
    assert recipe.CONFIGURATION_SHA256 == (
        "31de8dd2c78067cba2a333527d1af123e9d497443c4acdf7c5dea2da115b166f"
    )
    assert package_version(upstream_release=None, input_sha256=recipe.MATRIX_SHA256) == "92484a97"


def test_a_breakpoint_is_rounded_to_six_places() -> None:
    """The converter wrote ``round(value, 6)``; the atlas's committed mirror holds that spelling."""
    group = AssayGroup("g1", "7 week post conception", "heart", 3)
    rows = list(quantile_rows("E-MTAB-6814", group, [0.0, 1 / 3]))
    assert rows[50]["value"] == 0.166667
