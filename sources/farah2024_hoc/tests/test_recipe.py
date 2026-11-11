"""The Farah 2024 recipe on a synthetic pseudobulk: no network, no committed mirror.

Ported from the atlas's ``tests/unit/test_convert_farah2024.py`` (its direct
synthetic tests; the ones that read the committed mirror are left behind). The
join is no longer by HGNC id, so the colliding-alias cases became
``test_the_slice_keeps_registry_approved_symbols_only``.

The risks pinned here:

- the denominator of a CPM is the whole library, not the registry subset's;
- a donor is a heart (age and replicate), not an age;
- a quartile of fewer than three hearts is withheld;
- the slice keeps approved symbols only, so an alias that is another gene's
  approved symbol (``HGF`` is SOS1's alias) never reaches the package.

Each guard has a test that fails without it; the mutations were applied by hand
from an in-memory copy (results in the commit message).
"""

import gzip
import hashlib
import importlib.util
import json
from datetime import date
from pathlib import Path
from typing import Any

import pytest

from feeds.package import PackageRefused
from feeds.version import package_version
from sources.farah2024_hoc import recipe

BUILT_BY = {"host": "login01", "hvantk_release": "96ceafb4"}
RECIPE_COMMIT = "a" * 40
JOB_CODE_URL = "https://github.com/chd-atlas/chd-atlas.github.io/tree/" + "6" * 40 + "/scripts/hpc"
ATLAS_COMMIT = "f" * 40
MADE_ON = date(2026, 11, 2)
H5AD_SHA = "19a1d4961c2f38194b7589e86dc53cae45c272fc0f9d8a7e907edeefe3215c2c"
HEADER = "symbol\theart\tcell_type\tregion\tcounts\tcount_nonzero\tn_cells\tlibrary_size"
REGISTRY_HEADER = "hgnc_id\tsymbol\tname\taliases\tprev_symbols"
# Three rows' symbol columns transcribed from the atlas's mirrors/genes.tsv
# (2026-10-28): each carries an alias that is an unrelated gene's approved symbol.
REGISTRY_ROWS = [
    "HGNC:1\tAAA\tn\t\t",
    "HGNC:11187\tSOS1\tn\tHGF|GF1\tGINGF",
    "HGNC:1152\tPOPDC1\tn\tPOP1|HBVES\tBVES",
]


def _line(
    heart: str = "9pcw-1",
    *,
    counts: float = 50.0,
    library: float = 1_000_000.0,
    nonzero: int = 5,
    cells: int = 10,
    cell_type: str = "cardiomyocyte",
    region: str = "LV",
    symbol: str = "AAA",
) -> str:
    return f"{symbol}\t{heart}\t{cell_type}\t{region}\t{counts}\t{nonzero}\t{cells}\t{library}"


def _inputs(
    tmp_path: Path, lines: list[str], registry: list[str] | None = None
) -> tuple[Path, Path]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    pseudobulk = tmp_path / "pseudobulk.tsv"
    pseudobulk.write_text("\n".join([HEADER, *lines]) + "\n", encoding="utf-8")
    genes = tmp_path / "genes.tsv"
    genes.write_text(
        "\n".join([REGISTRY_HEADER, *(REGISTRY_ROWS if registry is None else registry)]) + "\n",
        encoding="utf-8",
    )
    return pseudobulk, genes


def _kwargs(tmp_path: Path, genes: Path) -> dict[str, Any]:
    return {
        "atlas_commit": ATLAS_COMMIT,
        "h5ad_sha256": H5AD_SHA,
        "registry_sha256": hashlib.sha256(genes.read_bytes()).hexdigest(),
        "pseudobulk_made_on": MADE_ON,
        "job_code_url": JOB_CODE_URL,
        "revision": 1,
        "built_by": BUILT_BY,
        "out": tmp_path / "pkg",
        "commit_of": lambda: RECIPE_COMMIT,
    }


def _build(tmp_path: Path, lines: list[str], **over: Any) -> Path:
    pseudobulk, genes = _inputs(tmp_path, lines, over.pop("registry", None))
    return recipe.build(pseudobulk, genes, **(_kwargs(tmp_path, genes) | over))


def _rows(pkg: Path) -> list[dict[str, str]]:
    text = gzip.decompress((pkg / "cell_type_profiles.hoc.tsv.gz").read_bytes()).decode()
    header, *body = [line.split("\t") for line in text.splitlines()]
    return [dict(zip(header, record, strict=True)) for record in body]


def _manifest(pkg: Path) -> dict[str, Any]:
    loaded: dict[str, Any] = json.loads((pkg / "manifest.json").read_text(encoding="utf-8"))
    return loaded


def test_cpm_divides_by_the_whole_library_and_the_columns_are_the_tables(tmp_path: Path) -> None:
    (row,) = _rows(_build(tmp_path, [_line(counts=50.0, library=1_000_000.0)]))
    assert row == {
        "dataset": "hoc",
        "gene": "AAA",
        "cell_type": "cardiomyocyte",
        "region": "LV",
        "median_abundance": "50.0",
        "unit": "cpm",
        "q25": "",
        "q75": "",
        "n_donors": "1",
        "n_cells": "10",
        "fraction_expressing": "0.5",
        "source": "farah2024_hoc",
    }
    halved = _build(tmp_path / "b", [_line(counts=50.0, library=2_000_000.0)])
    assert _rows(halved)[0]["median_abundance"] == "25.0"


def test_a_heart_is_age_and_replicate(tmp_path: Path) -> None:
    (row,) = _rows(_build(tmp_path, [_line("9pcw-1"), _line("9pcw-2")]))
    assert row["n_donors"] == "2"
    assert row["n_cells"] == "20"


def test_quartiles_need_three_hearts(tmp_path: Path) -> None:
    (two,) = _rows(_build(tmp_path / "a", [_line("9pcw-1", counts=10), _line("9pcw-2", counts=20)]))
    assert two["q25"] == "" and two["q75"] == ""
    assert two["median_abundance"] == "15.0"
    hearts = [_line("9pcw-1", counts=10), _line("10pcw-1", counts=20), _line("11pcw-1", counts=40)]
    (three,) = _rows(_build(tmp_path / "b", hearts))
    assert (three["q25"], three["median_abundance"], three["q75"]) == ("15.0", "20.0", "30.0")
    assert three["n_donors"] == "3"


def test_share_of_cells_pools_across_hearts(tmp_path: Path) -> None:
    lines = [_line("9pcw-1", nonzero=3, cells=10), _line("9pcw-2", nonzero=1, cells=30)]
    (row,) = _rows(_build(tmp_path, lines))
    assert float(row["fraction_expressing"]) == pytest.approx(0.1)


def test_groups_are_split_by_cell_type_and_region_and_sorted(tmp_path: Path) -> None:
    lines = [
        _line(region="RV"),
        _line(region="LV"),
        _line(cell_type="aaa", region="LV"),
        _line(symbol="SOS1", region="LV"),
    ]
    rows = _rows(_build(tmp_path, lines))
    assert [(r["gene"], r["cell_type"], r["region"]) for r in rows] == [
        ("AAA", "aaa", "LV"),
        ("AAA", "cardiomyocyte", "LV"),
        ("AAA", "cardiomyocyte", "RV"),
        ("SOS1", "cardiomyocyte", "LV"),
    ]


def test_the_slice_keeps_registry_approved_symbols_only(tmp_path: Path) -> None:
    """HGF is SOS1's alias and hepatocyte growth factor's approved symbol; POP1 is
    POPDC1's alias and the ribonuclease P subunit's; BVES is POPDC1's previous
    symbol. None may reach the package, and the build must not stop on them."""
    lines = [
        _line(symbol="SOS1"),
        _line(symbol="HGF", counts=999),
        _line(symbol="POP1", counts=999),
        _line(symbol="BVES", counts=999),
        _line(symbol="ZZZ", counts=999),
    ]
    pkg = _build(tmp_path, lines)
    assert [r["gene"] for r in _rows(pkg)] == ["SOS1"]
    (entry,) = _manifest(pkg)["dropped"]
    assert entry == {
        "reason": "symbol is not an approved symbol of a registry gene",
        "count": 4,
        "examples": ["BVES", "HGF", "POP1", "ZZZ"],
    }
    assert _manifest(pkg)["input_counts"] == {"rows": 5, "genes": 5}


def test_the_manifest_records_the_slice_commit_and_registry_size(tmp_path: Path) -> None:
    pkg = _build(tmp_path, [_line()])
    manifest = _manifest(pkg)
    assert manifest["coverage"] == {
        "registry_slice": {"atlas_commit": ATLAS_COMMIT, "genes": len(REGISTRY_ROWS)}
    }
    assert manifest["ids"] == "symbol"
    assert manifest["package"] == "farah2024_hoc"
    (table,) = manifest["tables"]
    assert (table["name"], table["shard"], table["file"]) == (
        "cell_type_profiles",
        "hoc",
        "cell_type_profiles.hoc.tsv.gz",
    )
    assert manifest["dropped"] == []


def test_the_version_is_the_h5ad_prefix_and_both_inputs_are_recorded(tmp_path: Path) -> None:
    pkg = _build(tmp_path, [_line()])
    pseudobulk = tmp_path / "pseudobulk.tsv"
    manifest = _manifest(pkg)
    assert (
        manifest["version"]
        == "19a1d496"
        == package_version(upstream_release=None, input_sha256=H5AD_SHA)
    )
    h5ad, pb = manifest["inputs"]
    assert (h5ad["name"], h5ad["sha256"], h5ad["retrieved_on"]) == (
        "farah_2024.h5ad",
        H5AD_SHA,
        "2026-04-18",
    )
    assert h5ad["url"] == "https://cells.ucsc.edu/?ds=hoc+all-heart"
    assert pb["retrieved_on"] == "2026-11-02"
    assert pb["url"] == JOB_CODE_URL
    assert pb["sha256"] == hashlib.sha256(pseudobulk.read_bytes()).hexdigest()
    assert manifest["source"]["version"] == f"sha256:{H5AD_SHA}"
    assert manifest["source"]["redistribution"] == "permitted_with_attribution"


def test_the_job_code_url_date_and_build_commit_reach_the_manifest(tmp_path: Path) -> None:
    other_url = "https://github.com/chd-atlas/feeds/tree/" + "c" * 40 + "/sources/x/hpc"
    first = _manifest(_build(tmp_path / "a", [_line()]))
    second = _manifest(
        _build(
            tmp_path / "b",
            [_line()],
            pseudobulk_made_on=date(2026, 12, 1),
            job_code_url=other_url,
            commit_of=lambda: "b" * 40,
        )
    )
    one, two = first["inputs"][1], second["inputs"][1]
    assert (one["retrieved_on"], two["retrieved_on"]) == ("2026-11-02", "2026-12-01")
    assert (one["url"], two["url"]) == (JOB_CODE_URL, other_url)
    # The input's URL names the job's code; built_by names the code that ran this build.
    assert first["built_by"] == {**BUILT_BY, "recipe_commit": RECIPE_COMMIT}
    assert second["built_by"]["recipe_commit"] == "b" * 40


@pytest.mark.parametrize("commit", ["dev", "A" * 40, "a" * 39, ""])
def test_a_malformed_build_commit_is_refused(tmp_path: Path, commit: str) -> None:
    with pytest.raises(PackageRefused, match="40-hex"):
        _build(tmp_path, [_line()], commit_of=lambda: commit)


def test_a_supplied_recipe_commit_must_equal_the_builds_own(tmp_path: Path) -> None:
    same = {**BUILT_BY, "recipe_commit": RECIPE_COMMIT}
    assert _manifest(_build(tmp_path / "a", [_line()], built_by=same))["built_by"] == same
    with pytest.raises(PackageRefused, match="recipe_commit"):
        _build(tmp_path / "b", [_line()], built_by={**BUILT_BY, "recipe_commit": "b" * 40})


@pytest.mark.parametrize(
    ("over", "match"),
    [
        ({"built_by": {"host": "login01"}}, "hvantk_release"),
        ({"built_by": {"host": "login01", "hvantk_release": " "}}, "hvantk_release"),
        ({"job_code_url": "https://github.com/chd-atlas/feeds/tree/main/sources"}, "job_code_url"),
        (
            {"job_code_url": "https://github.com/other/feeds/tree/" + "a" * 40 + "/x"},
            "job_code_url",
        ),
        (
            {"job_code_url": "https://github.com/chd-atlas/feeds/blob/" + "a" * 40 + "/x"},
            "job_code_url",
        ),
        ({"registry_sha256": "b" * 64}, "sliced with"),
        ({"registry_sha256": "B" * 64}, "registry_sha256"),
    ],
)
def test_provenance_the_job_printed_is_required_and_checked(
    tmp_path: Path, over: dict[str, Any], match: str
) -> None:
    with pytest.raises(PackageRefused, match=match):
        _build(tmp_path, [_line()], **over)
    assert not (tmp_path / "pkg").exists()


def test_the_revision_reaches_the_version(tmp_path: Path) -> None:
    assert _manifest(_build(tmp_path, [_line()], revision=2))["version"] == "19a1d496-r2"


def test_two_builds_are_byte_identical(tmp_path: Path) -> None:
    lines = [_line("9pcw-1"), _line("9pcw-2", counts=70), _line(symbol="SOS1")]
    a = _build(tmp_path / "a", lines)
    b = _build(tmp_path / "b", lines)
    names = sorted(p.name for p in a.iterdir())
    assert names == sorted(p.name for p in b.iterdir())
    assert all((a / n).read_bytes() == (b / n).read_bytes() for n in names)


def test_a_heart_counted_twice_in_one_group_is_refused(tmp_path: Path) -> None:
    with pytest.raises(PackageRefused, match="heart 9pcw-1 appears twice"):
        _build(tmp_path, [_line("9pcw-1"), _line("9pcw-1")])


@pytest.mark.parametrize("symbol", ["AAA", "ZZZ"])
def test_a_group_under_ten_cells_stops_the_build(tmp_path: Path, symbol: str) -> None:
    """hvantk drops a group under `--min-cells`, and the job passes 10. A group
    under 10 cells means the job ran without it, so its whole output is
    suspect: refused, not dropped, whether or not the row's symbol is kept."""
    small = _line("10pcw-1", cells=9, nonzero=1, symbol=symbol)
    with pytest.raises(PackageRefused, match="10pcw-1/cardiomyocyte/LV has n_cells 9"):
        _build(tmp_path, [_line("9pcw-1"), small])


@pytest.mark.parametrize("padded", ["AAA ", " AAA", "ZZZ "])
@pytest.mark.parametrize("beside_bare", [True, False])
def test_a_padded_symbol_is_refused_by_name_not_stripped(
    tmp_path: Path, padded: str, beside_bare: bool
) -> None:
    lines = [_line(symbol="AAA"), _line(symbol=padded)] if beside_bare else [_line(symbol=padded)]
    with pytest.raises(PackageRefused, match=f"symbol {padded!r} has leading or trailing"):
        _build(tmp_path, lines)


@pytest.mark.parametrize(
    "line",
    [
        "AAA\t9pcw-1\tcardiomyocyte\tLV\tNaN\t5\t10\t1000000",
        "AAA\t9pcw-1\tcardiomyocyte\tLV\t50\t5\t10\t0",
        "AAA\t9pcw-1\tcardiomyocyte\tLV\t50\t11\t10\t1000000",
        "AAA\t9pcw-1\tcardiomyocyte\tLV\t50\t5\t10",
        "AAA\t9pcw-1\tcardiomyocyte\tleft atrium\t50\t5\t10\t1000000",
        # a row outside the slice is checked too: it would otherwise vanish unread
        "ZZZ\t9pcw-1\tcardiomyocyte\tleft atrium\t50\t5\t10\t1000000",
        "AAA\t9pcw-1\tcardiomyocyte\tLV\t-1\t5\t10\t1000000",
        "AAA\t9pcw-1\tcardiomyocyte\tLV\t50\t5\t10.5\t1000000",
        "\t9pcw-1\tcardiomyocyte\tLV\t50\t5\t10\t1000000",
    ],
)
def test_nothing_is_written_when_a_row_is_malformed(tmp_path: Path, line: str) -> None:
    good = "AAA\t9pcw-2\tcardiomyocyte\tLV\t50\t5\t10\t1000000"
    with pytest.raises(PackageRefused):
        _build(tmp_path, [good, line])
    assert not (tmp_path / "pkg").exists()


def test_a_wrong_header_is_refused(tmp_path: Path) -> None:
    pseudobulk, genes = _inputs(tmp_path, [_line()])
    pseudobulk.write_text("a\tb\n1\t2\n", encoding="utf-8")
    with pytest.raises(PackageRefused, match="header is"):
        recipe.build(pseudobulk, genes, **_kwargs(tmp_path, genes))


@pytest.mark.parametrize(
    ("over", "match"),
    [
        ({"h5ad_sha256": "19A1D496"}, "h5ad_sha256"),
        ({"h5ad_sha256": "a" * 63}, "h5ad_sha256"),
        ({"atlas_commit": "f094c28"}, "atlas_commit"),
        ({"atlas_commit": "F" * 40}, "atlas_commit"),
    ],
)
def test_a_malformed_digest_or_commit_is_refused(
    tmp_path: Path, over: dict[str, str], match: str
) -> None:
    with pytest.raises(PackageRefused, match=match):
        _build(tmp_path, [_line()], **over)


def test_a_registry_listing_a_gene_twice_or_an_empty_slice_is_refused(tmp_path: Path) -> None:
    with pytest.raises(PackageRefused, match="listed twice"):
        _build(tmp_path / "a", [_line()], registry=[*REGISTRY_ROWS, "HGNC:1\tAAB\tn\t\t"])
    with pytest.raises(PackageRefused, match="listed twice"):
        _build(tmp_path / "b", [_line()], registry=[*REGISTRY_ROWS, "HGNC:2\tAAA\tn\t\t"])
    with pytest.raises(PackageRefused, match="no pseudobulk row"):
        _build(tmp_path / "c", [_line(symbol="HGF")])


# --- the cluster job: loaded by path, since hpc/ is not a package -------------------------


def _hpc(name: str) -> Path:
    return Path(recipe.__file__).resolve().parent / "hpc" / name


def _extract() -> Any:
    spec = importlib.util.spec_from_file_location("farah_extract", _hpc("farah_extract.py"))
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_extract_rows_keeps_only_wanted_symbols_and_tags_the_region() -> None:
    extract = _extract()
    groups = [
        {"Age": 9, "Rep": 1, "major_cell_class": "cardiomyocyte", "Region": "LV", "n_cells": 12},
        {"Age": 10, "Rep": 2, "major_cell_class": "fibroblast", "Region": "RA", "n_cells": 30},
    ]
    genes = ["AAA", "BBB", "CCC"]
    sums = [[1.0, 2.0, 7.0], [4.0, 4.0, 2.0]]
    nonzero = [[1, 2, 5], [3, 3, 2]]
    rows = extract.extract_rows(groups, genes, sums, nonzero, {"AAA", "CCC"}, None)
    assert [(r["symbol"], r["heart"], r["counts"], r["library_size"]) for r in rows] == [
        ("AAA", "9pcw-1", 1.0, 10.0),
        ("CCC", "9pcw-1", 7.0, 10.0),
        ("AAA", "10pcw-2", 4.0, 10.0),
        ("CCC", "10pcw-2", 2.0, 10.0),
    ]
    assert [r["region"] for r in rows] == ["LV", "LV", "RA", "RA"]
    pooled = extract.extract_rows(groups, genes, sums, nonzero, {"BBB"}, "all")
    assert {r["region"] for r in pooled} == {"all"}


def test_the_extract_keeps_only_the_registry_s_approved_symbols(tmp_path: Path) -> None:
    registry = tmp_path / "genes.tsv"
    registry.write_text(
        "hgnc_id\tsymbol\taliases\tprev_symbols\n"
        "HGNC:11187\tSOS1\tHGF|GF1\tGINGF\n"
        "HGNC:1152\tPOPDC1\tPOP1|HBVES\tBVES\n",
        encoding="utf-8",
    )
    assert _extract().read_wanted(registry) == {"SOS1", "POPDC1"}


def test_the_cluster_files_name_no_partition_and_no_mail_user() -> None:
    sbatch = _hpc("farah_celltypes.sbatch").read_text(encoding="utf-8")
    directives = [x for x in sbatch.splitlines() if x.startswith("#SBATCH")]
    assert not [d for d in directives if "--partition" in d or "--mail-user" in d]
    assert "#SBATCH --time=0-00:30:00" in directives
    assert "hpc-env/13.1" in sbatch
    submit = _hpc("submit.sh").read_text(encoding="utf-8")
    assert "#SBATCH" not in submit  # submit.sh passes --partition on the command line


def test_the_job_requires_and_prints_its_code_permalink_and_the_registry_digest() -> None:
    sbatch = _hpc("farah_celltypes.sbatch").read_text(encoding="utf-8")
    submit = _hpc("submit.sh").read_text(encoding="utf-8")
    assert 'JOB_CODE_URL="${JOB_CODE_URL:?' in sbatch
    assert 'JOB_CODE_URL="${JOB_CODE_URL:?' in submit
    assert 'echo "job code url: $JOB_CODE_URL"' in sbatch
    assert 'echo "registry sha256: $(sha256sum "$REGISTRY"' in sbatch
    # the shell check and recipe.py's check are the same pattern
    assert recipe._JOB_CODE_URL.pattern.removesuffix(".*") in sbatch
