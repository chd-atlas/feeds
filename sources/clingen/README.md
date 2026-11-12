# clingen

ClinGen Gene-Disease Validity bulk CSV (CC0-1.0) as the `clingen` package: table
`clingen_validity` v1, `ids: hgnc_id`, `coverage: every_gene`. Version is the
banner's `FILE CREATED` date; `source.version` is the banner text.

Two entry points in `recipe.py`:

- `build(csv, retrieved_on=..., revision=..., built_by=..., out=...)` parses a live
  download from `https://search.clinicalgenome.org/kb/gene-validity/download`.
- `bootstrap(mirror, atlas_commit=..., ...)` packages the atlas's committed
  `mirrors/clingen_gene_validity.tsv` once (decision F21), because the 2026-08-03
  download cannot be fetched again. It accepts only the file with sha256
  `b4c89cb6...` (atlas 4b1f8f6, 3,653 rows, 3,017 genes).

## Build

Live (ClinGen serves only the current file, so the sha256 is whatever you got):

```bash
mkdir -p local/clingen
curl -L -o local/clingen/clingen.csv https://search.clinicalgenome.org/kb/gene-validity/download
```

```python
from datetime import date
from pathlib import Path
from sources.clingen import recipe

recipe.build(
    Path("local/clingen/clingen.csv"),
    retrieved_on=date.today(),
    revision=1,
    built_by={"host": "<machine>"},
    out=Path("dist/clingen"),
)
```

Bootstrap, from the atlas mirror at the pinned commit (`recipe.ATLAS_COMMIT`,
`4b1f8f68cbc30beaf932fbbc58e3617ebe401b6d`; any other `atlas_commit` is refused, because
the manifest's input URL is a permanent permalink to it):

```bash
mkdir -p local/clingen
git -C <atlas> show 4b1f8f6:mirrors/clingen_gene_validity.tsv > local/clingen/clingen_gene_validity.tsv
```

```python
recipe.bootstrap(
    Path("local/clingen/clingen_gene_validity.tsv"),
    atlas_commit=recipe.ATLAS_COMMIT,
    revision=1,
    built_by={"host": "<machine>"},
    out=Path("dist/clingen"),
)
```

## Measurements

- Live parity, 2026-11-12: `read_rows` over the 2026-11-11 download (3,730 rows)
  rendered as TSV is `cmp`-equal to the atlas's `scripts/convert_clingen.py`
  (archive of atlas 4b1f8f6) on the same file.
- That download repeats the (gene, disease, moi) key for FOXE3, which the table's
  primary key refuses, so `build` on it stops in `write_package`; parity is row-level.
- A bootstrap build of the real mirror into a scratch directory succeeded.

## Tests

`pytest sources/clingen`: synthetic CSV and mirror, no network.
