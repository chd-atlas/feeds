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

## Measurements

- Live parity, 2026-11-12: `read_rows` over the 2026-11-11 download (3,730 rows)
  rendered as TSV is `cmp`-equal to the atlas's `scripts/convert_clingen.py`
  (archive of atlas 4b1f8f6) on the same file.
- That download repeats the (gene, disease, moi) key for FOXE3, which the table's
  primary key refuses, so `build` on it stops in `write_package`; parity is row-level.
- A bootstrap build of the real mirror into a scratch directory succeeded.

## Tests

`pytest sources/clingen`: synthetic CSV and mirror, no network.
