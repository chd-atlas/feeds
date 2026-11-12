# gencc

GenCC submissions export (CC0-1.0) as the `gencc` package: table `gencc_submissions`
v1, `ids: hgnc_id`, `coverage: every_gene`. GenCC states no release, so the version
is the input's sha256[:8] and `source.version` is `retrieved <retrieved_on>`.

Every submission is kept; the key is `(gene, disease, moi, submitter, sgc_id)` (see
`recipe.py` for why `sgc_id` is in it). The export is quoted TSV, parsed with polars,
never by splitting lines.

Two entry points in `recipe.py`:

- `build(tsv, retrieved_on=..., revision=..., built_by=..., out=...)` parses a live
  download from `https://search.thegencc.org/download/action/submissions-export-tsv?format=new`.
- `bootstrap(mirror, atlas_commit=..., ...)` packages the atlas's committed
  `mirrors/gencc_submissions.tsv` once (decision F21), because the 2026-08-03 export
  cannot be fetched again. It accepts only the file with sha256 `88c77272...` (atlas
  4b1f8f6, 30,410 rows, 6,113 genes); the package version is `88c77272`.

## Measurements

- Live parity, 2026-11-12: `read_rows` over the 2026-11-11 export gives 30,376
  records (46,968 physical lines) and `dropped == []`; rendered as TSV it is
  `cmp`-equal to the atlas's `scripts/convert_gencc.py` (archive of atlas 4b1f8f6)
  on the same file.
- A bootstrap build of the real mirror into a scratch directory succeeded; its table
  is byte-identical to the mirror.

## Tests

`pytest sources/gencc`: synthetic export and mirror, no network.
