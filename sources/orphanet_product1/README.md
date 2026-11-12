# orphanet_product1

Orphadata product 1 (`en_product1.xml`, CC-BY-4.0) as the `orphanet_product1` package:
table `orphanet_omim` v1, `ids: hgnc_id` (no gene column), `coverage: every_gene`.

Every `ExternalReference` whose Source is `OMIM` is a row, whatever its relation or
validation status; references of other sources are neither rows nor drops. The build
accepts only the download with sha256 `df8d562a...` (retrieved 2026-10-25); the package
version is `df8d562a`, and `source.file_version` / `file_date` are the file's own
`<JDBOR>` header.

## Measurements

Real build, 2026-11-12, from the pinned file: 8,745 OMIM references read, 8,745 rows,
nothing dropped. The table is byte-identical to the atlas's committed
`mirrors/orphanet_omim.tsv` (atlas 4b1f8f6).

## Tests

`pytest sources/orphanet_product1`: synthetic XML, no network. One test reads
`local/orphanet/en_product1.xml` and skips, naming it, when absent.
