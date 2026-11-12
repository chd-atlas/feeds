# orphanet_product4

Orphadata product 4 (`en_product4.xml`, CC-BY-4.0) as the `orphanet_product4` package:
table `orphanet_phenotypes` v1, `ids: hgnc_id` (no gene column), `coverage: every_gene`.

The package carries **every** HPO annotation. The atlas's converter kept only HPO ids
under HP:0001626; the atlas now filters at ingest (owner decision F23). Exact-duplicate
annotations collapse into one `dropped` entry; rows differing only in frequency band are
both kept. The build accepts only the download with sha256 `4f44e8a6...` (retrieved
2026-10-25); the package version is `4f44e8a6`, and `source.file_version` / `file_date`
are the file's own `<JDBOR>` header.

## Measurements

Real build, 2026-11-12, from the pinned file: 116,664 annotations read, 116,640 rows,
`dropped` 24 (first example `ORPHA:101011 HP:0001288`). Every row of the atlas's committed
`mirrors/orphanet_phenotypes.tsv` (7,902 rows at atlas 4b1f8f6) appears byte-for-byte in the
package table, in the same relative order.

## Tests

`pytest sources/orphanet_product4`: synthetic XML, no network. One test reads
`local/orphanet/en_product4.xml` and skips, naming it, when absent.
