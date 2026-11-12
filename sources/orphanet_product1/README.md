# orphanet_product1

Orphadata product 1 (`en_product1.xml`, CC-BY-4.0) as the `orphanet_product1` package:
table `orphanet_omim` v1, `ids: hgnc_id` (no gene column), `coverage: every_gene`.

Every `ExternalReference` whose Source is `OMIM` is a row, whatever its relation or
validation status; references of other sources are neither rows nor drops. The build
accepts only the download with sha256 `df8d562a...` (retrieved 2026-10-25); the package
version is `df8d562a`, and `source.file_version` / `file_date` are the file's own
`<JDBOR>` header.

## Build

Download (not done by the recipe; `local/` is git-ignored):

```bash
mkdir -p local/orphanet
curl -L -o local/orphanet/en_product1.xml https://www.orphadata.com/data/xml/en_product1.xml
```

```python
from pathlib import Path
from sources.orphanet_product1 import recipe

recipe.build(
    Path("local/orphanet/en_product1.xml"),
    revision=1,
    built_by={"host": "<machine>"},  # recipe_commit is read from the clean checkout
    out=Path("dist/orphanet_product1"),
)
```
The build refuses any other file by sha256, so a download made after the pin may be refused; 
see Conventions in the root README for `revision`.

## Measurements

Real build, 2026-11-12, from the pinned file: 8,745 OMIM references read, 8,745 rows,
nothing dropped. The table is byte-identical to the atlas's committed
`mirrors/orphanet_omim.tsv` (atlas 4b1f8f6).

## Tests

`pytest sources/orphanet_product1`: synthetic XML, no network. One test reads
`local/orphanet/en_product1.xml` and skips, naming it, when absent.
