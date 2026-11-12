# mondo_sssom

The Mondo SSSOM mapping file (`mondo.sssom.tsv`, release `v2026-10-06`, CC-BY-4.0) as the
`mondo_sssom` package: table `mondo_xrefs` v1, `ids: hgnc_id` (no gene column),
`coverage: every_gene`, package version `2026-10-06`.

A row is a `skos:exactMatch` from a `MONDO:` subject to an `Orphanet:` or `OMIM:` object
(`Orphanet:<n>` is written `ORPHA:<n>`), de-duplicated and sorted by `(mondo, target)`.
Other predicates and targets are neither rows nor drops. The package carries **every**
MONDO term; `source.name` keeps the atlas's wording ("for terms the validity mirrors
name"), which describes what the atlas mirrors after its ingest filter (F23), not what
this package holds. The build accepts only the download with sha256 `ec15afdb...`
(retrieved 2026-10-25).

## Build

Download (the release tag is in the URL, so it is stable; `local/` is git-ignored):

```bash
mkdir -p local/mondo
curl -L -o local/mondo/mondo.sssom.tsv https://raw.githubusercontent.com/monarch-initiative/mondo/v2026-10-06/src/ontology/mappings/mondo.sssom.tsv
```

```python
from pathlib import Path
from sources.mondo_sssom import recipe

recipe.build(
    Path("local/mondo/mondo.sssom.tsv"),
    revision=1,
    built_by={"host": "<machine>"},  # recipe_commit is read from the clean checkout
    out=Path("dist/mondo_sssom"),
)
```
The build refuses any other file by sha256.

## Measurements

Real build, 2026-11-12, from the pinned file: 19,832 carried matches read, 19,832 rows,
nothing dropped. Every row of the atlas's committed `mirrors/mondo_xrefs.tsv` (10,011
rows, atlas 4b1f8f6) appears byte for byte in the table, in the same relative order; the
other 9,821 rows (over 8,810 distinct MONDO terms, none of them a term the mirror carries)
are matches for terms no validity mirror names, which the atlas filters at ingest.

## Tests

`pytest sources/mondo_sssom`: synthetic files, no network. One test reads
`local/mondo/mondo.sssom.tsv` and skips, naming it, when absent.
