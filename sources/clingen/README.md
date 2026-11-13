# clingen

ClinGen Gene-Disease Validity bulk CSV (CC0-1.0) as the `clingen` package: table
`clingen_validity` v2, `ids: hgnc_id`, `coverage: every_gene`. Version is the
banner's `FILE CREATED` date; `source.version` is the banner text.

`recipe.build(csv, retrieved_on=..., revision=..., built_by=..., out=...)` parses a
live download from `https://search.clinicalgenome.org/kb/gene-validity/download`.

**v2 (F28): the key is (gene, disease, moi, gcep).** ClinGen's file lists FOXE3 /
MONDO:0019625 / AD twice, by two GCEPs (SOP4 Moderate 2016, SOP12 Limited 2026); both
rows are kept, sorted by panel. One GCEP twice on one triple is refused, never
deduplicated.

**No `bootstrap`.** The first release, `clingen/2026-08-03` (table v1), was
bootstrapped from the atlas's committed mirror at atlas 4b1f8f6 (decision F21). Under
v2 that entry point could no longer reproduce its own release, so it was deleted; the
release stays published, and git history holds the code that built it.

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

## Measurements

- Live parity, 2026-11-12: `read_rows` over the 2026-11-11 download (3,730 rows)
  rendered as TSV is `cmp`-equal to the atlas's `scripts/convert_clingen.py`
  (archive of atlas 4b1f8f6) on the same file, under v1's sort. v2 also sorts by panel,
  which reorders only FOXE3's two rows (v1's stable sort left them in file order, HTAAD
  first).
- That download repeats the v1 (gene, disease, moi) key for FOXE3, which v1's primary
  key refused, so `build` on it stopped in `write_package`; parity is row-level.

## Tests

`pytest sources/clingen`: synthetic CSV, no network.
