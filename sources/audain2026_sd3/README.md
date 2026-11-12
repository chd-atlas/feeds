# audain2026_sd3

Audain et al. 2026, npj Genomic Medicine 11, Supplementary Data 3: per-gene rare
variant burden for the whole exome, three strata (`aCHD_vs_controls`,
`sCHD_vs_controls`, `nsCHD_vs_controls`) against UK Biobank controls. Published
as one table, `burden` (one row per sheet row), with `ids: symbol` and
`coverage: every_gene`. The atlas's `chd-atlas ingest` resolves the symbols
(approved, then alias, then previous; owner decision F22) and keeps its
registry's genes, which recovers CRIPTO from `TDGF1`, POPDC1 from `BVES` and
CIROZ from `C1orf127`.

`source.yaml` holds the name, the Springer URL, the licence (CC-BY-4.0,
`redistribution: permitted_with_attribution`) and its evidence, as the atlas's
`mirrors/sources.yaml` records them.

## Fetch and build

The workbook is downloaded by hand from the URL in `source.yaml` into
`local/audain2026/` (ignored by git). `build` refuses any file whose sha256 is
not `8ebb6fdceb38b0099d6229a1ebf4d77db9f9492a62cfafb569d2865f9645ee2b`, so the
version is `8ebb6fdc` at revision 1.

```python
from pathlib import Path
from sources.audain2026_sd3 import recipe

recipe.build(
    Path("local/audain2026/41525_2026_582_MOESM4_ESM.xlsx"),
    revision=1,
    built_by={"host": "<machine>"},  # recipe_commit is read from the clean checkout
    out=Path("dist/audain2026_sd3"),
)
```

## Rules

Moved from the atlas's `scripts/convert_burden.py`, whose docstring (kept in
`recipe.py`) has the measurements behind each rule.

- **`gene` is the workbook's `SYMBOL` verbatim**; no registry join here.
- A symbol shaped `\d{4}-\d{2}-\d{2}` (an Excel date) is published in `dropped`
  as "Excel date corruption" with a row count and up to five distinct examples,
  never guessed.
- Each sheet is its stratum and carries its own carrier and denominator columns;
  `Infinity` odds ratios publish no effect and `effect_bound: unbounded_above`;
  `syn` rows are kept as the negative control.
- `input_counts.rows` is the sheet rows read (dropped ones included) and
  `genes` the distinct symbols read.
- Rows are sorted on the table's sort key with `maf_max` compared as a number.

## Measurements

Measured 2026-11-12 by `recipe.build` on the pinned workbook (sha256 `8ebb6fdc...`),
revision 1, run in about 20 s:

- `input_counts`: 138,609 sheet rows read, 17,155 distinct symbols.
- `dropped`: "Excel date corruption", 209 rows (examples `2001-03-01`,
  `2001-09-01`, `2001-12-01`, `2002-03-01`, `2002-09-01`).
- `burden`: 138,400 rows for 17,130 genes (138,609 - 209).
- Parity with the atlas's committed `mirrors/burden.tsv`: each of its 1,438
  `PMID:42230622` rows matched exactly one package row on (gene, stratum,
  consequence, `maf_max`), where the package symbol was mapped to an HGNC id
  through the approved `symbol` column of the atlas's `mirrors/genes.tsv`; all
  other columns were byte-identical in 1,438 of 1,438.

## Tests

`pytest sources/audain2026_sd3`: a synthetic three-sheet workbook; no network,
no committed mirror.
