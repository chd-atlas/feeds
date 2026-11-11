# sierant2025_ds6

Sierant et al. 2025, PNAS 122(13) e2420343122, Dataset S6: the gene-level
meta-analysis over all 248 prespecified MIPseq panel genes of the Pediatric
Cardiac Genomics Consortium. Published as two tables, `burden` (6 rows per gene:
three consequence classes x a de novo arm and a case-control arm) and
`burden_combined` (3 per gene: the study's own JL-FDR combination of the two
arms), with `ids: symbol` and `coverage: every_gene`. The atlas's
`chd-atlas ingest` resolves the symbols to HGNC ids and keeps its registry's
genes.

`source.yaml` holds the name, the pnas.org URL, the licence (CC-BY-4.0,
`redistribution: permitted_with_attribution`) and the licence evidence, as the
atlas's `mirrors/sources.yaml` records them (confirmed there on 2026-08-05 three
ways: the article's permissions statement, Crossref's licence record, NCBI's OA
service returning `license="CC BY"` for PMC12002227).

## Fetch and build

```python
from pathlib import Path
from sources.sierant2025_ds6 import fetch, recipe

workbook = fetch.fetch(Path("data"))  # standard library only; the server is slow (minutes)
recipe.build(
    workbook,
    built_by={"recipe_commit": "<40 hex>", "host": "<machine>"},
    out=Path("dist/sierant2025_ds6"),
)
```

`fetch` downloads the Europe PMC open-access supplementary bundle for PMC12002227
(pnas.org answers 403 to a scripted request, per the atlas's `sources.yaml`),
reads only `pnas.2420343122.sd06.xlsx` from it, refuses a member path that
escapes the destination, and refuses (writing nothing) a file whose sha256 is not
`d8e9f7777b329c810adde8b9f4a0f91a4d1a6211b760daeb9329ae6f36606daa`. `build`
refuses any workbook with another sha256 too: the column map is positional and
was written against that file. The version is therefore `d8e9f777`
(`-r2`, ... after a recipe change that keeps the same input).

## Rules

Moved from the atlas's `scripts/convert_sierant2025.py`, whose docstring (kept in
`recipe.py`, attributed) has the measurements behind each rule.

- **`_check_layout` and `_check_bands`** read the workbook's header row and the
  row naming each test, and refuse a workbook where a column is not where the
  map says (a column inserted would relabel one consequence's statistic as
  another's, with every number still plausible).
- **`gene` is the panel symbol verbatim.** A symbol shaped like an Excel date
  (`2001-03-01`, `1-Mar`) is published in `dropped` as "Excel date corruption"
  with a count and up to five examples, never guessed.
- **A symbol on two panel rows refuses the build.**
- **`count_unit` per arm**: de novo rows count `de_novo_mutations` against
  3,887 trios; case-control rows count `alleles` against a per-gene allele total.
- **The effect is the study's printed value**; `Inf` becomes an empty effect with
  `effect_bound: unbounded_above`.
- **Rows are sorted** on the schema's sort key (gene, study, stratum, lesion
  group, variant class, consequence class, origin, `maf_max`), empty smallest.
- `input_counts.rows` is the panel rows read (dropped ones included) and
  `genes` the distinct panel symbols.

## Measurements

Measured 2026-11-11 on the workbook `fetch.fetch` downloaded from Europe PMC
(sha256 `d8e9f777...` as pinned; 233,828 bytes):

- 248 panel rows, 248 distinct symbols; **0 Excel-date symbols**, so `dropped`
  is `[]` on the real workbook and the drop is exercised by the tests alone.
  `burden` 1,488 rows (248 x 6) and `burden_combined` 744 rows (248 x 3).
- Through the atlas's own `ingest_package` (atlas `dev` at bd185f4, validation on,
  against a stand-in HGNC package built from the atlas's committed
  `mirrors/genes.tsv`) the resulting `burden.tsv` and `burden_combined.tsv` equal
  the committed mirrors byte for byte **except for one gene**: panel symbol
  `C1orf127` is the previous symbol of CIROZ (`HGNC:26730`, registry gene), which
  the atlas's resolver reaches (approved, then alias, then previous symbol) and
  `scripts/convert_sierant2025.py` (approved symbols only) did not. That adds 6
  `burden` and 3 `burden_combined` rows (252 and 126 rows for 42 genes, against
  the 246 and 123 for 41 the committed mirrors hold). Validation: 0 errors,
  3 warnings (two SCP003, one SCP007), the atlas's baseline.

## Tests

`pytest sources/sierant2025_ds6`: a synthetic workbook written by the same
writer the atlas's converter test uses, and a local `file://` bundle for `fetch`.
No network, no committed mirror.
