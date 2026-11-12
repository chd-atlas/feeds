# audain2021_s6

Audain et al. 2021, PLoS Genetics 17(7) e1009679: S6 Table Sheet A (gene-based
case/control CNV-deletion permutation testing, every gene PLINK tested) as one
`burden` table, one row per S6 gene, with `ids: symbol` and
`coverage: every_gene`. S8 Table (the 528 significant genes, with the CHD-only
and TAA-only contributions) is the package's second input: its "only CHD"
`case` is transcribed into `n_chd_case_carriers` for the genes it lists. The
atlas's `chd-atlas ingest` resolves the symbols (approved, then alias, then
previous; owner decision F22) and keeps its registry's genes, which recovers
CIROZ from `C1orf127`.

`source.yaml` is the atlas's `audain2021_s6` entry (name, URL, CC0-1.0,
`redistribution: permitted`, the licence evidence). The atlas's `audain2021_s8`
entry stays hand-written there; the manifest records S8 as its second input
(URL, sha256) so the atlas can compare them.

## Fetch and build

Download both files by hand into `local/audain2021/` (ignored by git): S6 from
the URL in `source.yaml` as `s011.xlsx`, S8 from the `s013` URL in
`recipe.S8_URL` as `s013.xlsx`. `build` refuses any file whose sha256 is not
`e8bd50c0a454d3f5563a6c1ee52c62206d98840cc696c3c8db4c0e0eccb63e62` (S6) or
`d070834994949311b3cff7a8a8470944c6ddd14c4644351062b0243c57db0ebc` (S8). The
version is `e8bd50c0` (S6's prefix) at revision 1.

```python
from pathlib import Path
from sources.audain2021_s6 import recipe

recipe.build(
    Path("local/audain2021/s011.xlsx"),
    Path("local/audain2021/s013.xlsx"),
    revision=1,
    built_by={"host": "<machine>"},  # recipe_commit is read from the clean checkout
    out=Path("dist/audain2021_s6"),
)
```

## Rules

Moved from the atlas's `scripts/convert_audain2021.py`, whose docstring (kept in
`recipe.py`) has the measurements behind each rule.

- **`gene` is S6's `symbol` verbatim**; no registry join here.
- S8's "CHD+TAA" `case` must equal S6's `case` for every gene S8 lists, or the
  build is refused; so is an S8 band or `case` column that is not found exactly
  once from the two header rows, a gene S8 lists twice, and a gene S6 lists
  twice.
- A blank S8 "only CHD" cell stays blank (not 0); genes S8 does not list carry
  none.
- Every row carries 6,746 cases and 14,024 controls (the filtered set the
  Methods ran the test on) and the `method_note` saying the case set includes a
  thoracic aortic aneurysm series.
- `input_counts.rows` is the S6 rows read and `genes` the distinct symbols.
  Nothing is dropped.

## Measurements

Measured 2026-11-12 by `recipe.build` on the pinned workbooks, revision 1 (run
from a clean checkout at the commit that added the recipe; the package was
built into a scratch directory, not `dist/`):

- `input_counts`: 11,515 S6 rows read, 11,515 distinct symbols; `dropped`: none.
- `burden`: 11,515 rows; 521 carry an S8 `n_chd_case_carriers` (S8 lists 528
  genes and leaves the cell blank for 7).
- Parity with the atlas's committed `mirrors/burden.tsv`: each of its 121
  `PMID:34324492` rows matched exactly one package row, where the package symbol
  was mapped to an HGNC id through the approved `symbol` column of the atlas's
  `mirrors/genes.tsv`; all other columns were byte-identical in 121 of 121 (0
  cells differing), and 12 of the 121 carry an S8 count.

## Tests

`pytest sources/audain2021_s6`: synthetic S6 and S8 workbooks written with
`tests/xlsx_fixture.py`; no network, no mirror.
