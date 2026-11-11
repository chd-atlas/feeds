# hgnc

Every approved HGNC gene, as the atlas's `genes` v2 table (`ids: hgnc_id`,
`coverage: every_gene`). The atlas's `chd-atlas ingest` writes only the rows of
its own registry genes from it, and resolves every other package's gene ids
against it (approved symbol, then alias, then previous symbol).

`source.yaml` holds the name, URL, licence (CC0-1.0, `redistribution:
permitted`) and the licence evidence, which quotes
<https://www.genenames.org/about/license/> as the atlas's `mirrors/sources.yaml`
recorded it on 2026-10-25.

## Build

1. On the login node, stage the complete set (compute nodes have no internet),
   and note the file's `last-modified` header: that date is the release.

   ```bash
   mkdir -p "$RAW_DIR"
   curl -sSLD - -o "$RAW_DIR/hgnc_complete_set.txt" \
     https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
   ```

2. Submit `hpc/hgnc_build.sbatch` with `RAW_DIR`, `OUT_DIR`, `HVANTK_REPO` and
   `HVANTK_PYTHON` set. It runs hvantk's `hgnc:lookup` (`--skip-download`) and
   exports `hgnc_lookup.tsv.gz`. The job log prints the hvantk checkout and
   commit, the host and the job id: put them in `built_by`. The sbatch has not
   been run; its time and memory are placeholders until a first run's `seff`.

3. Build the package into an empty directory:

   ```python
   from datetime import date
   from pathlib import Path
   from sources.hgnc import recipe

   recipe.build(
       Path("hgnc_lookup.tsv.gz"),
       Path("hgnc_complete_set.txt"),
       raw_url="https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
       retrieved_on=date(2026, 11, 11),
       release="2026-10-16",
       revision=1,
       built_by={
           "recipe_commit": "<40 hex>",
           "host": "<node>",
           "job": "<id>",
           "hvantk_release": "<hvantk commit>",
       },
       out=Path("dist/hgnc"),
   )
   ```

   The version is `release` (`-r2`, `-r3` ... after a recipe change that keeps the
   same input). The writer refuses before writing a byte.

## Rules, and why

Moved from the atlas's `scripts/convert_hgnc.py`, whose docstring has the
measurements (taken there on 2026-10-25 over the 2026-10-16 complete set).

- **`mane_select`** is HGNC's `ENST...|NM_...` pair. The column holds the RefSeq
  half, chosen by its `NM_` prefix rather than position, so a pair without one
  yields null and never an Ensembl id in a column read as RefSeq.
- **`uniprot`** is HGNC's `uniprot_ids` verbatim, a pipe-joined list included.
  Exception: RBM10 (`HGNC:9896`) is pinned to `P98175`. HGNC lists
  `P0DW28|P98175` alphabetically, not ranked; P98175 is RBM10_HUMAN (930 aa) and
  P0DW28 is MINAS-60 (130 aa), a microprotein encoded in an alternative reading
  frame of RBM10 transcripts (UniProt REST, 2026-10-25, recorded in `recipe.py`).
  A pin that HGNC no longer lists for the gene, or whose gene is absent from the
  table, refuses the build. This differs from the atlas converter, which refused
  every other multi-accession cell: the package writes it verbatim and the atlas
  decides at ingest (the atlas's decision 3).
- **Dropped rows.** hvantk keeps `Approved` genes only. The complete set's other
  rows are published in `dropped` by status, with count and up to five example
  ids; the table must hold exactly the complete set's Approved ids, or the build
  refuses.
- **Columns are read by name**, as strings with no quoting
  (`quote_char=None`), so no cell is retyped. A missing column is named. The
  complete set is read with polars' default quoting, because its multi-value
  cells are quoted.
- Rows are sorted on the `hgnc_id` string, as `genes` v2 declares.

## Tests

`pytest sources/hgnc`: five synthetic genes, no network. The expected rows were
also compared byte for byte with `scripts/convert_hgnc.py` run on a four-gene raw
file (the two outputs were identical).
