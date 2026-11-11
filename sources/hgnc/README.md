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

1. On the login node, stage a dated monthly archive (compute nodes have no
   internet). The live file's bytes change under the same URL, so a manifest
   naming it would carry a sha256 its URL no longer returns; an archive file is
   fixed to one date. List the archive, pick a date, download it:

   ```bash
   BASE=https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv
   curl -s "https://storage.googleapis.com/storage/v1/b/public-download-files/o?prefix=hgnc/archive/archive/monthly/tsv/&fields=items(name)" \
     | grep -o 'hgnc_complete_set_[0-9-]*\.txt' | sort | tail -3
   mkdir -p "$RAW_DIR" && curl -sSL -o "$RAW_DIR/hgnc_complete_set_<YYYY-MM-DD>.txt" \
     "$BASE/hgnc_complete_set_<YYYY-MM-DD>.txt"
   ```

   `release` is that date and `raw_url` is that URL; `build` refuses a `release`
   that is not `YYYY-MM-DD` or that differs from the date in the file's name.

2. Submit `hpc/hgnc_build.sbatch` (`sbatch -p <the cluster's CPU partition>`; the
   file names no site partition) with `RAW_FILE` (or a `RAW_DIR` holding only that
   file), `OUT_DIR`, `HVANTK_REPO` and `HVANTK_PYTHON` set. It needs a `java` on
   PATH (it refuses without one and warns if it is not 11) and sets hvantk's
   Spark settings itself. It runs hvantk's `hgnc:lookup` (`--skip-download`) and
   exports `hgnc_lookup.tsv.gz`. The job log prints the hvantk checkout and
   commit, the host and the job id: put them in `built_by`. The sbatch has not
   been run; its time and memory are placeholders until a first run's `seff`.

   Local fallback, measured by the review of this recipe: hvantk's own poetry
   environment with Java 21 and hail 0.2.137 ran the same command in 24 s. This
   is a fallback the owner may allow or refuse.

3. Build the package into an empty directory:

   ```python
   from datetime import date
   from pathlib import Path
   from sources.hgnc import recipe

   recipe.build(
       Path("hgnc_lookup.tsv.gz"),
       Path("hgnc_complete_set_2026-10-06.txt"),
       raw_url=".../monthly/tsv/hgnc_complete_set_2026-10-06.txt",  # the full URL
       retrieved_on=date(2026, 11, 11),
       release="2026-10-06",
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
measurements (its RBM10 and MANE figures were taken on 2026-10-25 over the 2026-10-16 complete set; its `prev_symbol` figures over the 2026-07-31 set).

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
  refuses. Measured 2026-11-11: the live set (45,387 rows) and the 2026-10-06
  archive (45,233 rows) hold only Approved rows, so `dropped` is `[]` on both and
  this guard is exercised by fixtures alone.
- **Columns are read by name**, as strings with no quoting
  (`quote_char=None`), so no cell is retyped. A missing column is named. The
  complete set is read with polars' default quoting, because its multi-value
  cells are quoted.
- Rows are sorted on the `hgnc_id` string, as `genes` v2 declares.

## Tests

`pytest sources/hgnc`: five synthetic genes, no network. On real data (review,
2026-11-11), hvantk at 4c90a241 run on the live set (45,387 rows) and on the
2026-10-06 archive (45,233 rows), through this recipe and then the atlas's own
ingest, gave a `genes.tsv` byte-identical to `scripts/convert_hgnc.py`'s.
