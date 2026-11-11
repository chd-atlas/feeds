# farah2024_hoc

Farah et al. 2024's developing human heart (UCSC Cell Browser dataset
`hoc/all-heart`, 142,946 cells), aggregated to one pseudobulk per heart, as the
atlas's `cell_type_profiles` table (shard `hoc`). A **registry slice**: `ids:
symbol`, and `coverage` is `{"registry_slice": {"atlas_commit", "genes"}}`.

Moved from the atlas's `scripts/convert_farah2024.py` and its three cluster
files (`scripts/hpc/farah_celltypes.sbatch`, `farah_extract.py`,
`submit_farah_celltypes.sh`, now `hpc/submit.sh`).

`source.yaml` holds what the atlas's `mirrors/sources.yaml` records: UCSC's
`hoc/desc.json` and `hoc/dataset.json` (fetched 2026-11-02) state **no licence**;
by an owner decision of 2026-11-02 the data are treated as attributable to Farah
et al. 2024, whose article (PMID:38480880) is CC BY 4.0, hence
`redistribution: permitted_with_attribution`. This is not a grant by the data's
submitter.

## Build

The 3.3 GB h5ad stays on the cluster; the recipe cannot hash it. The version is
the first 8 hex digits of the h5ad's sha256, which the job log prints.

1. On the login node, stage the inputs in one directory: the h5ad
   (`farah_2024.h5ad`), `hpc/farah_celltypes.sbatch`, `hpc/farah_extract.py`,
   `hpc/submit.sh`, and the atlas's gene registry at the commit the slice is cut
   at:

   ```bash
   git -C <atlas checkout> show <atlas_commit>:mirrors/genes.tsv > genes.tsv
   ```

2. Submit (compute nodes have no git and no internet; the files name no
   partition, no mail address or site path). `JOB_CODE_URL` is required: a GitHub
   tree permalink, at the full commit, to the `hpc/` scripts you staged (for a
   feeds checkout, `https://github.com/chd-atlas/feeds/tree/<40 hex>/sources/farah2024_hoc/hpc`).
   The job refuses anything else, prints it, and also prints the registry's
   sha256:

   ```bash
   CPU_PARTITION=<the cluster CPU partition> INPUT_H5AD=<h5ad> \
     HVANTK_REPO=<flat hvantk checkout> HVANTK_PYTHON=<hvantk python> \
     JOB_CODE_URL=<permalink> OUT_DIR=<dir> bash submit.sh
   ```

   or `sbatch -p <partition> farah_celltypes.sbatch` with the same variables and
   `REGISTRY` set. The job runs `hvantk expression summarize` twice (by Age,
   Rep, cell class and Region; then without Region) with `--min-cells 10`, and
   `farah_extract.py` flattens both into `$OUT_DIR/pseudobulk.tsv`. Measured on
   the first run (job 21367175, hvantk 96ceafb4, 2026-11-02, per the atlas's
   `sbatch` comments): 2 m 30 s, MaxRSS 6.39 GB, 81,993 rows. At hvantk
   a7ce8302 (job 21476247, 2026-11-11), 1 m 02 s and MaxRSS 12.10 GB, with
   the same `pseudobulk.tsv` bytes; a 12G request was killed for memory, so
   the sbatch asks for 48G. The time limit
   must stay at or under the partition's MaxTime (2 h on the atlas's cluster).

3. Read from the log the h5ad's digest (`input sha256: ...`), the registry's
   (`registry sha256: ...`), the code permalink (`job code url: ...`) and the
   hvantk commit (`hvantk HEAD: ...`). Then, from a **clean checkout of the
   `origin/dev` commit you will publish**, build the package into an empty
   directory. `build` reads the commit from git and refuses uncommitted changes
   under `src`, `sources`, `schemas`, `pyproject.toml` and `uv.lock`:

   ```python
   from datetime import date
   from pathlib import Path
   from sources.farah2024_hoc import recipe

   recipe.build(
       Path("pseudobulk.tsv"),
       Path("genes.tsv"),
       atlas_commit="<40 hex>",  # the commit genes.tsv was read at
       h5ad_sha256="<64 hex, from the job log>",
       registry_sha256="<64 hex, from the job log; must equal genes.tsv's>",
       pseudobulk_made_on=date(2026, 11, 2),  # the job log's `date:` line, or sacct
       job_code_url="<the job log's `job code url:` line>",
       revision=1,
       built_by={
           "host": "<node>",
           "job": "<id>",
           "hvantk_release": "<hvantk HEAD from the job log>",
       },
       out=Path("dist/farah2024_hoc"),
   )
   ```

   `build` refuses an `h5ad_sha256` or `registry_sha256` that is not 64 lowercase
   hex, a `registry_sha256` that is not the digest of the `genes.tsv` it is given,
   an `atlas_commit` that is not a 40-hex sha, a `job_code_url` that is not a
   `chd-atlas/feeds` or `chd-atlas/chd-atlas.github.io` tree permalink at a 40-hex
   commit, and a `built_by` without `hvantk_release`. The version is `19a1d496`
   for the h5ad of `sources.yaml`'s `sha256:19a1d4961c2f...` at revision 1 (see
   Conventions in the root README for `-r2`, ...).

`inputs` in the manifest holds two entries: the h5ad (the UCSC URL, the digest
the log printed, retrieved_on 2026-04-18, the h5ad's own build date) and the
pseudobulk (its sha256, computed by `build`; `retrieved_on` is `pseudobulk_made_on`, the
job's date). The pseudobulk file is a cluster-local intermediate and is not published:
only its sha256 is recorded. The schema requires a `url`, so it is `job_code_url`,
the permalink the job printed to the code that made it. That is not necessarily
`built_by["recipe_commit"]`, which names the code that ran `build`.

For the cached stand-in run (job 21367175, 2026-11-02, hvantk 96ceafb4), the code
that ran was the atlas's `scripts/hpc`, not feeds, so its `job_code_url` is
`https://github.com/chd-atlas/chd-atlas.github.io/tree/<sha>/scripts/hpc`, where
`<sha>` is atlas commit 6d5c08a expanded to its full 40 hex digits
(`git -C <atlas checkout> rev-parse 6d5c08a`).

## Rules, and why

- **`input_counts` describes the job's already-sliced pseudobulk.** The cut from
  the source's 33,538 genes to the registry's happens in the job
  (`hpc/farah_extract.py`), so the manifest shows only what survived it. Measured
  2026-11-11 with the real pseudobulk: `input_counts = {rows: 81,993, genes: 181}`
  and `dropped: []`; nothing in `build` is left out on real data, and the `dropped`
  entry is exercised by fixtures alone.
- **Per (gene, cell type, region)** the package carries the median across hearts
  of counts per million, quartiles when at least three hearts support them, the
  number of hearts, the cells pooled and the share of them expressing the gene.
  The CPM denominator is the heart group's whole library (summed over every gene
  the source measured), not the registry's. A heart is `{age}pcw-{replicate}`:
  two hearts of one age are two donors. Quartiles use R's type-7 estimator,
  copied from the atlas's `convert_cardoso_moreira.py`.
- **The slice is approved symbols only** (the atlas's Decision 8). The source
  names genes by symbol and has no Ensembl id, so an alias cannot be told from
  another gene's approved symbol. The atlas converter measured this on
  2026-10-28 against its `mirrors/genes.tsv`: `HGF` is SOS1's alias and
  hepatocyte growth factor's symbol, `PKD1` is PRKD1's alias and polycystin 1's,
  and `POP1` is POPDC1's alias and the ribonuclease P subunit's. A source symbol
  that is not a registry row's approved symbol is left out and counted in
  `dropped`. A registry gene whose approved symbol the source lacks has no rows:
  CRIPTO, because the source uses its previous symbol TDGF1.
- **The build stops** on a registry that is not the one the job sliced with (the
  `registry_sha256` check), on a source symbol with leading or trailing whitespace (named, never
  stripped), a heart counted twice in one group, a heart's group with fewer than
  10 cells (hvantk drops such a group at `--min-cells`, so one arriving means the
  job ran without it; checked on every row, in the slice or not), a registry
  listing a gene or symbol twice, and any malformed input row (wrong header or
  field count, blank symbol/heart/cell type, unknown region, a non-finite or
  negative count, a non-positive library, a non-integer or out-of-range cell
  count). The 10-cell and 3-heart values are copied from the atlas's
  `tables.py`; the atlas enforces them too (CTP009, CTP010).

## Tests

`pytest sources/farah2024_hoc`: synthetic input only, no network. Parity with
the atlas, measured on 2026-11-11 with the real pseudobulk (sha256
`3414fe0c...59d9`), the registry at atlas commit f094c28 (185 rows) and
`chd-atlas ingest`: `mirrors/cell_type_profiles/hoc.tsv` byte-identical to the
committed file (12,308 rows over 181 genes).
