# expression_atlas_e_mtab_6814

EMBL-EBI Expression Atlas experiment E-MTAB-6814, the reprocessing of
Cardoso-Moreira et al. 2019 (human developmental transcriptomes, seven organs):
the per-assay-group TPM matrix and its configuration. Published as two sharded
tables (shard `E-MTAB-6814`), with `ids: ensembl_gene_id` and
`coverage: every_gene`:

- `profiles`: every one of the matrix's 58,735 genes, one row per measured
  (tissue, stage), by Ensembl id. The atlas's `chd-atlas ingest` resolves the ids
  against its genes table and keeps its registry's genes.
- `profile_quantiles`: 101 breakpoints per (tissue, stage), computed over the
  whole matrix, with `n_genes` naming the population.

The values are the Atlas's reprocessing (iRAP/HISAT2/FeatureCounts, TPM), not the
paper's RPKM; the docstring in `recipe.py` has the detail and the rules.

## Where it runs, and why

On the workstation, not the cluster. The atlas's old `scripts/hpc/cardoso_mirrors.sbatch`
was an 11-second standard-library pass over a public download; there is nothing
for a cluster to do. It is not ported. `write_package` holds the 6,989,465
profile rows in memory and checks each cell in Python, measured at 4.49 GB peak RSS on the workstation (see Measurements).

## Build

Download both files from the EBI FTP (`MATRIX_URL`, `CONFIGURATION_URL` in
`recipe.py`) into `local/e-mtab-6814/`; `build` refuses any bytes but the pinned
sha256s (matrix `92484a97...`, version `92484a97`).

```python
from pathlib import Path
from sources.expression_atlas_e_mtab_6814 import recipe

recipe.build(
    Path("local/e-mtab-6814/E-MTAB-6814-tpms.tsv"),
    Path("local/e-mtab-6814/E-MTAB-6814-configuration.xml"),
    revision=1,
    built_by={"host": "workstation"},
    out=Path("dist/expression_atlas_e_mtab_6814"),
)
```

## Measurements

Measured 2026-11-12 by `recipe.build` on the pinned matrix and configuration,
revision 1, on the workstation, into a scratch directory (not `dist/`); the final
review's rebuild is the source of the time and memory:

- `input_counts`: 58,735 genes (each with at least one measured cell); `dropped`:
  none.
- `profiles`: 6,989,465 rows; `profile_quantiles`: 12,019 rows.
- Cost: about 3:00 wall (179.5 s) and 4.49 GB peak RSS (`resource.getrusage`
  `ru_maxrss`).
- Parity with the atlas's committed mirrors at atlas 4b1f8f6: the quantile shard
  is byte-identical to `mirrors/profile_quantiles`'s, and all 22,015 committed
  profile rows are present in the package.
