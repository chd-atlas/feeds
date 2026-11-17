# mane_select

MANE v1.5 Select transcripts and their exons, cut +-50 bp from GRCh38.p14, as the
`mane_select` package: tables `mane_transcripts` v1 (one row per gene: transcript and
protein sequence, `cds_start`, `cds_end`) and `mane_exons` v1 (one row per exon: transcript
and genomic coordinates, and the forward-strand genome window of the exon plus 50 bases
either side). `ids: hgnc_id`, `coverage: every_gene`, version `1.5`. Licence: NCBI states no
restrictions on use or distribution; recorded as `permitted_with_attribution` (see
`source.yaml`).

The package is the "frame" the atlas reads variants against. The atlas checks it itself
(MAN001: the exon cores, joined, are the transcript; MAN002: the coding sequence translates
to the protein; MAN003/4 against its gene registry). This recipe runs MAN001 and MAN002
before writing, with the atlas's two MAN002 allowances (`U` against `TGA`, and `M` at
residue 1 against a codon other than `ATG`), so a transcript that would fail them is a
`dropped` entry with its reason instead of an atlas error.

## Inputs (all five pinned by sha256; any other bytes are refused)

| file | sha256 |
|---|---|
| `MANE.GRCh38.v1.5.summary.txt.gz` | `d10ace27...` |
| `MANE.GRCh38.v1.5.refseq_genomic.gtf.gz` | `59e5ff1c...` |
| `MANE.GRCh38.v1.5.refseq_rna.fna.gz` | `7280512e...` |
| `MANE.GRCh38.v1.5.refseq_protein.faa.gz` | `b12307f9...` |
| `GCF_000001405.40_GRCh38.p14_genomic.fna.gz` (972,898,531 bytes) | `11912a45...` |

The MANE files are under
`https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/`; the genome is NCBI's
`GCF_000001405.40_GRCh38.p14` FTP directory. All five were downloaded on 2026-11-17.
`README_versions.txt` of the release: MANE 1.5, NCBI RefSeq Annotation Release
`GCF_000001405.40-RS_2025_08`, Ensembl release 116.

## Build

```bash
mkdir -p local/mane   # local/ is git-ignored; download the five files above into it
```

```python
from datetime import date
from pathlib import Path
from sources.mane_select import recipe

L = Path("local/mane")
recipe.build(
    L / recipe.SUMMARY,
    L / recipe.GTF,
    L / recipe.RNA,
    L / recipe.PROTEIN,
    L / recipe.GENOME_NAME,
    retrieved_on=date(2026, 11, 17),
    revision=1,
    built_by={"host": "<machine>"},  # recipe_commit is read from the clean checkout
    out=Path("dist/mane_select"),
)
```

See Conventions in the root README for `revision`. The genome is streamed once, one
chromosome at a time (only the chromosomes the summary names are kept in memory).

## How the rows are made

* **Rows**: summary rows whose `MANE_status` is exactly `MANE Select`. The GTF is read for
  `exon` and `CDS` features that carry `tag "MANE Select";`, are on the chromosome and
  strand the summary names, and whose `transcript_id` is a kept summary transcript.
* **Chromosome names**: `NC_0000NN` 1-22 is `1`..`22`, `NC_000023` is `X`, `NC_000024` is
  `Y`, `NC_012920` is `MT`.
* **Exons** are numbered in transcript order (ascending on `+`, descending on `-`; on `-`
  exon 1 holds the highest genomic coordinates), and `transcript_start` / `transcript_end`
  are cumulated from 1.
* **`cds_start` / `cds_end`**: the transcript position of the first `CDS` base, and of the
  last `CDS` base plus 3 (the GTF `CDS` excludes the stop codon). `start_codon` /
  `stop_codon` features are not used: measured, the GTF has 19,400 `start_codon` and
  19,388 `stop_codon` features for 19,437 transcripts, so they are neither one per
  transcript nor always whole (a codon can be split by an intron).
* **Windows** are `[genomic_start - 50, genomic_end + 50]` clipped to the chromosome;
  `window_start` / `window_end` are the truth, never a fixed 50. Soft-masked bases are
  upper-cased.
* **Ambiguity codes** in a flank are written as `N`; one inside an exon core would drop the
  transcript (the transcript sequence is `ACGT` only).
* More than 50 failures (`MAX_FAILURES`) raise `PackageRefused`: that many would mean the
  parser is wrong, not MANE. Only genuine self-check failures count; the three by-design
  exclusions below are listed in `dropped` but do not count.

## Measurements (real build, 2026-11-17, workstation, pinned inputs)

Method: `recipe.read_rows` over `local/mane/`, and `recipe.build` into a scratch directory
(`/usr/bin/time -v` for time and peak memory).

* 19,363 `MANE Select` rows read (`input_counts.rows`); **19,179 transcripts and 201,933
  exons kept**; 184 listed in `dropped`:
  * 70: non-coding (`NR_`) MANE Select transcripts, which have no protein;
  * 52: `NM_` transcripts of genes with a blank `HGNC_ID` in the summary (`LOC...` genes);
  * 61: transcripts the summary places on a non-primary accession (`NT_`, `NW_`);
  * 1: `HGNC:42967` `NM_001424184.1` (TMEM247), MAN002 `cds length 659 is not a multiple of
    3`. Its GTF `CDS` features are annotated `exception "annotated by transcript or
    proteomic data"` and their lengths sum to 656 bases plus the 3-base stop. The exons do
    reproduce the RNA sequence (MAN001 passes). It is dropped, never repaired.
* Every one of the atlas's 186 registry genes (`mirrors/genes.tsv`, atlas 2850343) is in the
  package, and every registry `mane_select` equals the package's transcript (MAN003 and MAN004
  would report nothing).
* Ambiguity codes: **0 flank bases carry an IUPAC code other than `N`** (so none were
  rewritten) and **1 flank base is already `N`** in GRCh38.p14; **0 transcripts were dropped
  for a non-ACGT base in an exon core**. No MANE RNA sequence differed from the
  genome-assembled one (0 MAN001 failures among the 19,180 candidates).
* `read_rows` took 68 s and `build` took 208 s (3 min 28 s) in all, peak memory
  1.1 GB (it hashes the five inputs and validates every cell against the vendored schemas).
* Package: `mane_transcripts` 19,179 rows over 19,179 genes (28,964,465 bytes gzipped),
  `mane_exons` 201,933 rows over 19,179 genes (31,820,381 bytes), `manifest.json` 4,926 bytes.
* Cross-check with the atlas's own code, read-only: the two tables, un-gzipped into a scratch
  `mirrors/` beside the atlas's `mirrors/genes.tsv` (atlas 2850343), went through the atlas's
  `read_mane_frame` and `validate_mane`: 19,179 transcripts read, 0 read issues, 0
  MAN001-MAN004 issues.

## Drift

`scripts/drift.py` probes `current/README_versions.txt` on NCBI's FTP (`FilesSource`;
sha256 and size, baseline `drift.json`, captured 2026-11-17 and equal to the release-1.5
file). It reports the day MANE publishes a newer release; it does not re-read the pinned
files.

## Tests

`pytest sources/mane_select`: a synthetic genome and four synthetic MANE files, no network.
`test_the_real_v1_5_files_build_when_present` needs `local/mane/` and `FEEDS_REAL=1` (about
70 s) and skips otherwise; it asserts the counts above, TBX5 `NM_181486.4` on `-` with
`cds_start` 556 and residue 279 `R`, and NR2F2 `NM_021005.4` with `cds_start` 1407.
