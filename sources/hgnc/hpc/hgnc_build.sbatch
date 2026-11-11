#!/bin/bash
#SBATCH --job-name=feeds-hgnc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=0-00:30:00
#SBATCH --partition=all_cpu.p
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --mail-type=END,FAIL
# --mail-user is deliberately NOT a #SBATCH line: Slurm reads directives
# literally and does not expand env vars, so a `--mail-user=$VAR` directive
# silently mails nobody. Pass it on the sbatch command line.
#
# Turns HGNC's complete set into hvantk's `hgnc:lookup` table, exported as
# gzipped TSV for sources/hgnc/recipe.py. Not run yet: the time and memory
# above are placeholders to replace with `seff` from the first real run.
# Keep the time at or under all_cpu.p's MaxTime (2 h): that site sets
# EnforcePartLimits=NO, so a longer request is accepted and then never starts
# (measured 2026-11-02, atlas scripts/hpc/farah_celltypes.sbatch).
#
# Compute nodes have no git and no internet. Before submitting, on the login node:
#   mkdir -p "$RAW_DIR" && curl -sSL -o "$RAW_DIR/hgnc_complete_set.txt" \
#     https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
# Required environment: RAW_DIR (holds hgnc_complete_set.txt), OUT_DIR, HVANTK_REPO
# (a flat hvantk checkout at the commit to record), HVANTK_PYTHON.
#
# `fields` names the columns the export keeps: the TSV export refuses `lsdb`
# (array<struct>), which the builder otherwise keeps. The list is exactly
# HVANTK_FIELDS in recipe.py; a test compares the two.

set -euo pipefail

HPC_ENV_MODULE="${HPC_ENV_MODULE:-hpc-env/13.1}"
RAW_DIR="${RAW_DIR:?Set RAW_DIR to the directory holding hgnc_complete_set.txt}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to a writable output directory}"
HVANTK_REPO="${HVANTK_REPO:?Set HVANTK_REPO to the hvantk checkout}"
HVANTK_PYTHON="${HVANTK_PYTHON:?Set HVANTK_PYTHON to the hvantk environment python}"
FIELDS="hgnc_id,gene_symbol,gene_name,alias_symbols,prev_symbols,ensembl_gene_id,entrez_id,location,uniprot_ids,mane_select"
OUTPUT="$OUT_DIR/hgnc_lookup.tsv.gz"

for f in "$RAW_DIR/hgnc_complete_set.txt"; do
  [[ -r "$f" ]] || { echo "ERROR: not readable: $f" >&2; exit 1; }
done
[[ -x "$HVANTK_PYTHON" ]] || { echo "ERROR: HVANTK_PYTHON not executable: $HVANTK_PYTHON" >&2; exit 1; }
mkdir -p "$OUT_DIR"
: > "$OUT_DIR/.writetest" 2>/dev/null || { echo "ERROR: $OUT_DIR not writable" >&2; exit 1; }
rm -f "$OUT_DIR/.writetest"

if command -v ml >/dev/null 2>&1; then
  ml purge || true
  ml "$HPC_ENV_MODULE" || echo "WARNING: could not load $HPC_ENV_MODULE" >&2
fi
# A cluster module must not shadow the venv's own packages.
unset PYTHONPATH
PYTHONPATH="$HVANTK_REPO" "$HVANTK_PYTHON" -c "import hvantk" \
  || { echo "ERROR: cannot import hvantk from $HVANTK_REPO" >&2; exit 1; }

# The compute nodes have no git: read the commit from the checkout's own files.
hvantk_head() {
  local git_dir="$HVANTK_REPO/.git" head ref
  head="$(<"$git_dir/HEAD")"
  if [[ "$head" != ref:* ]]; then echo "$head"; return; fi
  ref="${head#ref: }"
  if [[ -r "$git_dir/$ref" ]]; then
    cat "$git_dir/$ref"
  else
    awk -v r="$ref" '$2 == r { print $1 }' "$git_dir/packed-refs"
  fi
}
HVANTK_HEAD="$(hvantk_head)"
[[ "$HVANTK_HEAD" =~ ^[0-9a-f]{40}$ ]] \
  || { echo "ERROR: could not read the hvantk commit from $HVANTK_REPO" >&2; exit 1; }

echo "== feeds hgnc =="
echo "date:         $(date -Iseconds)"
echo "host:         $(hostname)"
echo "job:          ${SLURM_JOB_ID:-<none>}"
echo "python:       $("$HVANTK_PYTHON" -V 2>&1)"
echo "hvantk repo:  $HVANTK_REPO"
echo "hvantk HEAD:  $HVANTK_HEAD"
echo "raw:          $RAW_DIR/hgnc_complete_set.txt"
echo "raw sha256:   $(sha256sum "$RAW_DIR/hgnc_complete_set.txt" | cut -d' ' -f1)"
echo "fields:       $FIELDS"
echo "out:          $OUTPUT"

PYTHONPATH="$HVANTK_REPO" "$HVANTK_PYTHON" -m hvantk reprocess hgnc:lookup \
  --skip-download --raw-dir "$RAW_DIR" \
  --plugin-arg "fields=$FIELDS" \
  --output "$OUTPUT"

echo "-- output --"
zcat "$OUTPUT" | wc -l
sha256sum "$OUTPUT"
echo "done: $(date -Iseconds)"
