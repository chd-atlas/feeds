#!/bin/bash
#SBATCH --job-name=feeds-hgnc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=0-00:30:00
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
# Keep the time at or under the partition's MaxTime (2 h on the atlas's cluster): that site sets
# EnforcePartLimits=NO, so a longer request is accepted and then never starts
# (measured 2026-11-02, atlas scripts/hpc/farah_celltypes.sbatch).
#
# Compute nodes have no git and no internet: stage the dated archive file on the
# login node first (README.md step 1).
# Submit with `sbatch -p <the cluster CPU partition> ...`: this file is public and
# names no site partition. Required environment: OUT_DIR, HVANTK_REPO (a flat
# hvantk checkout at the commit to record), HVANTK_PYTHON, and the raw file as
# RAW_FILE (the dated archive, hgnc_complete_set_<YYYY-MM-DD>.txt) or RAW_DIR
# holding exactly one *.txt or *.tsv (hvantk's own rule).
#
# Hail needs a Java on PATH. The job refuses without one and warns when it is not
# 11 (hvantk's image carries 11; hail 0.2.137 was measured running under Java 21
# locally, so the cluster's version is decided by the first run). The Spark
# settings below follow hvantk's skills/hpc-hail-job template: local[cpus], a
# driver capped at 80% of the memory request, node-local scratch.
#
# `fields` names the columns the export keeps. At hvantk 8683258d the export flattens
# the struct-array column `lsdb` into `lsdb.*` array fields rather than refusing it; the
# selection stays so the package input carries only what recipe.py reads. The list is
# exactly HVANTK_FIELDS in recipe.py; a test compares the two.

set -euo pipefail

HPC_ENV_MODULE="${HPC_ENV_MODULE:-hpc-env/13.1}"
RAW_DIR="${RAW_DIR:-}"
RAW_FILE="${RAW_FILE:-}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to a writable output directory}"
HVANTK_REPO="${HVANTK_REPO:?Set HVANTK_REPO to the hvantk checkout}"
HVANTK_PYTHON="${HVANTK_PYTHON:?Set HVANTK_PYTHON to the hvantk environment python}"
FIELDS="hgnc_id,gene_symbol,gene_name,alias_symbols,prev_symbols,ensembl_gene_id,entrez_id,location,uniprot_ids,mane_select"
OUTPUT="$OUT_DIR/hgnc_lookup.tsv.gz"

if [[ -z "$RAW_FILE" ]]; then
  [[ -n "$RAW_DIR" ]] || { echo "ERROR: set RAW_FILE or RAW_DIR" >&2; exit 1; }
  shopt -s nullglob
  candidates=("$RAW_DIR"/*.txt "$RAW_DIR"/*.tsv)
  shopt -u nullglob
  [[ ${#candidates[@]} -eq 1 ]] \
    || { echo "ERROR: RAW_DIR must hold exactly one *.txt or *.tsv; found ${#candidates[@]}" >&2; exit 1; }
  RAW_FILE="${candidates[0]}"
fi
[[ -r "$RAW_FILE" ]] || { echo "ERROR: not readable: $RAW_FILE" >&2; exit 1; }
# hvantk reads the single *.txt/*.tsv in --raw-dir, so the directory must hold only this file.
RAW_DIR="$(dirname "$RAW_FILE")"
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
command -v java >/dev/null 2>&1 \
  || { echo "ERROR: no java on PATH; hail cannot start (see the header)" >&2; exit 1; }
java -version 2>&1 | sed 's/^/java: /'
if ! java -version 2>&1 | head -1 | grep -q '"11\.'; then
  echo "WARNING: java is not version 11 (hvantk's image carries 11); continuing" >&2
fi
DRIVER_MB=$(( ${SLURM_MEM_PER_NODE:-12288} * 8 / 10 ))
export PYSPARK_SUBMIT_ARGS="--master local[${SLURM_CPUS_PER_TASK:-4}] --driver-memory ${DRIVER_MB}m --conf spark.local.dir=${TMPDIR:-/tmp}/feeds-hgnc-spark-${SLURM_JOB_ID:-0} pyspark-shell"
export JAVA_TOOL_OPTIONS="-Duser.language=en -Duser.country=US"
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

export HVANTK_BUILDER_COMMIT="$HVANTK_HEAD"

echo "== feeds hgnc =="
echo "date:         $(date -Iseconds)"
echo "host:         $(hostname)"
echo "job:          ${SLURM_JOB_ID:-<none>}"
echo "python:       $("$HVANTK_PYTHON" -V 2>&1)"
echo "hvantk repo:  $HVANTK_REPO"
echo "hvantk HEAD:  $HVANTK_HEAD"
echo "raw name:     $(basename "$RAW_FILE")"
echo "raw sha256:   $(sha256sum "$RAW_FILE" | cut -d' ' -f1)"
echo "spark args:   $PYSPARK_SUBMIT_ARGS"
echo "fields:       $FIELDS"
echo "out:          $OUTPUT"

PYTHONPATH="$HVANTK_REPO" "$HVANTK_PYTHON" -m hvantk reprocess hgnc:lookup \
  --skip-download --no-check-drift --raw-dir "$RAW_DIR" \
  --plugin-arg "fields=$FIELDS" \
  --output "$OUTPUT"

echo "-- output --"
zcat "$OUTPUT" | wc -l
sha256sum "$OUTPUT"
echo "done: $(date -Iseconds)"
