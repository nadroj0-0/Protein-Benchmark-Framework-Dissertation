#!/usr/bin/env bash
# Compare preserved narrow-policy rebuilds with the accepted benchmark CSVs.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=24G
#$ -l tscratch=12G
#$ -l scratch0free=20G
#$ -l h_rt=8:0:0
#$ -pe smp 1
#$ -t 1-2
#$ -tc 2
#$ -N go_rel_cmp
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

TASK_ID="${SGE_TASK_ID:-}"
case "$TASK_ID" in
  1)
    PROFILE="2025_01_to_2026_02_supervisor"
    BROAD_DIR="/SAN/bioinf/bmpfp/benchmarks/contemporary/2025_01_to_2026_02_supervisor"
    PRESERVED_RUN="/SAN/bioinf/bmpfp/diagnostics/go_relationship_policy/2025_01_to_2026_02_supervisor/7138009_20260806_103030"
    ;;
  2)
    PROFILE="2025_01_to_2026_02_supervisor_nk_lk"
    BROAD_DIR="/SAN/bioinf/bmpfp/benchmarks/contemporary/2025_01_to_2026_02_supervisor_nk_lk/7128984_20260801_205146/outputs"
    PRESERVED_RUN="/SAN/bioinf/bmpfp/diagnostics/go_relationship_policy/2025_01_to_2026_02_supervisor_nk_lk/7138010_20260806_103030"
    ;;
  *) die "SGE_TASK_ID must be 1 or 2" ;;
esac

NARROW_DIR="$PRESERVED_RUN/outputs"
OBO_FILE="${OBO_FILE:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2025-02-06/go-basic.obo}"
RESULTS_ROOT="${RESULTS_ROOT:-$PRESERVED_RUN/comparison_retries}"
JOB_TOKEN="${JOB_ID:-manual}_${TASK_ID}"
RUN_TAG="${JOB_TOKEN}_$(date -u +%Y%m%dT%H%M%SZ)"
WORK="/scratch0/go_relationship_policy_compare_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
STAGED_BROAD="$WORK/broad"
STAGED_NARROW="$WORK/narrow"
STAGED_OBO="$WORK/go-basic.obo"
SCRATCH_OUTPUT="$WORK/output"
FINAL_OUTPUT="$RESULTS_ROOT/$RUN_TAG"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
WORK_OWNED=0

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  if [[ "$WORK_OWNED" == 1 && "$WORK" == /scratch0/go_relationship_policy_compare_* && ! -L "$WORK" ]]; then
    cd "$HOME"
    rm -rf -- "$WORK"
  fi
  exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

copy_with_retry() {
  local source="$1" destination="$2" attempt
  for attempt in 1 2 3 4 5 6; do
    rm -f -- "${destination}.partial" "$destination"
    if cp -p "$source" "${destination}.partial" && mv "${destination}.partial" "$destination"; then
      return 0
    fi
    echo "Copy failed ($attempt/6), retrying: $source" >&2
    sleep 10
  done
  return 1
}

[[ -d "$BROAD_DIR" ]] || die "Accepted benchmark directory is missing: $BROAD_DIR"
[[ -d "$NARROW_DIR" ]] || die "Preserved narrow rebuild is missing: $NARROW_DIR"
[[ -f "$OBO_FILE" ]] || die "Frozen ontology is missing: $OBO_FILE"
[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
[[ ! -e "$FINAL_OUTPUT" ]] || die "Output already exists: $FINAL_OUTPUT"
mkdir -p "$STAGED_BROAD" "$STAGED_NARROW" "$RESULTS_ROOT"
WORK_OWNED=1

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a framework checkout"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"

echo "Profile          : $PROFILE"
echo "Accepted CSVs    : $BROAD_DIR"
echo "Narrow CSVs      : $NARROW_DIR"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "Final output     : $FINAL_OUTPUT"

for aspect in bp cc mf; do
  for split in training validation test; do
    name="${aspect}-${split}.csv"
    copy_with_retry "$BROAD_DIR/$name" "$STAGED_BROAD/$name" || die "Could not stage $BROAD_DIR/$name"
    copy_with_retry "$NARROW_DIR/$name" "$STAGED_NARROW/$name" || die "Could not stage $NARROW_DIR/$name"
  done
done
copy_with_retry "$OBO_FILE" "$STAGED_OBO" || die "Could not stage ontology"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

"$PYTHON_BIN" scripts/diagnostics/compare_go_relationship_policy.py \
  --broad-label accepted-all-relationships \
  --broad-dir "$STAGED_BROAD" \
  --narrow-label diagnostic-is-a-plus-part-of \
  --narrow-dir "$STAGED_NARROW" \
  --obo-file "$STAGED_OBO" \
  --output-dir "$SCRATCH_OUTPUT"

STAGING="${FINAL_OUTPUT}.staging-${JOB_TOKEN}"
cp -a "$SCRATCH_OUTPUT" "$STAGING"
mv "$STAGING" "$FINAL_OUTPUT"
echo "Completed relationship-policy comparison: $FINAL_OUTPUT"
