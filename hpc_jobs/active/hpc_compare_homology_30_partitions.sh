#!/usr/bin/env bash
# Compare the independently generated framework and Daniel 30% UniRef50 partitions.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=8G
#$ -l tscratch=16G
#$ -l scratch0free=40G
#$ -l h_rt=8:0:0
#$ -pe smp 2
#$ -N hom30_cmp
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }

LEFT_ASSIGNMENTS="${LEFT_ASSIGNMENTS:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/mmseqs_cluster_assignments/uniref50_sensitivity_4/uniref50_2026_02/identity_30/contract_e3e44179dd219d65/cluster_assignments.tsv.gz}"
RIGHT_ASSIGNMENTS="${RIGHT_ASSIGNMENTS:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/supervisor_daniel_buchan/mmseqs_cluster_assignments/uniref50_sensitivity_4/identity_30/raw/cluster_assignments.tsv.gz}"
RESULTS_ROOT="${RESULTS_ROOT:-/SAN/bioinf/bmpfp/diagnostics/homology_cluster_comparison/2026_02/identity_30/framework_vs_daniel}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date -u +%Y%m%dT%H%M%SZ)"
WORK="/scratch0/homology_cluster_comparison_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/framework"
SCRATCH_OUTPUT="$WORK/output"
SORT_SCRATCH="$WORK/sort"
FINAL_OUTPUT="$RESULTS_ROOT/$RUN_TAG"
FAILED_OUTPUT="${FINAL_OUTPUT}.failed"
PUBLISHED=0

publish() {
  local status="$1" destination="$FINAL_OUTPUT" staging="${FINAL_OUTPUT}.staging-${JOB_TOKEN}"
  [[ "$PUBLISHED" == 0 ]] || return 0
  if [[ "$status" != 0 ]]; then
    destination="$FAILED_OUTPUT"
    staging="${FAILED_OUTPUT}.staging-${JOB_TOKEN}"
  fi
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging"
  [[ ! -d "$SCRATCH_OUTPUT" ]] || cp -a "$SCRATCH_OUTPUT/." "$staging/"
  if [[ "$status" != 0 ]]; then
    rm -f "$staging/RUN_COMPLETE.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" > "$staging/WORKFLOW_FAILED.json"
  fi
  mv "$staging" "$destination"
  PUBLISHED=1
  echo "Published comparison: $destination"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish "$status" || publish_status=$?
  if [[ -d "$WORK" && ! -L "$WORK" && "$WORK" == /scratch0/homology_cluster_comparison_* ]]; then
    cd "$HOME"
    rm -rf -- "$WORK"
  else
    echo "Refusing unsafe scratch cleanup: $WORK" >&2
    [[ "$status" != 0 ]] || status=1
  fi
  if [[ "$status" == 0 && "$publish_status" != 0 ]]; then status=$publish_status; fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ -f "$LEFT_ASSIGNMENTS" ]] || die "Framework assignment file is missing: $LEFT_ASSIGNMENTS"
[[ -f "$RIGHT_ASSIGNMENTS" ]] || die "Daniel assignment file is missing: $RIGHT_ASSIGNMENTS"
[[ "$RESULTS_ROOT" == /SAN/* ]] || die "RESULTS_ROOT must be on SAN"
[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK" "$SCRATCH_OUTPUT" "$SORT_SCRATCH" "$RESULTS_ROOT"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a framework checkout"
  FRAMEWORK_COMMIT="$(git -C "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"

echo "Host             : $(hostname -f 2>/dev/null || hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "Framework input  : $LEFT_ASSIGNMENTS"
echo "Daniel input     : $RIGHT_ASSIGNMENTS"
echo "Final output     : $FINAL_OUTPUT"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git -C "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"

python3 "$FRAMEWORK_DIR/scripts/diagnostics/compare_homology_cluster_assignments.py" \
  --left "$LEFT_ASSIGNMENTS" \
  --right "$RIGHT_ASSIGNMENTS" \
  --left-label framework-mmseqs-18-8cc5c \
  --right-label daniel-supervisor-generated \
  --output-dir "$SCRATCH_OUTPUT" \
  --scratch-dir "$SORT_SCRATCH" \
  --sort-parallel "${NSLOTS:-2}" \
  --sort-memory 12G

[[ -f "$SCRATCH_OUTPUT/RUN_COMPLETE.json" ]] || die "Comparison completion marker is missing"
publish 0
echo "Finished         : $(date -Is)"
echo "Read first       : $FINAL_OUTPUT/summary.md"
