#!/usr/bin/env bash
# Compare Daniel Buchan's and the framework's 25%, 20% and 15% UniRef50 partitions.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=8G
#$ -l tscratch=16G
#$ -l scratch0free=40G
#$ -l h_rt=12:0:0
#$ -pe smp 2
#$ -t 1-3
#$ -tc 3
#$ -N hom_df_cmp
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

TASK_ID="${SGE_TASK_ID:-${1:-}}"
[[ "$TASK_ID" =~ ^[1-3]$ ]] || die "SGE_TASK_ID must be 1-3"
IDENTITIES=(25 20 15)
CONTRACTS=(a6fd69df69a78cb0 fe783953dbe1c275 fb6fcfd60847780e)
INDEX=$((TASK_ID - 1))
IDENTITY="${IDENTITIES[$INDEX]}"
CONTRACT="${CONTRACTS[$INDEX]}"

FRAMEWORK_BASE="${FRAMEWORK_BASE:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/mmseqs_cluster_assignments/uniref50_sensitivity_4/uniref50_2026_02}"
DANIEL_BASE="${DANIEL_BASE:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/supervisor_daniel_buchan/mmseqs_cluster_assignments/uniref50_sensitivity_4}"
LEFT_ASSIGNMENTS="${LEFT_ASSIGNMENTS:-$FRAMEWORK_BASE/identity_${IDENTITY}/contract_${CONTRACT}/cluster_assignments.tsv.gz}"
RIGHT_ASSIGNMENTS="${RIGHT_ASSIGNMENTS:-$DANIEL_BASE/identity_${IDENTITY}/raw/cluster_assignments.tsv.gz}"
RESULTS_ROOT="${RESULTS_ROOT:-/SAN/bioinf/bmpfp/diagnostics/homology_cluster_comparison/2026_02/identity_${IDENTITY}/framework_vs_daniel}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
PYTHON_BIN="${PYTHON_BIN:-/share/apps/miniforge3_mamba/bin/python}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
JOB_TOKEN="${JOB_ID:-manual_$$}_${TASK_ID}"
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
[[ -x "$PYTHON_BIN" ]] || die "Shared Python executable is unavailable: $PYTHON_BIN"
[[ "$RESULTS_ROOT" == /SAN/* ]] || die "RESULTS_ROOT must be on SAN"
[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK" "$SCRATCH_OUTPUT" "$SORT_SCRATCH" "$RESULTS_ROOT"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a framework checkout"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"

echo "Host             : $(hostname -f 2>/dev/null || hostname)"
echo "Job/task ID      : ${JOB_ID:-manual}/${TASK_ID}"
echo "Identity         : ${IDENTITY}%"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "Framework input  : $LEFT_ASSIGNMENTS"
echo "Daniel input     : $RIGHT_ASSIGNMENTS"
echo "Final output     : $FINAL_OUTPUT"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"

"$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/diagnostics/compare_homology_cluster_assignments.py" \
  --left "$LEFT_ASSIGNMENTS" \
  --right "$RIGHT_ASSIGNMENTS" \
  --left-label "framework-identity-${IDENTITY}" \
  --right-label "daniel-identity-${IDENTITY}" \
  --title "Framework versus Daniel UniRef50 ${IDENTITY}% full-partition comparison" \
  --interpretation-boundary "Both partitions cover the same frozen UniRef50 identifier universe and nominal identity threshold. Agreement measures implementation-level partition stability; it does not establish exhaustive pairwise sequence-identity separation or identify the cause of any disagreement." \
  --output-dir "$SCRATCH_OUTPUT" \
  --scratch-dir "$SORT_SCRATCH" \
  --sort-parallel "${NSLOTS:-2}" \
  --sort-memory 12G

[[ -f "$SCRATCH_OUTPUT/RUN_COMPLETE.json" ]] || die "Comparison completion marker is missing"
publish 0
echo "Finished         : $(date -Is)"
echo "Read first       : $FINAL_OUTPUT/summary.md"
