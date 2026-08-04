#!/usr/bin/env bash
# Resolve the accepted 30% homology reuse plan to exact cache archive members.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=12G
#$ -l tscratch=8G
#$ -l scratch0free=12G
#$ -l h_rt=8:0:0
#$ -pe smp 1
#$ -N hom30_src
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

COARSE_LEDGER_ROOT="${COARSE_LEDGER_ROOT:-/SAN/bioinf/bmpfp/embedding_states/homology/2026_02/identity_30/cluster-count-random/reuse_ledger/7128474_20260801_045130}"
COARSE_PLAN_DIR="${COARSE_PLAN_DIR:-$COARSE_LEDGER_ROOT/plan}"
RESULTS_ROOT="${RESULTS_ROOT:-/SAN/bioinf/bmpfp/embedding_states/homology/2026_02/identity_30/cluster-count-random/reuse_ledger/source_resolved}"

CONTEMPORARY_ARCHIVE="${CONTEMPORARY_ARCHIVE:-/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/variants/text-cutoff-2025-03-08__ppi-paper-faithful/finalized_pfp_cache/contemporary_embedding_cache.tar.gz}"
CONTEMPORARY_ARCHIVE_SHA256="${CONTEMPORARY_ARCHIVE_SHA256:-8c579d492a9e9ee93a3539f722e479da9f917c6aada2f0238893b263066aa70e}"
CAFA3_ARCHIVE="${CAFA3_ARCHIVE:-/SAN/bioinf/bmpfp/diagnostics/cafa3_embedding_hydration_comparison/7125100_20260730T154052Z/artifacts/cafa3_reproduction_hydrated_cache.tar.gz}"
CAFA3_ARCHIVE_SHA256="${CAFA3_ARCHIVE_SHA256:-c6cdcafd00b0cb871a50beb8cd649ce5c13d5882fcde9c3663a19a86905b9e87}"
CONTEMPORARY_EMBEDDED_BENCHMARK="${CONTEMPORARY_EMBEDDED_BENCHMARK:-contemporary_hydrated_population}"
CAFA3_EMBEDDED_BENCHMARK="${CAFA3_EMBEDDED_BENCHMARK:-cafa3_hydrated_population}"
CONTEMPORARY_TEXT_REUSE_POLICY="${CONTEMPORARY_TEXT_REUSE_POLICY:-source-current}"
CAFA3_TEXT_REUSE_POLICY="${CAFA3_TEXT_REUSE_POLICY:-source-current}"

JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date -u +%Y%m%dT%H%M%SZ)"
WORK="/scratch0/homology_embedding_source_ledger_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
SCRATCH_OUTPUT="$WORK/source_resolved"
WORKFLOW_LOG="$WORK/workflow.log"
FINAL_OUTPUT="$RESULTS_ROOT/$RUN_TAG"
FAILED_OUTPUT="${FINAL_OUTPUT}.failed"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
WORK_OWNED=0
PUBLISHED=0

publish() {
  local status="$1" destination="$FINAL_OUTPUT"
  local staging="${FINAL_OUTPUT}.staging-${JOB_TOKEN}" copy_status=0
  [[ "$PUBLISHED" == 0 ]] || return 0
  if [[ "$status" != 0 ]]; then
    destination="$FAILED_OUTPUT"
    staging="${FAILED_OUTPUT}.staging-${JOB_TOKEN}"
  fi
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging/logs" || return 1
  [[ ! -d "$SCRATCH_OUTPUT" ]] || cp -a "$SCRATCH_OUTPUT/." "$staging/" || copy_status=$?
  [[ ! -f "$WORKFLOW_LOG" ]] || cp -p "$WORKFLOW_LOG" "$staging/logs/workflow.log" || copy_status=$?
  if [[ "$status" == 0 ]]; then
    [[ -f "$staging/RUN_COMPLETE.json" ]] || copy_status=1
    [[ -f "$staging/output_manifest.json" ]] || copy_status=1
    [[ -f "$staging/resolved_embedding_pairs.tsv.gz" ]] || copy_status=1
  else
    rm -f "$staging/RUN_COMPLETE.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" \
      > "$staging/WORKFLOW_FAILED.json" || copy_status=$?
  fi
  if [[ "$copy_status" == 0 ]]; then mv "$staging" "$destination" || copy_status=$?; fi
  if [[ "$copy_status" == 0 ]]; then
    PUBLISHED=1
    echo "Published source-resolved ledger: $destination"
  elif [[ -d "$staging" && ! -L "$staging" ]]; then
    rm -rf -- "$staging"
  fi
  return "$copy_status"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish "$status" || publish_status=$?
  if [[ "$WORK_OWNED" == 1 && "$WORK" == /scratch0/homology_embedding_source_ledger_* && ! -L "$WORK" ]]; then
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

[[ -d "$COARSE_PLAN_DIR" ]] || die "Coarse plan is missing: $COARSE_PLAN_DIR"
[[ -f "$CONTEMPORARY_ARCHIVE" ]] || die "Contemporary archive is missing"
[[ -f "$CAFA3_ARCHIVE" ]] || die "CAFA3 archive is missing"
[[ "$CONTEMPORARY_EMBEDDED_BENCHMARK" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || \
  die "Unsafe contemporary embedded-benchmark name"
[[ "$CAFA3_EMBEDDED_BENCHMARK" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || \
  die "Unsafe CAFA3 embedded-benchmark name"
[[ "$CONTEMPORARY_TEXT_REUSE_POLICY" =~ ^(never|same-role|source-current)$ ]] || \
  die "Invalid contemporary text reuse policy"
[[ "$CAFA3_TEXT_REUSE_POLICY" =~ ^(never|same-role|source-current)$ ]] || \
  die "Invalid CAFA3 text reuse policy"
[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
[[ "$RESULTS_ROOT" == /SAN/* ]] || die "RESULTS_ROOT must be on SAN"
mkdir -p "$WORK" "$RESULTS_ROOT"
WORK_OWNED=1

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a framework checkout"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"

echo "Host                 : $(hostname -f 2>/dev/null || hostname)"
echo "Job ID               : ${JOB_ID:-manual}"
echo "Framework commit     : $FRAMEWORK_COMMIT"
echo "Coarse plan          : $COARSE_PLAN_DIR"
echo "Contemporary source  : $CONTEMPORARY_ARCHIVE"
echo "Contemporary label   : $CONTEMPORARY_EMBEDDED_BENCHMARK"
echo "Contemporary text    : $CONTEMPORARY_TEXT_REUSE_POLICY"
echo "CAFA3 source         : $CAFA3_ARCHIVE"
echo "CAFA3 label          : $CAFA3_EMBEDDED_BENCHMARK"
echo "CAFA3 text           : $CAFA3_TEXT_REUSE_POLICY"
echo "Final output         : $FINAL_OUTPUT"
echo "Started              : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind /SAN/bioinf/bmpfp
add_mmfp_singularity_bind "$WORK"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

set +e
"$PYTHON_BIN" scripts/embeddings/resolve_embedding_reuse_sources.py \
  --coarse-plan-dir "$COARSE_PLAN_DIR" \
  --cache-source "contemporary_paper_faithful=${CONTEMPORARY_EMBEDDED_BENCHMARK}=${CONTEMPORARY_ARCHIVE}=${CONTEMPORARY_ARCHIVE_SHA256}" \
  --cache-source "cafa3_regenerated_hydrated=${CAFA3_EMBEDDED_BENCHMARK}=${CAFA3_ARCHIVE}=${CAFA3_ARCHIVE_SHA256}" \
  --source-text-policy "contemporary_paper_faithful=${CONTEMPORARY_TEXT_REUSE_POLICY}" \
  --source-text-policy "cafa3_regenerated_hydrated=${CAFA3_TEXT_REUSE_POLICY}" \
  --output-dir "$SCRATCH_OUTPUT" \
  2>&1 | tee "$WORKFLOW_LOG"
status=${PIPESTATUS[0]}
set -e
[[ "$status" == 0 ]] || exit "$status"

publish 0
echo "Finished             : $(date -Is)"
echo "Read first           : $FINAL_OUTPUT/summary.md"
