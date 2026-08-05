#!/usr/bin/env bash
# Audit published model inputs and retained clusters across six identity thresholds.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=4G
#$ -l h_rt=8:0:0
#$ -pe smp 1
#$ -N hom_thr_audit
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

BENCHMARK_BASE="${BENCHMARK_BASE:-/SAN/bioinf/bmpfp/benchmarks/homology/2026_02}"
ROOT_30="${ROOT_30:-$BENCHMARK_BASE/uniref50_sensitivity_4_daniel_aligned_cached_random_resplit/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_88ae24587086/run_runtime-7128983/job_7128983/task_1_identity_30/benchmark/source_sprot-and-trembl/framework_88ae24587086/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_30/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_25="${ROOT_25:-$BENCHMARK_BASE/uniref50_sensitivity_4_daniel_aligned_cached_random_resplit/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_88ae24587086/run_runtime-7128983/job_7128983/task_2_identity_25/benchmark/source_sprot-and-trembl/framework_88ae24587086/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_25/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_20="${ROOT_20:-$BENCHMARK_BASE/uniref50_sensitivity_4_daniel_aligned_cached_random_resplit/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_88ae24587086/run_runtime-7128983/job_7128983/task_3_identity_20/benchmark/source_sprot-and-trembl/framework_88ae24587086/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_20/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_15="${ROOT_15:-$BENCHMARK_BASE/uniref50_sensitivity_4_daniel_aligned_12core/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_e08c2dc2733e/run_runtime-7127410/job_7127410/task_4_identity_15/benchmark/source_sprot-and-trembl/framework_e08c2dc2733e/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_15/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_10="${ROOT_10:-$BENCHMARK_BASE/uniref50_sensitivity_4_daniel_aligned_12core/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_e08c2dc2733e/run_runtime-7127410/job_7127410/task_5_identity_10/benchmark/source_sprot-and-trembl/framework_e08c2dc2733e/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_10/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_05="${ROOT_05:-$BENCHMARK_BASE/uniref50_sensitivity_4_daniel_aligned_12core/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_e08c2dc2733e/run_runtime-7127410/job_7127410/task_6_identity_5/benchmark/source_sprot-and-trembl/framework_e08c2dc2733e/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_05/cluster-count-random/annotated-only/seed_0/min_count_50}"
CLUSTER_BASE="${CLUSTER_BASE:-/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/mmseqs_cluster_assignments/uniref50_sensitivity_4/uniref50_2026_02}"
FULL_30="${FULL_30:-$CLUSTER_BASE/identity_30/contract_e3e44179dd219d65/cluster_assignments.tsv.gz}"
FULL_25="${FULL_25:-$CLUSTER_BASE/identity_25/contract_a6fd69df69a78cb0/cluster_assignments.tsv.gz}"
FULL_20="${FULL_20:-$CLUSTER_BASE/identity_20/contract_fe783953dbe1c275/cluster_assignments.tsv.gz}"
FULL_15="${FULL_15:-$CLUSTER_BASE/identity_15/contract_fb6fcfd60847780e/cluster_assignments.tsv.gz}"
FULL_10="${FULL_10:-$CLUSTER_BASE/identity_10/contract_9206700890e1ccce/cluster_assignments.tsv.gz}"
FULL_05="${FULL_05:-$CLUSTER_BASE/identity_05/contract_abd0c5c4beaae966/cluster_assignments.tsv.gz}"
RESULTS_ROOT="${RESULTS_ROOT:-/SAN/bioinf/bmpfp/diagnostics/homology_threshold_progression/2026_02/framework_stream/benchmark_inputs}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
PYTHON_BIN="${PYTHON_BIN:-/share/apps/miniforge3_mamba/bin/python}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date -u +%Y%m%dT%H%M%SZ)"
WORK="/scratch0/homology_threshold_audit_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/framework"
SCRATCH_OUTPUT="$WORK/output"
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
  echo "Published threshold audit: $destination"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish "$status" || publish_status=$?
  if [[ -d "$WORK" && ! -L "$WORK" && "$WORK" == /scratch0/homology_threshold_audit_* ]]; then
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

for root in "$ROOT_30" "$ROOT_25" "$ROOT_20" "$ROOT_15" "$ROOT_10" "$ROOT_05"; do
  [[ -f "$root/RUN_COMPLETE.json" ]] || die "Completed benchmark is unavailable: $root"
done
for assignment in "$FULL_30" "$FULL_25" "$FULL_20" "$FULL_15" "$FULL_10" "$FULL_05"; do
  [[ -f "$assignment" ]] || die "Full cluster assignment is unavailable: $assignment"
done
[[ -x "$PYTHON_BIN" ]] || die "Shared Python executable is unavailable: $PYTHON_BIN"
[[ "$RESULTS_ROOT" == /SAN/* ]] || die "RESULTS_ROOT must be on SAN"
[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK" "$RESULTS_ROOT"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a framework checkout"
  FRAMEWORK_COMMIT="$(git -C "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"

echo "Host             : $(hostname -f 2>/dev/null || hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "Final output     : $FINAL_OUTPUT"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"

"$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/diagnostics/audit_homology_threshold_progression.py" \
  --benchmark "30=$ROOT_30" \
  --benchmark "25=$ROOT_25" \
  --benchmark "20=$ROOT_20" \
  --benchmark "15=$ROOT_15" \
  --benchmark "10=$ROOT_10" \
  --benchmark "05=$ROOT_05" \
  --full-assignment "30=$FULL_30" \
  --full-assignment "25=$FULL_25" \
  --full-assignment "20=$FULL_20" \
  --full-assignment "15=$FULL_15" \
  --full-assignment "10=$FULL_10" \
  --full-assignment "05=$FULL_05" \
  --output-dir "$SCRATCH_OUTPUT"

[[ -f "$SCRATCH_OUTPUT/RUN_COMPLETE.json" ]] || die "Audit completion marker is missing"
publish 0
echo "Finished         : $(date -Is)"
echo "Read first       : $FINAL_OUTPUT/summary.md"
