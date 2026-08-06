#!/usr/bin/env bash
# Audit direct evidence-code composition without changing homology benchmarks.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=4G
#$ -l h_rt=4:0:0
#$ -pe smp 1
#$ -N hom_ev_audit
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

BENCHMARK_BASE="${BENCHMARK_BASE:-/SAN/bioinf/bmpfp/benchmarks/homology/2026_02}"
FRAMEWORK_PREFIX="$BENCHMARK_BASE/uniref50_sensitivity_4_daniel_aligned_cached_random_resplit/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_88ae24587086/run_runtime-7128983/job_7128983"
FRAMEWORK_OLD_PREFIX="$BENCHMARK_BASE/uniref50_sensitivity_4_daniel_aligned_12core/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_e08c2dc2733e/run_runtime-7127410/job_7127410"
SUFFIX_BASE="benchmark/source_sprot-and-trembl"
ROOT_30="${ROOT_30:-$FRAMEWORK_PREFIX/task_1_identity_30/$SUFFIX_BASE/framework_88ae24587086/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_30/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_25="${ROOT_25:-$FRAMEWORK_PREFIX/task_2_identity_25/$SUFFIX_BASE/framework_88ae24587086/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_25/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_20="${ROOT_20:-$FRAMEWORK_PREFIX/task_3_identity_20/$SUFFIX_BASE/framework_88ae24587086/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_20/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_15="${ROOT_15:-$FRAMEWORK_OLD_PREFIX/task_4_identity_15/$SUFFIX_BASE/framework_e08c2dc2733e/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_15/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_10="${ROOT_10:-$FRAMEWORK_OLD_PREFIX/task_5_identity_10/$SUFFIX_BASE/framework_e08c2dc2733e/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_10/cluster-count-random/annotated-only/seed_0/min_count_50}"
ROOT_05="${ROOT_05:-$FRAMEWORK_OLD_PREFIX/task_6_identity_5/$SUFFIX_BASE/framework_e08c2dc2733e/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_05/cluster-count-random/annotated-only/seed_0/min_count_50}"
DANIEL_30="${DANIEL_30:-$BENCHMARK_BASE/supervisor_daniel_buchan/uniref50_sensitivity_4/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_ffe3c038bb78/run_runtime-7132993/job_7132993/task_1_identity_30/benchmark/source_sprot-and-trembl/framework_ffe3c038bb78/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_30/cluster-count-random/annotated-only/seed_0/min_count_50}"
RESULTS_ROOT="${RESULTS_ROOT:-/SAN/bioinf/bmpfp/diagnostics/homology_evidence_policy/2026_02}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
PYTHON_BIN="${PYTHON_BIN:-/share/apps/miniforge3_mamba/bin/python}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
JOB_TOKEN="${JOB_ID:-manual_$$}"
RUN_TAG="${JOB_TOKEN}_$(date -u +%Y%m%dT%H%M%SZ)"
WORK="/scratch0/homology_evidence_audit_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/framework"
SCRATCH_OUTPUT="$WORK/output"
FINAL_OUTPUT="$RESULTS_ROOT/$RUN_TAG"
PUBLISHED=0

publish() {
  local status="$1" destination="$FINAL_OUTPUT" staging="${FINAL_OUTPUT}.staging-${JOB_TOKEN}"
  [[ "$PUBLISHED" == 0 ]] || return 0
  if [[ "$status" != 0 ]]; then destination="${FINAL_OUTPUT}.failed"; staging="${destination}.staging-${JOB_TOKEN}"; fi
  [[ ! -e "$destination" && ! -e "$staging" ]] || return 1
  mkdir -p "$staging"
  [[ ! -d "$SCRATCH_OUTPUT" ]] || cp -a "$SCRATCH_OUTPUT/." "$staging/"
  if [[ "$status" != 0 ]]; then
    rm -f "$staging/RUN_COMPLETE.json"
    printf '{"complete":false,"workflow_exit_status":%s}\n' "$status" > "$staging/WORKFLOW_FAILED.json"
  fi
  mv "$staging" "$destination"
  PUBLISHED=1
  echo "Published evidence audit: $destination"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish "$status" || publish_status=$?
  if [[ -d "$WORK" && ! -L "$WORK" && "$WORK" == /scratch0/homology_evidence_audit_* ]]; then
    cd "$HOME" && rm -rf -- "$WORK"
  else
    echo "Refusing unsafe scratch cleanup: $WORK" >&2
    [[ "$status" != 0 ]] || status=1
  fi
  if [[ "$status" == 0 && "$publish_status" != 0 ]]; then status=$publish_status; fi
  exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

for root in "$ROOT_30" "$ROOT_25" "$ROOT_20" "$ROOT_15" "$ROOT_10" "$ROOT_05" "$DANIEL_30"; do
  [[ -f "$root/qualifying_annotations.tsv.gz" ]] || die "Missing qualifying annotations: $root"
done
[[ -x "$PYTHON_BIN" ]] || die "Missing Python: $PYTHON_BIN"
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

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"

"$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/diagnostics/audit_homology_evidence_policy.py" \
  --benchmark "framework-30=$ROOT_30" \
  --benchmark "framework-25=$ROOT_25" \
  --benchmark "framework-20=$ROOT_20" \
  --benchmark "framework-15=$ROOT_15" \
  --benchmark "framework-10=$ROOT_10" \
  --benchmark "framework-05=$ROOT_05" \
  --benchmark "daniel-30=$DANIEL_30" \
  --output-dir "$SCRATCH_OUTPUT"

[[ -f "$SCRATCH_OUTPUT/RUN_COMPLETE.json" ]] || die "Completion marker is missing"
publish 0
echo "Read first       : $FINAL_OUTPUT/summary.md"
