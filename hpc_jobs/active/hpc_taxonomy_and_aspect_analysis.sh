#!/usr/bin/env bash
# Analyse taxonomy and exact GO-aspect membership for accepted evidence-12 benchmarks.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=8G
#$ -l tscratch=8G
#$ -l scratch0free=10G
#$ -l h_rt=4:0:0
#$ -pe smp 1
#$ -N pfp_taxasp
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

E12="${E12_ROOT:-/SAN/bioinf/bmpfp/reruns/fixed-evidence-codes/evidence12-20260824T072711Z}"
RND="${RANDOM_ROOT:-/SAN/bioinf/bmpfp/reruns/fixed-evidence-codes/random-controls-evidence12-20260907T052655Z}"
HBASE="$E12/benchmarks/homology/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_bb47beee0e57/run_evidence12-20260824T072711Z/job_7228712"
GLOBAL_ROOT="$E12/benchmarks/global_nk/benchmark"
NKLK_ROOT="$E12/benchmarks/nk_lk/benchmark"
HOMOLOGY_ROOT="$HBASE/task_1_identity_30/benchmark/source_sprot-and-trembl/framework_bb47beee0e57/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_30/cluster-count-random/annotated-only/seed_0/min_count_50"
RANDOM_BENCHMARK="$RND/random_h_homology30"
TAXONOMY_DIR="${TAXONOMY_DIR:-/SAN/bioinf/bmpfp/frozen_inputs/taxonomy/ncbi/2026-07-01}"
CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
PYTHON_BIN="${PYTHON_BIN:-/share/apps/miniforge3_mamba/bin/python}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/pfp_taxonomy_aspects_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/framework"
SCRATCH_OUTPUT="$WORK/output"
FINAL_OUTPUT="$CAMPAIGN_ROOT/taxonomy_and_aspects"
FAILED_OUTPUT="$CAMPAIGN_ROOT/taxonomy_and_aspects.failed-${JOB_TOKEN}"
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
  echo "Published: $destination"
}

cleanup() {
  local status=$? publish_status=0
  trap - EXIT
  set +e
  publish "$status" || publish_status=$?
  if [[ -d "$WORK" && ! -L "$WORK" && "$WORK" == /scratch0/pfp_taxonomy_aspects_* ]]; then
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
trap 'exit 130' INT TERM

[[ -n "$CAMPAIGN_ROOT" && "$CAMPAIGN_ROOT" == /SAN/bioinf/bmpfp/reruns/fixed-evidence-codes/* ]] || die "CAMPAIGN_ROOT must be an explicit fixed-evidence-codes SAN path"
[[ -x "$PYTHON_BIN" ]] || die "Python is unavailable: $PYTHON_BIN"
[[ -d "$TAXONOMY_DIR" ]] || die "NCBI taxonomy directory is unavailable: $TAXONOMY_DIR"
[[ ! -e "$FINAL_OUTPUT" && ! -e "$FAILED_OUTPUT" ]] || die "Output already exists"
mkdir -p "$CAMPAIGN_ROOT"
available_kb="$(df -Pk "$CAMPAIGN_ROOT" | awk 'NR==2 {print $4}')"
[[ "$available_kb" =~ ^[0-9]+$ && "$available_kb" -ge 5242880 ]] || die "SAN has less than 5 GiB free"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from the framework checkout or pass FRAMEWORK_COMMIT"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"
mkdir -p "$WORK"
export MPLCONFIGDIR="$WORK/matplotlib"
mkdir -p "$MPLCONFIGDIR"

echo "Host             : $(hostname -f 2>/dev/null || hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "Campaign root    : $CAMPAIGN_ROOT"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"

"$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/diagnostics/analyse_taxonomy_and_aspects.py" \
  --taxonomy-dir "$TAXONOMY_DIR" \
  --global-root "$GLOBAL_ROOT" \
  --nklk-root "$NKLK_ROOT" \
  --homology-root "$HOMOLOGY_ROOT" \
  --random-root "$RANDOM_BENCHMARK" \
  --skip-plots \
  --output-dir "$SCRATCH_OUTPUT"

[[ -f "$SCRATCH_OUTPUT/RUN_COMPLETE.json" ]] || die "Analysis completion marker is missing"
publish 0
echo "Finished         : $(date -Is)"
