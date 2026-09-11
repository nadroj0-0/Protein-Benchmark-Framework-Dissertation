#!/usr/bin/env bash
# Analyse supervised cluster sparsity and Random-H same-cluster exposure.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=8G
#$ -l tscratch=8G
#$ -l scratch0free=10G
#$ -l h_rt=8:0:0
#$ -pe smp 1
#$ -N pfp_hmech
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

E12="${E12_ROOT:-/SAN/bioinf/bmpfp/reruns/fixed-evidence-codes/evidence12-20260824T072711Z}"
RND="${RANDOM_ROOT:-/SAN/bioinf/bmpfp/reruns/fixed-evidence-codes/random-controls-evidence12-20260907T052655Z}"
HBASE="$E12/benchmarks/homology/runtime_array/source_sprot-and-trembl/uniref50_sensitivity_4/framework_bb47beee0e57/run_evidence12-20260824T072711Z/job_7228712"
SUFFIX="cluster-count-random/annotated-only/seed_0/min_count_50"
T30="$HBASE/task_1_identity_30/benchmark/source_sprot-and-trembl/framework_bb47beee0e57/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_30/$SUFFIX"
T25="$HBASE/task_2_identity_25/benchmark/source_sprot-and-trembl/framework_bb47beee0e57/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_25/$SUFFIX"
T20="$HBASE/task_3_identity_20/benchmark/source_sprot-and-trembl/framework_bb47beee0e57/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_20/$SUFFIX"
T15="$HBASE/task_4_identity_15/benchmark/source_sprot-and-trembl/framework_bb47beee0e57/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_15/$SUFFIX"
T10="$HBASE/task_5_identity_10/benchmark/source_sprot-and-trembl/framework_bb47beee0e57/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_10/$SUFFIX"
T05="$HBASE/task_6_identity_5/benchmark/source_sprot-and-trembl/framework_bb47beee0e57/uniref50_sensitivity_4/mmseqs_daniel-aligned-defaults/identity_05/$SUFFIX"
RANDOM_BENCHMARK="$RND/random_h_homology30"
RANDOM_PREDICTIONS="$RND/model_runs/random_h/full_retry/run/evaluation/prediction_artifacts"
STANDARD_PREDICTIONS="$E12/model_runs/homology/identity_30/full/run/evaluation/prediction_artifacts"
GO_OBO="${GO_OBO:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2026-06-19/go-basic.obo}"
CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
PYTHON_BIN="${PYTHON_BIN:-/share/apps/miniforge3_mamba/bin/python}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/pfp_homology_mechanism_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/framework"
SCRATCH_OUTPUT="$WORK/output"
FINAL_OUTPUT="$CAMPAIGN_ROOT/homology_plateau_mechanism"
FAILED_OUTPUT="$CAMPAIGN_ROOT/homology_plateau_mechanism.failed-${JOB_TOKEN}"
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
  if [[ -d "$WORK" && ! -L "$WORK" && "$WORK" == /scratch0/pfp_homology_mechanism_* ]]; then
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

"$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/diagnostics/analyse_homology_plateau_mechanism.py" \
  --threshold-dir "30=$T30" \
  --threshold-dir "25=$T25" \
  --threshold-dir "20=$T20" \
  --threshold-dir "15=$T15" \
  --threshold-dir "10=$T10" \
  --threshold-dir "5=$T05" \
  --random-root "$RANDOM_BENCHMARK" \
  --random-predictions "$RANDOM_PREDICTIONS" \
  --standard-predictions "$STANDARD_PREDICTIONS" \
  --go-obo "$GO_OBO" \
  --bootstrap-replicates 2000 \
  --output-dir "$SCRATCH_OUTPUT"

[[ -f "$SCRATCH_OUTPUT/RUN_COMPLETE.json" ]] || die "Analysis completion marker is missing"
publish 0
echo "Finished         : $(date -Is)"
