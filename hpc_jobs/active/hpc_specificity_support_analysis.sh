#!/usr/bin/env bash
# Test Xu specificity against positive training support on an accepted PFP run.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=8G
#$ -l tscratch=4G
#$ -l h_rt=2:0:0
#$ -pe smp 1
#$ -N pfp_specsup
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

E12="${E12_ROOT:-/SAN/bioinf/bmpfp/reruns/fixed-evidence-codes/evidence12-20260824T072711Z}"
RUN_ROOT="${RUN_ROOT:-$E12/model_runs/homology/identity_30/full/run}"
PREDICTION_MANIFEST="${PREDICTION_MANIFEST:-$RUN_ROOT/evaluation/prediction_artifacts/prediction_artifact_manifest.json}"
PREPARED_DATA_DIR="${PREPARED_DATA_DIR:-$RUN_ROOT/prepared_data}"
SPECIFICITY_DIR="${SPECIFICITY_DIR:-$E12/analyses/homology/identity_30/diagnostics/analysis/specificity}"
OUTPUT_DIR="${OUTPUT_DIR:-}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/mmfp/bin/python}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/pfp_specificity_support_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/framework"
SCRATCH_OUTPUT="$WORK/output"
FINAL_OUTPUT="$OUTPUT_DIR"
FAILED_OUTPUT="${OUTPUT_DIR}.failed-${JOB_TOKEN}"
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
  if [[ -d "$WORK" && ! -L "$WORK" && "$WORK" == /scratch0/pfp_specificity_support_* ]]; then
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

case "$OUTPUT_DIR" in
  "$HOME"/pfp-runs/specificity-support-*|/SAN/bioinf/bmpfp/reruns/fixed-evidence-codes/specificity-support-*) ;;
  *) die "OUTPUT_DIR must be a new specificity-support path under HOME/pfp-runs or fixed-evidence-codes SAN" ;;
esac
[[ -x "$PYTHON_BIN" ]] || die "Python is unavailable: $PYTHON_BIN"
[[ -f "$PREDICTION_MANIFEST" ]] || die "Prediction manifest is unavailable"
[[ -d "$PREPARED_DATA_DIR" ]] || die "Prepared data are unavailable"
[[ -d "$SPECIFICITY_DIR" ]] || die "Specificity analysis is unavailable"
[[ ! -e "$FINAL_OUTPUT" && ! -e "$FAILED_OUTPUT" ]] || die "Output already exists"
mkdir -p "$(dirname "$OUTPUT_DIR")"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from the framework checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout has uncommitted changes"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"

mkdir -p "$WORK"
export MPLCONFIGDIR="$WORK/matplotlib"
mkdir -p "$MPLCONFIGDIR"

echo "Host             : $(hostname -f 2>/dev/null || hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "Output directory : $OUTPUT_DIR"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || die "Framework checkout differs from requested commit"

"$PYTHON_BIN" -c 'import matplotlib, numpy, scipy'
"$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/diagnostics/analyze_pfp_specificity_support.py" \
  --prediction-manifest "$PREDICTION_MANIFEST" \
  --prepared-data-dir "$PREPARED_DATA_DIR" \
  --specificity-dir "$SPECIFICITY_DIR" \
  --candidate-bins 4 \
  --minimum-cell-terms 10 \
  --output-dir "$SCRATCH_OUTPUT"

[[ -f "$SCRATCH_OUTPUT/RUN_COMPLETE.json" ]] || die "Analysis completion marker is missing"
publish 0
echo "Finished         : $(date -Is)"
