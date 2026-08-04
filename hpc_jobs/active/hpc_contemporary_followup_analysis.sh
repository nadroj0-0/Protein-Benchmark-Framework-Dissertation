#!/usr/bin/env bash
# Run specificity or validation-fitted calibration for the corrected contemporary full model.

#$ -S /bin/bash
#$ -l tmem=32G
#$ -l scratch0free=20G
#$ -l tscratch=20G
#$ -l h_rt=24:0:0
#$ -pe smp 2
#$ -j y
#$ -N ct25_followup
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  cat <<'EOF'
Usage:
  qsub hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
    --analysis specificity --output-dir /SAN/.../specificity

  qsub hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
    --analysis specificity --source-run /SAN/.../model_run \
    --obo /SAN/.../go.obo --source-label cafa3-regenerated-full \
    --output-dir /SAN/.../specificity

  qsub -hold_jid CAPTURE_JOB hpc_jobs/active/hpc_contemporary_followup_analysis.sh \
    --analysis calibration --capture-pair-dir /SAN/.../capture_pair \
    --output-dir /SAN/.../calibration

Both modes are CPU-only. Specificity uses the accepted corrected full-model
test arrays. Calibration fits on a newly captured validation split and applies
the fitted model once to its paired test split.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

ANALYSIS=""
CAPTURE_PAIR_DIR=""
OUTPUT_DIR=""
SOURCE_RUN_OVERRIDE="${SOURCE_RUN:-}"
OBO_FILE_OVERRIDE="${OBO_FILE:-}"
SOURCE_LABEL="${SOURCE_LABEL:-contemporary-corrected-text-paper-faithful-ppi-full}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --analysis) require_value "$@"; ANALYSIS="$2"; shift 2 ;;
    --capture-pair-dir) require_value "$@"; CAPTURE_PAIR_DIR="$2"; shift 2 ;;
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    --source-run) require_value "$@"; SOURCE_RUN_OVERRIDE="$2"; shift 2 ;;
    --obo) require_value "$@"; OBO_FILE_OVERRIDE="$2"; shift 2 ;;
    --source-label) require_value "$@"; SOURCE_LABEL="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "$ANALYSIS" == "specificity" || "$ANALYSIS" == "calibration" ]] || \
  die "--analysis must be specificity or calibration"
[[ "$OUTPUT_DIR" == /SAN/* ]] || die "--output-dir must be an absolute SAN path"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"
[[ "$SOURCE_LABEL" =~ ^[A-Za-z0-9._-]+$ ]] || \
  die "--source-label must use only letters, digits, dot, underscore or hyphen"
if [[ "$ANALYSIS" == "calibration" ]]; then
  [[ "$CAPTURE_PAIR_DIR" == /SAN/* ]] || \
    die "Calibration requires --capture-pair-dir on SAN"
fi

SOURCE_RUN="${SOURCE_RUN_OVERRIDE:-/SAN/bioinf/bmpfp/model_runs/contemporary/2025_01_to_2026_02_supervisor/variants/text-cutoff-2025-03-08__ppi-paper-faithful/full/7118745_20260728_164527}"
OBO_FILE="${OBO_FILE_OVERRIDE:-/SAN/bioinf/bmpfp/frozen_inputs/ontology/2025-02-06/go-basic.obo}"
[[ "$SOURCE_RUN" == /SAN/* ]] || die "--source-run must be an absolute SAN path"
[[ "$OBO_FILE" == /SAN/* ]] || die "--obo must be an absolute SAN path"
SOURCE_TEST_MANIFEST="$SOURCE_RUN/evaluation/prediction_artifacts/prediction_artifact_manifest.json"
[[ -f "$SOURCE_TEST_MANIFEST" ]] || die "Accepted test prediction manifest is missing"
[[ -f "$OBO_FILE" ]] || die "Frozen ontology is missing: $OBO_FILE"
PROVENANCE_SOURCE_RUN="$SOURCE_RUN"
PROVENANCE_TEST_MANIFEST="$SOURCE_TEST_MANIFEST"
PROVENANCE_VALIDATION_MANIFEST=""

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/ct25_followup_${ANALYSIS}_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
ANALYSIS_OUTPUT="$WORK/analysis"
LOG_FILE="$WORK/${ANALYSIS}.log"
PUBLISH_STAGE="${OUTPUT_DIR}.staging-${JOB_TOKEN}"
PUBLISH_LOCK="${OUTPUT_DIR}.publish-lock"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
LOCK_HELD=0

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  if [[ "$LOCK_HELD" == "1" && -d "$PUBLISH_LOCK" && ! -L "$PUBLISH_LOCK" ]]; then
    rmdir -- "$PUBLISH_LOCK"
  fi
  [[ ! -d "$PUBLISH_STAGE" || -L "$PUBLISH_STAGE" ]] || rm -rf -- "$PUBLISH_STAGE"
  [[ ! -d "$WORK" || -L "$WORK" || "$WORK" != /scratch0/ct25_followup_* ]] || \
    rm -rf -- "$WORK"
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
[[ ! -e "$PUBLISH_STAGE" ]] || die "Publication stage already exists"
[[ ! -e "$PUBLISH_LOCK" ]] || die "Publication lock already exists"
mkdir -p "$WORK" "$(dirname "$OUTPUT_DIR")"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a clean framework checkout"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || \
    die "Submission checkout has uncommitted changes"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || \
  die "FRAMEWORK_COMMIT must be a full commit"

if [[ "$ANALYSIS" == "calibration" ]]; then
  [[ -f "$CAPTURE_PAIR_DIR/WORKFLOW_COMPLETE.json" ]] || \
    die "Paired capture is incomplete: $CAPTURE_PAIR_DIR"
  VALIDATION_MANIFEST="$CAPTURE_PAIR_DIR/valid/evaluation/prediction_artifacts/prediction_artifact_manifest.json"
  TEST_MANIFEST="$CAPTURE_PAIR_DIR/test/evaluation/prediction_artifacts/prediction_artifact_manifest.json"
  [[ -f "$VALIDATION_MANIFEST" && -f "$TEST_MANIFEST" ]] || \
    die "Paired validation/test prediction manifests are missing"
  PROVENANCE_SOURCE_RUN="$CAPTURE_PAIR_DIR"
  PROVENANCE_TEST_MANIFEST="$TEST_MANIFEST"
  PROVENANCE_VALIDATION_MANIFEST="$VALIDATION_MANIFEST"
fi

echo "Host             : $(hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Analysis         : $ANALYSIS"
echo "Source label     : $SOURCE_LABEL"
echo "Source run       : $SOURCE_RUN"
echo "Ontology         : $OBO_FILE"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "SAN output       : $OUTPUT_DIR"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || \
  die "Framework checkout differs from requested commit"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
add_mmfp_singularity_bind /SAN/bioinf/bmpfp
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

if [[ "$ANALYSIS" == "specificity" ]]; then
  "$PYTHON_BIN" scripts/diagnostics/evaluate_pfp_information_content.py \
    --prediction-manifest "$SOURCE_TEST_MANIFEST" \
    --obo "$OBO_FILE" \
    --specificity-measure all_separate \
    --positive-bins 4 \
    --bootstrap-replicates 2000 \
    --bootstrap-seed 42 \
    --output-dir "$ANALYSIS_OUTPUT" \
    >"$LOG_FILE" 2>&1
else
  "$PYTHON_BIN" scripts/diagnostics/calibrate_pfp_predictions.py \
    --validation-prediction-manifest "$VALIDATION_MANIFEST" \
    --test-prediction-manifest "$TEST_MANIFEST" \
    --obo "$OBO_FILE" \
    --positive-ia-bins 4 \
    --reliability-bins 10 \
    --output-dir "$ANALYSIS_OUTPUT" \
    >"$LOG_FILE" 2>&1
fi

echo "==> Publishing $ANALYSIS analysis atomically"
mkdir -p "$PUBLISH_STAGE/logs"
cp -a "$ANALYSIS_OUTPUT" "$PUBLISH_STAGE/analysis"
cp -p "$LOG_FILE" "$PUBLISH_STAGE/logs/"
"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py write \
  --root "$PUBLISH_STAGE" --include-nested-control-files
MANIFEST_SHA256="$(
  "$PYTHON_BIN" -c 'import hashlib,pathlib,sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' \
    "$PUBLISH_STAGE/output_manifest.json"
)"
"$PYTHON_BIN" - "$PUBLISH_STAGE/WORKFLOW_COMPLETE.json" "$ANALYSIS" \
  "$SOURCE_LABEL" "$PROVENANCE_SOURCE_RUN" "$PROVENANCE_TEST_MANIFEST" \
  "$OBO_FILE" "$MANIFEST_SHA256" "$PROVENANCE_VALIDATION_MANIFEST" <<'PY'
import hashlib
import json
import pathlib
import sys

def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

payload = {
    "complete": True,
    "analysis_kind": sys.argv[2],
    "source_label": sys.argv[3],
    "source_run": str(pathlib.Path(sys.argv[4]).resolve()),
    "prediction_manifest": str(pathlib.Path(sys.argv[5]).resolve()),
    "prediction_manifest_sha256": sha256(sys.argv[5]),
    "obo_file": str(pathlib.Path(sys.argv[6]).resolve()),
    "obo_sha256": sha256(sys.argv[6]),
    "mode": "full",
    "manifest": "output_manifest.json",
    "manifest_sha256": sys.argv[7],
}
if sys.argv[8]:
    payload["validation_prediction_manifest"] = str(
        pathlib.Path(sys.argv[8]).resolve()
    )
    payload["validation_prediction_manifest_sha256"] = sha256(sys.argv[8])
pathlib.Path(sys.argv[1]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY
"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py verify \
  --root "$PUBLISH_STAGE" --include-nested-control-files
mkdir "$PUBLISH_LOCK"
LOCK_HELD=1
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory appeared during publication"
mv -T "$PUBLISH_STAGE" "$OUTPUT_DIR"
rmdir "$PUBLISH_LOCK"
LOCK_HELD=0
echo "Published $ANALYSIS analysis: $OUTPUT_DIR"
