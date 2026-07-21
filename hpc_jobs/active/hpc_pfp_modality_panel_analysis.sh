#!/usr/bin/env bash
# Compare a completed five-condition PFP modality panel and its prediction artifacts.

#$ -l tmem=32G
#$ -l scratch0free=20G
#$ -l tscratch=20G
#$ -l h_rt=24:0:0
#$ -pe smp 2
#$ -j y
#$ -N pfp_modality_analysis
#$ -V
#$ -notify

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_pfp_modality_panel_analysis.sh \
  --obo-file FILE --output-dir DIR \
  --run MODE=TRAIN_EVAL_RUN_OR_RESULTS_ROOT \
  --prediction-run MODE=CAPTURE_RUN_OR_RESULTS_ROOT [repeat each five times] \
  [--allow-framework-commit-drift]

Canonical --run inputs must be completed train-eval runs. Prediction inputs
may be those same runs or separately captured eval-only runs, but exact
checkpoint and scientific-contract binding is mandatory. Each path may be a
completed run or a root containing exactly one completed run.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

OBO_FILE=""
OUTPUT_DIR=""
RUN_SPECS=()
PREDICTION_SPECS=()
RUN_COUNT=0
PREDICTION_COUNT=0
ALLOW_FRAMEWORK_COMMIT_DRIFT=0
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --obo-file) require_value "$@"; OBO_FILE="$2"; shift 2 ;;
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    --run) require_value "$@"; RUN_SPECS+=("$2"); RUN_COUNT=$((RUN_COUNT + 1)); shift 2 ;;
    --prediction-run) require_value "$@"; PREDICTION_SPECS+=("$2"); PREDICTION_COUNT=$((PREDICTION_COUNT + 1)); shift 2 ;;
    --allow-framework-commit-drift) ALLOW_FRAMEWORK_COMMIT_DRIFT=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

absolute_path() {
  "$HOST_PYTHON" -c 'import os,sys; p=sys.argv[1]; base=sys.argv[2]; sys.stdout.write(os.path.realpath(p if os.path.isabs(p) else os.path.join(base,p))+"\n")' "$1" "$SUBMISSION_DIR"
}

HOST_PYTHON="$(command -v python3 || command -v python || true)"
[[ -n "$HOST_PYTHON" ]] || die "Python is required to normalize submission paths"
[[ -n "$OBO_FILE" ]] || die "--obo-file is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
OBO_FILE="$(absolute_path "$OBO_FILE")"
OUTPUT_DIR="$(absolute_path "$OUTPUT_DIR")"
[[ -f "$OBO_FILE" ]] || die "GO OBO file is missing: $OBO_FILE"
[[ "$RUN_COUNT" -eq 5 ]] || die "Exactly five --run specifications are required"
[[ "$PREDICTION_COUNT" -eq 5 ]] || die "Exactly five --prediction-run specifications are required"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"

declare -A RUN_INPUTS=()
for specification in "${RUN_SPECS[@]}"; do
  [[ "$specification" == *=* ]] || die "--run must use MODE=PATH"
  mode="${specification%%=*}"
  path="${specification#*=}"
  [[ "$mode" =~ ^(full|sequence-only|sequence-text|sequence-structure|sequence-ppi)$ ]] || \
    die "Unsupported modality mode: $mode"
  [[ -z "${RUN_INPUTS[$mode]+x}" ]] || die "Repeated modality mode: $mode"
  path="$(absolute_path "$path")"
  [[ -d "$path" ]] || die "PFP run path is missing: $path"
  RUN_INPUTS[$mode]="$path"
done
declare -A PREDICTION_INPUTS=()
for specification in "${PREDICTION_SPECS[@]}"; do
  [[ "$specification" == *=* ]] || die "--prediction-run must use MODE=PATH"
  mode="${specification%%=*}"
  path="${specification#*=}"
  [[ "$mode" =~ ^(full|sequence-only|sequence-text|sequence-structure|sequence-ppi)$ ]] || \
    die "Unsupported prediction modality mode: $mode"
  [[ -z "${PREDICTION_INPUTS[$mode]+x}" ]] || die "Repeated prediction modality mode: $mode"
  path="$(absolute_path "$path")"
  [[ -d "$path" ]] || die "PFP prediction run path is missing: $path"
  PREDICTION_INPUTS[$mode]="$path"
done
for mode in full sequence-only sequence-text sequence-structure sequence-ppi; do
  [[ -n "${RUN_INPUTS[$mode]+x}" ]] || die "Missing --run for $mode"
  [[ -n "${PREDICTION_INPUTS[$mode]+x}" ]] || die "Missing --prediction-run for $mode"
done

resolve_completed_run() {
  local input="$1"
  if [[ -f "$input/WORKFLOW_COMPLETE.json" ]]; then
    printf '%s\n' "$input"
    return
  fi
  local completions=()
  while IFS= read -r marker; do completions+=("$marker"); done < <(
    find "$input" -mindepth 2 -maxdepth 2 -type f -name WORKFLOW_COMPLETE.json | sort
  )
  [[ "${#completions[@]}" -eq 1 ]] || \
    die "$input must contain exactly one completed PFP run; found ${#completions[@]}"
  printf '%s\n' "${completions[0]%/WORKFLOW_COMPLETE.json}"
}

declare -A RUNS=()
declare -A PREDICTION_RUNS=()
for mode in full sequence-only sequence-text sequence-structure sequence-ppi; do
  RUNS[$mode]="$(resolve_completed_run "${RUN_INPUTS[$mode]}")"
  PREDICTION_RUNS[$mode]="$(resolve_completed_run "${PREDICTION_INPUTS[$mode]}")"
  [[ -f "${RUNS[$mode]}/reports/run_report.json" ]] || \
    die "Run report is missing for $mode: ${RUNS[$mode]}"
  [[ -f "${PREDICTION_RUNS[$mode]}/evaluation/prediction_artifacts/prediction_artifact_manifest.json" ]] || \
    die "Prediction artifacts are missing for $mode: ${PREDICTION_RUNS[$mode]}"
done

JOB_TOKEN="${JOB_ID:-manual}"
WORK="/scratch0/pfp_modality_panel_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
RESULTS_STAGE="$WORK/results"
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
  [[ ! -d "$WORK" || -L "$WORK" || "$WORK" != /scratch0/pfp_modality_panel_* ]] || \
    rm -rf -- "$WORK"
  exit "$status"
}
trap cleanup EXIT

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
[[ ! -e "$PUBLISH_STAGE" ]] || die "Publication stage already exists: $PUBLISH_STAGE"
[[ ! -e "$PUBLISH_LOCK" ]] || die "Publication lock already exists: $PUBLISH_LOCK"
mkdir -p "$WORK" "$RESULTS_STAGE" "$(dirname "$OUTPUT_DIR")"

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from the framework checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || \
    die "Submission checkout has uncommitted changes"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be a full commit"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || \
  die "Cloned framework commit differs from requested revision"
[[ -z "$(git_in_dir "$FRAMEWORK_DIR" status --porcelain)" ]] || \
  die "Cloned framework checkout is not clean"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
add_mmfp_singularity_bind "$(dirname "$OUTPUT_DIR")"
add_mmfp_singularity_bind "$(dirname "$OBO_FILE")"
for mode in full sequence-only sequence-text sequence-structure sequence-ppi; do
  add_mmfp_singularity_bind "${RUNS[$mode]}"
  add_mmfp_singularity_bind "${PREDICTION_RUNS[$mode]}"
done
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

SENSITIVITY_REPORTS=()
RUN_REPORTS=()
PREDICTION_MANIFESTS=()
for mode in sequence-only sequence-text sequence-structure sequence-ppi full; do
  RUN_REPORTS+=(--run-report "${RUNS[$mode]}/reports/run_report.json")
  PREDICTION_MANIFESTS+=(
    --prediction-manifest "$mode=${PREDICTION_RUNS[$mode]}/evaluation/prediction_artifacts/prediction_artifact_manifest.json"
  )
done

echo "==> Compare canonical modality results"
COMPARISON_COMMAND=(
  "$PYTHON_BIN" scripts/diagnostics/compare_pfp_modality_runs.py
  "${RUN_REPORTS[@]}"
  "${PREDICTION_MANIFESTS[@]}"
  --output-dir "$RESULTS_STAGE/canonical_comparison"
)
if [[ "$ALLOW_FRAMEWORK_COMMIT_DRIFT" == "1" ]]; then
  COMPARISON_COMMAND+=(--allow-framework-commit-drift)
fi
"${COMPARISON_COMMAND[@]}"

for mode in sequence-only sequence-text sequence-structure sequence-ppi full; do
  echo "==> Root-only sensitivity: $mode"
  "$PYTHON_BIN" scripts/diagnostics/evaluate_pfp_label_sensitivity.py \
    --prediction-manifest "${PREDICTION_RUNS[$mode]}/evaluation/prediction_artifacts/prediction_artifact_manifest.json" \
    --obo-file "$OBO_FILE" \
    --output-dir "$RESULTS_STAGE/sensitivity/$mode"
  SENSITIVITY_REPORTS+=(
    --report "$RESULTS_STAGE/sensitivity/$mode/root_exclusion_sensitivity.json"
  )
done

echo "==> Compare prediction sensitivities"
"$PYTHON_BIN" scripts/diagnostics/compare_pfp_label_sensitivity.py \
  "${SENSITIVITY_REPORTS[@]}" \
  --output-dir "$RESULTS_STAGE/sensitivity_comparison"

mkdir -p "$PUBLISH_STAGE"
cp -a "$RESULTS_STAGE/." "$PUBLISH_STAGE/"
"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py write \
  --root "$PUBLISH_STAGE" --include-nested-control-files
MANIFEST_SHA256="$(
  "$PYTHON_BIN" -c 'import hashlib,pathlib,sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' \
    "$PUBLISH_STAGE/output_manifest.json"
)"
"$PYTHON_BIN" -c 'import json,pathlib,sys; pathlib.Path(sys.argv[1]).write_text(json.dumps({"complete":True,"manifest":"output_manifest.json","manifest_sha256":sys.argv[2]},indent=2)+"\n")' \
  "$PUBLISH_STAGE/WORKFLOW_COMPLETE.json" "$MANIFEST_SHA256"
"$PYTHON_BIN" scripts/model_execution/manage_output_manifest.py verify \
  --root "$PUBLISH_STAGE" --include-nested-control-files
mkdir "$PUBLISH_LOCK"
LOCK_HELD=1
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory appeared during analysis: $OUTPUT_DIR"
mv -T "$PUBLISH_STAGE" "$OUTPUT_DIR"
rmdir "$PUBLISH_LOCK"
LOCK_HELD=0
echo "Published modality panel analysis: $OUTPUT_DIR"
