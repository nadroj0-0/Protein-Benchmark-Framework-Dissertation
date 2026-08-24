#!/usr/bin/env bash
# Analyse one completed NK+LK model against its nested global-NK comparator.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=32G
#$ -l tscratch=12G
#$ -l scratch0free=12G
#$ -l h_rt=24:0:0
#$ -pe smp 2
#$ -N nk_lk_analysis
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

PREDICTION_MANIFEST=""
MEMBERSHIP_TSV=""
NK_EVALUATION_SUMMARY=""
OUTPUT_DIR=""
BOOTSTRAP_REPLICATES=2000
BOOTSTRAP_SEED=20260805
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prediction-manifest) require_value "$@"; PREDICTION_MANIFEST="$2"; shift 2 ;;
    --membership-tsv) require_value "$@"; MEMBERSHIP_TSV="$2"; shift 2 ;;
    --nk-evaluation-summary) require_value "$@"; NK_EVALUATION_SUMMARY="$2"; shift 2 ;;
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    --bootstrap-replicates) require_value "$@"; BOOTSTRAP_REPLICATES="$2"; shift 2 ;;
    --bootstrap-seed) require_value "$@"; BOOTSTRAP_SEED="$2"; shift 2 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

for value in PREDICTION_MANIFEST MEMBERSHIP_TSV NK_EVALUATION_SUMMARY OUTPUT_DIR; do
  [[ -n "${!value}" ]] || die "$value is required"
  [[ "${!value}" == /SAN/* ]] || die "$value must be an absolute SAN path"
done
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"
for path in "$PREDICTION_MANIFEST" "$MEMBERSHIP_TSV" "$NK_EVALUATION_SUMMARY"; do
  [[ -f "$path" && ! -L "$path" ]] || die "Required input is missing or unsafe: $path"
done

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/nk_lk_completed_analysis_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
SCRATCH_OUTPUT="$WORK/output"
STAGING="${OUTPUT_DIR}.staging-${JOB_TOKEN}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  [[ ! -d "$STAGING" || -L "$STAGING" ]] || rm -rf -- "$STAGING"
  [[ ! -d "$WORK" || -L "$WORK" || "$WORK" != /scratch0/nk_lk_completed_analysis_* ]] || rm -rf -- "$WORK"
  exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

[[ ! -e "$WORK" && ! -e "$STAGING" ]] || die "Scratch or staging path already exists"
mkdir -p "$WORK" "$(dirname "$OUTPUT_DIR")"
if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a framework checkout or pass FRAMEWORK_COMMIT"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout is dirty"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
add_mmfp_singularity_bind /SAN/bioinf/bmpfp
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

"$PYTHON_BIN" scripts/diagnostics/evaluate_nk_lk_completed_run.py \
  --framework-root "$FRAMEWORK_DIR" \
  --prediction-manifest "$PREDICTION_MANIFEST" \
  --membership-tsv "$MEMBERSHIP_TSV" \
  --nk-evaluation-summary "$NK_EVALUATION_SUMMARY" \
  --output-dir "$SCRATCH_OUTPUT" \
  --bootstrap-replicates "$BOOTSTRAP_REPLICATES" \
  --bootstrap-seed "$BOOTSTRAP_SEED"

[[ -f "$SCRATCH_OUTPUT/ANALYSIS_COMPLETE.json" ]] || die "Analysis completion marker is missing"
cp -a "$SCRATCH_OUTPUT" "$STAGING"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory appeared during publication"
mv "$STAGING" "$OUTPUT_DIR"
echo "Published NK+LK completed-run analysis: $OUTPUT_DIR"
