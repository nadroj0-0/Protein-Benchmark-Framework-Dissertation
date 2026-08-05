#!/usr/bin/env bash
# Compare complete IA/Xu specificity analyses on UCL Grid Engine.

#$ -S /bin/bash
#$ -l tmem=8G
#$ -l scratch0free=4G
#$ -l tscratch=4G
#$ -l h_rt=2:0:0
#$ -pe smp 1
#$ -j y
#$ -N spec_3way
#$ -V
#$ -notify

set -Eeuo pipefail
export PYTHONDONTWRITEBYTECODE=1

usage() {
  cat <<'EOF'
Usage: qsub hpc_jobs/active/hpc_specificity_comparison.sh \
  --source cafa3=/SAN/.../specificity \
  --source global-nk=/SAN/.../specificity \
  --source nk-lk=/SAN/.../specificity \
  --output-dir /SAN/.../three-way-comparison
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
require_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

SOURCES=()
OUTPUT_DIR=""
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source) require_value "$@"; SOURCES+=("$2"); shift 2 ;;
    --output-dir) require_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "${#SOURCES[@]}" -ge 2 ]] || die "At least two --source values are required"
[[ "$OUTPUT_DIR" == /SAN/* ]] || die "--output-dir must be an absolute SAN path"
[[ ! -e "$OUTPUT_DIR" ]] || die "Output directory already exists: $OUTPUT_DIR"
for source in "${SOURCES[@]}"; do
  [[ "$source" == *=/SAN/* ]] || die "Each --source must be LABEL=/SAN/path"
done

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/specificity_comparison_${JOB_TOKEN}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  [[ ! -d "$WORK" || -L "$WORK" || "$WORK" != /scratch0/specificity_comparison_* ]] || rm -rf -- "$WORK"
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK" "$(dirname "$OUTPUT_DIR")"
if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Submit from a clean framework checkout"
  [[ -z "$(git_in_dir "$SUBMISSION_DIR" status --porcelain)" ]] || die "Submission checkout has uncommitted changes"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be a full commit"

echo "Host             : $(hostname)"
echo "Job ID           : ${JOB_ID:-manual}"
echo "Framework commit : $FRAMEWORK_COMMIT"
echo "SAN output       : $OUTPUT_DIR"
printf 'Source           : %s\n' "${SOURCES[@]}"
echo "Started          : $(date -Is)"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || die "Framework checkout differs from requested commit"

ARGS=()
for source in "${SOURCES[@]}"; do
  ARGS+=(--source "$source")
done
cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind "$WORK"
add_mmfp_singularity_bind "$(dirname "$OUTPUT_DIR")"
for source in "${SOURCES[@]}"; do
  add_mmfp_singularity_bind "${source#*=}"
done
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

"$PYTHON_BIN" scripts/diagnostics/compare_pfp_specificity_runs.py \
  "${ARGS[@]}" \
  --output-dir "$OUTPUT_DIR"

[[ -f "$OUTPUT_DIR/RUN_COMPLETE.json" ]] || die "Completion marker was not published"
[[ -f "$OUTPUT_DIR/output_manifest.json" ]] || die "Output manifest was not published"
echo "Completed         : $(date -Is)"
echo "Published results : $OUTPUT_DIR"
