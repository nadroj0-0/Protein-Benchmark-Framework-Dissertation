#!/usr/bin/env bash
# Bind one hydrated PFP cache archive to the exact selected benchmark CSVs.

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l tmem=8G
#$ -l tscratch=4G
#$ -l scratch0free=8G
#$ -l h_rt=12:0:0
#$ -pe smp 1
#$ -N emb_bind
#$ -V
#$ -notify

set -Eeuo pipefail

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

: "${BENCHMARK_DIR:?Pass BENCHMARK_DIR with qsub -v}"
: "${BENCHMARK_ID:?Pass BENCHMARK_ID with qsub -v}"
: "${EMBEDDING_ARCHIVE:?Pass EMBEDDING_ARCHIVE with qsub -v}"
: "${RUN_CONFIG:?Pass RUN_CONFIG with qsub -v}"
: "${EVIDENCE_OUTPUT:?Pass EVIDENCE_OUTPUT with qsub -v}"

JOB_TOKEN="${JOB_ID:-manual_$$}"
WORK="/scratch0/embedding_evidence_bind_${JOB_TOKEN}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_COMMIT="${FRAMEWORK_COMMIT:-}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
SCRATCH_EVIDENCE="$WORK/evidence"
SUBMISSION_DIR="${SGE_O_WORKDIR:-$PWD}"
WORK_OWNED=0

wait_for_file() {
  local path="$1" attempts="${2:-30}" delay="${3:-10}"
  local attempt
  for ((attempt=1; attempt<=attempts; attempt++)); do
    [[ -f "$path" ]] && return 0
    echo "Waiting for SAN file ($attempt/$attempts): $path"
    sleep "$delay"
  done
  return 1
}

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  if [[ "$WORK_OWNED" == 1 && "$WORK" == /scratch0/embedding_evidence_bind_* && ! -L "$WORK" ]]; then
    cd "$HOME"
    rm -rf -- "$WORK"
  fi
  exit "$status"
}
trap cleanup EXIT
trap 'echo "Received termination signal"; exit 130' INT TERM

[[ "$EVIDENCE_OUTPUT" == /SAN/* ]] || die "EVIDENCE_OUTPUT must be on SAN"
[[ ! -e "$EVIDENCE_OUTPUT" ]] || die "Evidence output already exists: $EVIDENCE_OUTPUT"
wait_for_file "$EMBEDDING_ARCHIVE" || die "Embedding archive did not become visible"
[[ -d "$BENCHMARK_DIR" ]] || die "Benchmark directory is missing: $BENCHMARK_DIR"
[[ ! -e "$WORK" ]] || die "Scratch path already exists: $WORK"
mkdir -p "$WORK" "$(dirname "$EVIDENCE_OUTPUT")"
WORK_OWNED=1

if [[ -z "$FRAMEWORK_COMMIT" ]]; then
  [[ -d "$SUBMISSION_DIR/.git" ]] || die "Pass FRAMEWORK_COMMIT outside a framework checkout"
  FRAMEWORK_COMMIT="$(git_in_dir "$SUBMISSION_DIR" rev-parse HEAD)"
fi
[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || die "FRAMEWORK_COMMIT must be complete"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
add_mmfp_singularity_bind /SAN/bioinf/bmpfp
add_mmfp_singularity_bind "$WORK"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

CONFIG_PATH="$RUN_CONFIG"
[[ "$CONFIG_PATH" == /* ]] || CONFIG_PATH="$FRAMEWORK_DIR/$CONFIG_PATH"
COMMAND=(
  "$PYTHON_BIN" scripts/embeddings/bind_embedding_archive_evidence.py
  --benchmark-dir "$BENCHMARK_DIR"
  --benchmark-id "$BENCHMARK_ID"
  --archive "$EMBEDDING_ARCHIVE"
  --config "$CONFIG_PATH"
  --framework-commit "$FRAMEWORK_COMMIT"
  --output-dir "$SCRATCH_EVIDENCE"
)
if [[ -n "${EMBEDDING_ARCHIVE_SHA256:-}" ]]; then
  COMMAND+=(--archive-sha256 "$EMBEDDING_ARCHIVE_SHA256")
fi
"${COMMAND[@]}"

[[ -f "$SCRATCH_EVIDENCE/RUN_COMPLETE.json" ]] || die "Evidence binder did not complete"
STAGING="${EVIDENCE_OUTPUT}.staging-${JOB_TOKEN}"
[[ ! -e "$STAGING" ]] || die "Evidence staging path already exists: $STAGING"
cp -a "$SCRATCH_EVIDENCE" "$STAGING"
mv "$STAGING" "$EVIDENCE_OUTPUT"
echo "Published strict embedding evidence: $EVIDENCE_OUTPUT"
