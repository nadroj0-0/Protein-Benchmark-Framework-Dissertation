#!/usr/bin/env bash

#$ -l tmem=48G
#$ -l tscratch=200G
#$ -l scratch0free=200G
#$ -l h_rt=72:0:0
#$ -j y
#$ -N cafa3_org_replay
#$ -V

set -Eeuo pipefail

die() { echo "ERROR: $*" >&2; exit 2; }
git_in_dir() { local directory="$1"; shift; (cd "$directory" && git "$@"); }

: "${FRAMEWORK_COMMIT:?Submit with the exact pushed FRAMEWORK_COMMIT}"
: "${RESULTS_ROOT:?Submit with an absolute RESULTS_ROOT}"

[[ "$FRAMEWORK_COMMIT" =~ ^[0-9a-f]{40}$ ]] || die "Invalid FRAMEWORK_COMMIT"
[[ "$RESULTS_ROOT" == /* && "$RESULTS_ROOT" != / ]] || die "RESULTS_ROOT must be absolute"

JOB_KEY="${JOB_ID:-manual}"
RUN_TAG="${RUN_TAG:-${JOB_KEY}_$(date -u +%Y%m%dT%H%M%SZ)}"
[[ "$RUN_TAG" =~ ^[A-Za-z0-9._-]+$ && "$RUN_TAG" =~ [A-Za-z0-9] ]] || die "Unsafe RUN_TAG"

WORK="/scratch0/cafa3_organizer_replay_${JOB_KEY}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
STAGING="$RESULTS_ROOT/.${RUN_TAG}.staging-${JOB_KEY}"
FINAL="$RESULTS_ROOT/$RUN_TAG"
FAILED="$RESULTS_ROOT/${RUN_TAG}.failed"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
STATUS=1

cleanup() {
  local status=$?
  trap - EXIT
  set +e
  if [[ -d "$WORK" && ! -L "$WORK" && "$WORK" == /scratch0/cafa3_organizer_replay_* ]]; then
    cd "$HOME"
    rm -rf -- "$WORK"
  fi
  exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

for path in "$WORK" "$STAGING" "$FINAL" "$FAILED"; do
  [[ ! -e "$path" ]] || die "Refusing to overwrite existing path: $path"
done
mkdir -p "$WORK" "$RESULTS_ROOT"

git clone --no-checkout "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git_in_dir "$FRAMEWORK_DIR" checkout --detach "$FRAMEWORK_COMMIT"
[[ "$(git_in_dir "$FRAMEWORK_DIR" rev-parse HEAD)" == "$FRAMEWORK_COMMIT" ]] || \
  die "Framework checkout does not match the requested commit"

cd "$FRAMEWORK_DIR"
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"

export CAFA3_ORGANIZER_REPLAY_MATRIX=1
export HISTORICAL_TRAINING_SNAPSHOT=february-2017-legacy
export TARGET_UNIVERSE_POLICY=official-cafa3-targets
export HISTORICAL_TEST_SOURCE=raw-goa
export HISTORICAL_T1_ENDPOINT_POLICY=assigned-date-proxy
export HISTORICAL_BACKFILL_POLICY=exclude-pre-t0
export HISTORICAL_BENCHMARK_ONTOLOGY=february-go
export INCLUDE_TREMBL_TARGETS=1
export DECOMPRESS_GOA=1
export USE_PIGZ=1
export KEEP_SCRATCH=0
export PYTHON_BIN
export SCRATCH_BASE="$WORK"
export TIMESTAMP="$RUN_TAG"
export REPORT_COPY_DIR="$STAGING"

echo "Framework commit : $FRAMEWORK_COMMIT"
echo "Results root     : $RESULTS_ROOT"
echo "Run tag          : $RUN_TAG"
echo "Started at       : $(date -u +%Y-%m-%dT%H:%M:%SZ)"

if bash scripts/validation/run_cafa3_historical_validation.sh; then
  [[ -f "$STAGING/generated/organizer_replay/REPLAY_MATRIX_COMPLETE.json" ]] || \
    die "Replay script returned success without a matrix completion marker"
  mv "$STAGING" "$FINAL"
  STATUS=0
  echo "Published replay matrix: $FINAL"
else
  status=$?
  if [[ -e "$STAGING" ]]; then mv "$STAGING" "$FAILED"; fi
  echo "Replay matrix failed; retained diagnostics at: $FAILED" >&2
  exit "$status"
fi

exit "$STATUS"
