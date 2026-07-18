#!/usr/bin/env bash
# UCL Grid Engine wrapper for persistent SAN frozen-input acquisition.

#$ -l tmem=8G
#$ -l tscratch=120G
#$ -l scratch0free=120G
#$ -l h_rt=72:0:0
#$ -j y
#$ -N san_frozen_inputs
#$ -V

set -euo pipefail

WORK="/scratch0/san_frozen_inputs_${JOB_ID:-manual}"
FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_DIR="$WORK/Protein-Benchmark-Framework-Dissertation"
SAN_ROOT="${SAN_ROOT:-/SAN/bioinf/bmpfp}"
SAN_INPUT_PROFILES="${SAN_INPUT_PROFILES:-all}"
SAN_INPUT_RESERVE_GB="${SAN_INPUT_RESERVE_GB:-40}"
HOMOLOGY_CACHE_WORK_DIR="${HOMOLOGY_CACHE_WORK_DIR:-$WORK/homology-common-cache}"

cleanup() {
    local status=$?
    set +e
    echo "Cleaning job-owned scratch directory: $WORK"
    cd "$HOME"
    rm -rf "$WORK"
    exit "$status"
}
trap cleanup EXIT
trap 'echo "Received kill signal"; exit 130' SIGINT SIGTERM

echo "Host        : $(hostname)"
echo "Job ID      : ${JOB_ID:-manual}"
echo "SAN root    : $SAN_ROOT"
echo "Profiles    : $SAN_INPUT_PROFILES"
echo "Reserve GiB : $SAN_INPUT_RESERVE_GB"
echo "Scratch     : $WORK"
echo "Cache work  : $HOMOLOGY_CACHE_WORK_DIR"

[[ -d "$SAN_ROOT" && -w "$SAN_ROOT" ]] || {
    echo "SAN root is unavailable or not writable: $SAN_ROOT" >&2
    exit 1
}

mkdir -p "$WORK"
git clone "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
cd "$FRAMEWORK_DIR"
echo "Framework revision: $(git rev-parse HEAD)"

PROFILE_ARGUMENTS=()
IFS=',' read -r -a profile_values <<< "$SAN_INPUT_PROFILES"
for profile in "${profile_values[@]}"; do
    [[ -n "$profile" ]] || {
        echo "SAN_INPUT_PROFILES contains an empty profile" >&2
        exit 2
    }
    PROFILE_ARGUMENTS+=(--profile "$profile")
done

HOMOLOGY_CACHE_WORK_DIR="$HOMOLOGY_CACHE_WORK_DIR" \
bash scripts/data_acquisition/populate_san_frozen_inputs.sh \
    --root "$SAN_ROOT" \
    --reserve-gb "$SAN_INPUT_RESERVE_GB" \
    "${PROFILE_ARGUMENTS[@]}" \
    "$@"

echo "Finished successfully: $(date)"
