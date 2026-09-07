#!/usr/bin/env bash
# Build one of the two paired random controls on UCL Grid Engine.

#$ -l tmem=64G
#$ -l tscratch=20G
#$ -l scratch0free=20G
#$ -l h_rt=72:0:0
#$ -j y
#$ -N random_baseline
#$ -V

set -euo pipefail

MODE="${RANDOM_BASELINE_MODE:-}"
case "$MODE" in
    random-h|random-t) ;;
    *) echo "RANDOM_BASELINE_MODE must be random-h or random-t" >&2; exit 2 ;;
esac

FRAMEWORK_REPO_URL="${FRAMEWORK_REPO_URL:-https://github.com/nadroj0-0/Protein-Benchmark-Framework-Dissertation.git}"
FRAMEWORK_REVISION="${FRAMEWORK_REVISION:-}"
[[ "$FRAMEWORK_REVISION" =~ ^[0-9a-f]{40}$ ]] || {
    echo "FRAMEWORK_REVISION must be an exact 40-character commit" >&2
    exit 2
}

GLOBAL_SOURCE="${GLOBAL_SOURCE:-}"
HOMOLOGY_SOURCE="${HOMOLOGY_SOURCE:-}"
GO_T0="${GO_T0:-}"
GO_T1="${GO_T1:-}"
UNIPROT_SPROT_T1="${UNIPROT_SPROT_T1:-}"
UNIPROT_TREMBL_T1="${UNIPROT_TREMBL_T1:-}"
GOA_T1="${GOA_T1:-}"
RESULTS_ROOT="${RESULTS_ROOT:-}"

[[ -n "$RESULTS_ROOT" && "$RESULTS_ROOT" == /SAN/* ]] || {
    echo "RESULTS_ROOT must be an absolute SAN path" >&2
    exit 2
}

if [[ "$MODE" == "random-h" ]]; then
    SOURCE_DIR="$HOMOLOGY_SOURCE"
    INPUT_MANIFEST="$HOMOLOGY_SOURCE/input_manifest.json"
    OUTPUT_DIR="$RESULTS_ROOT/random_h_homology30"
    BENCHMARK_ID="random-h-evidence12-homology30"
    REQUIRED_PATHS=(
        "$HOMOLOGY_SOURCE"
        "$HOMOLOGY_SOURCE/RUN_COMPLETE.json"
        "$HOMOLOGY_SOURCE/validation_report.json"
        "$INPUT_MANIFEST"
        "$GO_T1"
    )
else
    SOURCE_DIR="$GLOBAL_SOURCE"
    INPUT_MANIFEST="${GLOBAL_MANIFEST:-}"
    OUTPUT_DIR="$RESULTS_ROOT/random_t_t1"
    BENCHMARK_ID="random-t-evidence12-t1"
    REQUIRED_PATHS=(
        "$GLOBAL_SOURCE"
        "$INPUT_MANIFEST"
        "$GO_T0"
        "$GO_T1"
        "$UNIPROT_SPROT_T1"
        "$UNIPROT_TREMBL_T1"
        "$GOA_T1"
    )
fi

for path in "${REQUIRED_PATHS[@]}"; do
    [[ -e "$path" ]] || { echo "Required SAN input is missing: $path" >&2; exit 1; }
done
for filename in train_data_train.pkl train_data_valid.pkl test_data.pkl terms.pkl; do
    [[ -s "$SOURCE_DIR/$filename" ]] || {
        echo "Required source benchmark file is missing or empty: $SOURCE_DIR/$filename" >&2
        exit 1
    }
done
for aspect in bp cc mf; do
    for split in training validation test; do
        [[ -s "$SOURCE_DIR/$aspect-$split.csv" ]] || {
            echo "Required source CSV is missing or empty: $SOURCE_DIR/$aspect-$split.csv" >&2
            exit 1
        }
    done
done
[[ ! -e "$OUTPUT_DIR" ]] || { echo "Output already exists: $OUTPUT_DIR" >&2; exit 1; }

WORK="/scratch0/random_baseline_${JOB_ID:-manual}_${MODE}"
FRAMEWORK_DIR="$WORK/framework"
cleanup() {
    status=$?
    cd "$HOME"
    rm -rf "$WORK"
    exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

mkdir -p "$WORK" "$RESULTS_ROOT"
git clone "$FRAMEWORK_REPO_URL" "$FRAMEWORK_DIR"
git --git-dir="$FRAMEWORK_DIR/.git" --work-tree="$FRAMEWORK_DIR" checkout "$FRAMEWORK_REVISION"
ACTUAL_REVISION="$(git --git-dir="$FRAMEWORK_DIR/.git" --work-tree="$FRAMEWORK_DIR" rev-parse HEAD)"
[[ "$ACTUAL_REVISION" == "$FRAMEWORK_REVISION" ]] || {
    echo "Framework checkout mismatch: $ACTUAL_REVISION" >&2
    exit 1
}

cd "$FRAMEWORK_DIR"
# shellcheck disable=SC1091
source scripts/reproduction_common.sh
load_framework_paths "$FRAMEWORK_DIR"
activate_or_create_mmfp_env
export PYTHONPATH="$FRAMEWORK_DIR/benchmark_builders/contemporary_cafa/src${PYTHONPATH:+:$PYTHONPATH}"
export CAFA_BUILDER_GOA_PROGRESS_INTERVAL="${CAFA_BUILDER_GOA_PROGRESS_INTERVAL:-1000000}"

COMMAND=(
    python -m cafa_benchmark_builder.random_baseline
    --benchmark-id "$BENCHMARK_ID"
    --source-benchmark-dir "$SOURCE_DIR"
    --output-dir "$OUTPUT_DIR"
    --seed 0
    --input-manifest "$INPUT_MANIFEST"
)
if [[ "$MODE" == "random-h" ]]; then
    COMMAND+=(--mode prepared --go-obo "$GO_T1")
else
    COMMAND+=(
        --mode t1-snapshot
        --go-obo "$GO_T0"
        --source-go-obo "$GO_T1"
        --uniprot "$UNIPROT_SPROT_T1"
        --uniprot "$UNIPROT_TREMBL_T1"
        --goa "$GOA_T1"
    )
    for evidence_code in EXP IDA IPI IMP IGI IEP HTP HDA HMP HGI HEP IGC; do
        COMMAND+=(--evidence-code "$evidence_code")
    done
fi

echo "Host: $(hostname)"
echo "Mode: $MODE"
echo "Framework revision: $FRAMEWORK_REVISION"
echo "Source benchmark: $SOURCE_DIR"
echo "Output: $OUTPUT_DIR"
"${COMMAND[@]}"

[[ -s "$OUTPUT_DIR/RUN_COMPLETE.json" ]] || {
    echo "Builder returned without RUN_COMPLETE.json" >&2
    exit 1
}
echo "Completed: $OUTPUT_DIR"
