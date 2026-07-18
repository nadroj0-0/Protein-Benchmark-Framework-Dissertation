#!/usr/bin/env bash
# Real shell entrypoint for Daniel's frozen UniRef90/MMseqs2 homology benchmark.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRAMEWORK_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILDER_ROOT="$FRAMEWORK_ROOT/benchmark_builders/homology_cluster"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"
artifact_catalog_configure "$FRAMEWORK_ROOT" "${ARTIFACT_CATALOG:-}"

PYTHON_BIN="${PYTHON_BIN:-python3}"
IDENTITY="${IDENTITY:-30}"
SPLIT_POLICY="${SPLIT_POLICY:-sequence-balanced}"
TRAINING_POPULATION="${TRAINING_POPULATION:-annotated-only}"
UNIPROT_SOURCE_SCOPE="${UNIPROT_SOURCE_SCOPE:-}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$FRAMEWORK_ROOT/results/homology_cluster_benchmark}"
TEMP_DIR="${TEMP_DIR:-${TMPDIR:-/tmp}/homology-cluster-benchmark}"
THREADS="${THREADS:-1}"
SEED="${SEED:-0}"
MIN_COUNT="${MIN_COUNT:-50}"
MMSEQS_BIN="${MMSEQS_BIN:-mmseqs}"
EXPECTED_MMSEQS_VERSION="${EXPECTED_MMSEQS_VERSION:-}"
FROZEN_INPUT_MANIFEST="${FROZEN_INPUT_MANIFEST:-}"
ATTRITION_POLICY="${ATTRITION_POLICY:-}"
ATTRITION_OVERRIDE="${ATTRITION_OVERRIDE:-}"
FRAMEWORK_REVISION="${FRAMEWORK_REVISION:-}"
DIAGNOSTIC_PILOT="${DIAGNOSTIC_PILOT:-0}"
RUN_ID="${RUN_ID:-local}"
REQUESTED_SLOTS="${REQUESTED_SLOTS:-}"
ALLOCATED_SLOTS="${ALLOCATED_SLOTS:-}"
UNIPROT_RELEASE="${UNIPROT_RELEASE:-2026_02}"
GOA_RELEASE="${GOA_RELEASE:-234}"
ONTOLOGY_RELEASE="${ONTOLOGY_RELEASE:-releases/2026-06-15}"
GO_OBO_URL="${GO_OBO_URL:-https://release.geneontology.org/2026-06-19/ontology/go-basic.obo}"
NO_DOWNLOADS="${NO_DOWNLOADS:-0}"
DRY_RUN="${DRY_RUN:-0}"
KEEP_TEMP="${KEEP_TEMP:-0}"
FIXTURE_MODE="${FIXTURE_MODE:-0}"
SCRATCH_SAFETY_MULTIPLIER="${SCRATCH_SAFETY_MULTIPLIER:-1}"
MINIMUM_FREE_DISK_GB="${MINIMUM_FREE_DISK_GB:-0}"
PERSISTENT_RESULTS_ROOT="${PERSISTENT_RESULTS_ROOT:-}"
MMSEQS_WORK_MULTIPLIER="${MMSEQS_WORK_MULTIPLIER:-1}"
PUBLICATION_SAFETY_MULTIPLIER="${PUBLICATION_SAFETY_MULTIPLIER:-1}"
EXCLUDED_SAMPLE_PER_REASON="${EXCLUDED_SAMPLE_PER_REASON:-1000}"
LOG_FILE="${LOG_FILE:-}"
ACTIVE_CHILD_PID=""
SIGNAL_STATUS=0
SIGNAL_WATCHDOG_PID=""
SIGNAL_GRACE_SECONDS=10

UNIREF90_FASTA="$(resolve_artifact_path uniref90_t1 "${UNIREF90_FASTA:-}" || true)"
GO_OBO="$(resolve_artifact_path go_basic_t1 "${GO_OBO:-}" || true)"
resolve_common_preprocessing_cache() {
    local explicit="${1:-}"
    local candidate=""
    if [[ -n "$explicit" ]]; then
        if [[ -d "$explicit" && -s "$explicit/CACHE_COMPLETE.json" ]]; then
            (cd "$explicit" && pwd -P)
            return 0
        fi
        if [[ -s "$explicit" && "$(basename "$explicit")" == "CACHE_COMPLETE.json" ]]; then
            (cd "$(dirname "$explicit")" && pwd -P)
            return 0
        fi
        artifact_catalog_warn \
            "explicit homology common cache is incomplete; trying catalogue/raw fallback: $explicit"
    fi
    candidate="$(artifact_catalog_lookup homology_common_preprocessing_2026_02 2>/dev/null || true)"
    if [[ -s "$candidate" && "$(basename "$candidate")" == "CACHE_COMPLETE.json" ]]; then
        (cd "$(dirname "$candidate")" && pwd -P)
        return 0
    fi
    return 1
}
HOMOLOGY_COMMON_PREPROCESSING_CACHE="$(
    resolve_common_preprocessing_cache "${HOMOLOGY_COMMON_PREPROCESSING_CACHE:-}" || true
)"
if [[ -n "$HOMOLOGY_COMMON_PREPROCESSING_CACHE" ]]; then
    IDMAPPING="${IDMAPPING:-}"
    UNIPROT_SPROT_SEQUENCES="${UNIPROT_SPROT_SEQUENCES:-}"
    UNIPROT_TREMBL_SEQUENCES="${UNIPROT_TREMBL_SEQUENCES:-}"
    GOA="${GOA:-}"
else
    IDMAPPING="$(resolve_artifact_path idmapping_t1 "${IDMAPPING:-}" || true)"
    UNIPROT_SPROT_SEQUENCES="$(resolve_artifact_path uniprot_sprot_t1 "${UNIPROT_SPROT_SEQUENCES:-}" || true)"
    UNIPROT_TREMBL_SEQUENCES="$(resolve_artifact_path uniprot_trembl_t1 "${UNIPROT_TREMBL_SEQUENCES:-}" || true)"
    GOA="$(resolve_artifact_path goa_t1 "${GOA:-}" || true)"
fi
for catalog_input in "$UNIREF90_FASTA" "$IDMAPPING" "$UNIPROT_SPROT_SEQUENCES" \
    "$UNIPROT_TREMBL_SEQUENCES" "$GOA" "$GO_OBO" \
    "$HOMOLOGY_COMMON_PREPROCESSING_CACHE"; do
    [[ -z "$catalog_input" ]] || add_mmfp_singularity_bind "$(dirname "$catalog_input")"
done

forward_signal() {
    local signal_name="$1"
    local status="$2"
    SIGNAL_STATUS="$status"
    echo "Received $signal_name; forwarding to active benchmark process" >&2
    if [[ -n "$ACTIVE_CHILD_PID" ]]; then
        kill -s "$signal_name" "$ACTIVE_CHILD_PID" 2>/dev/null || true
        # Bash 3.2 asynchronous children can inherit ignored SIGINT. TERM is the portable
        # termination mechanism; this wrapper still exits 130 for an incoming INT.
        if [[ "$signal_name" == "INT" ]]; then
            kill -s TERM "$ACTIVE_CHILD_PID" 2>/dev/null || true
        fi
        local child_pid="$ACTIVE_CHILD_PID"
        (
            sleep "$SIGNAL_GRACE_SECONDS"
            kill -s KILL "$child_pid" 2>/dev/null || true
        ) &
        SIGNAL_WATCHDOG_PID=$!
    fi
}
trap 'forward_signal INT 130' INT
trap 'forward_signal TERM 143' TERM

require_local_or_url() {
    local label="$1"
    local local_path="$2"
    local source_url="$3"
    if [[ -n "$local_path" ]]; then
        if [[ ! -f "$local_path" ]]; then
            echo "Configured local $label does not exist: $local_path" >&2
            exit 1
        fi
        return
    fi
    if [[ -z "$source_url" ]]; then
        echo "Supply either the local path or frozen URL for $label" >&2
        exit 1
    fi
}

append_input() {
    local option="$1"
    local local_path="$2"
    local source_url="$3"
    local expected_sha="$4"
    if [[ -n "$local_path" ]]; then
        COMMAND+=("--$option" "$local_path")
    fi
    if [[ -n "$source_url" ]]; then
        COMMAND+=("--$option-url" "$source_url")
    fi
    if [[ -n "$expected_sha" ]]; then
        COMMAND+=("--$option-sha256" "$expected_sha")
    fi
}

if [[ "$TRAINING_POPULATION" != "annotated-only" ]]; then
    echo "TRAINING_POPULATION=$TRAINING_POPULATION is unsupported; only annotated-only is implemented." >&2
    echo "No zero-negative, homology-transfer, or representative-label policy is authorized." >&2
    exit 2
fi
if [[ -z "$UNIPROT_SOURCE_SCOPE" && "$FIXTURE_MODE" == "1" ]]; then
    UNIPROT_SOURCE_SCOPE="sprot-only"
fi
case "$UNIPROT_SOURCE_SCOPE" in
    sprot-only|trembl-only|sprot-and-trembl) ;;
    *) echo "UNIPROT_SOURCE_SCOPE must be explicitly set to sprot-only, trembl-only, or sprot-and-trembl" >&2; exit 2 ;;
esac

if [[ "$DRY_RUN" != "1" ]]; then
    require_local_or_url UniRef90 "${UNIREF90_FASTA:-}" "${UNIREF90_FASTA_URL:-}"
    if [[ -z "$HOMOLOGY_COMMON_PREPROCESSING_CACHE" ]]; then
        require_local_or_url idmapping_selected "${IDMAPPING:-}" "${IDMAPPING_URL:-}"
        if [[ "$UNIPROT_SOURCE_SCOPE" != "trembl-only" ]]; then
            require_local_or_url Swiss-Prot-DAT "${UNIPROT_SPROT_SEQUENCES:-}" "${UNIPROT_SPROT_SEQUENCES_URL:-}"
        fi
        if [[ "$UNIPROT_SOURCE_SCOPE" != "sprot-only" ]]; then
            require_local_or_url TrEMBL-DAT "${UNIPROT_TREMBL_SEQUENCES:-}" "${UNIPROT_TREMBL_SEQUENCES_URL:-}"
        fi
        require_local_or_url GOA "${GOA:-}" "${GOA_URL:-}"
    fi
    require_local_or_url GO-OBO "${GO_OBO:-}" "$GO_OBO_URL"
    if [[ "$FIXTURE_MODE" != "1" && -n "${CLUSTER_ASSIGNMENTS:-}" ]]; then
        echo "CLUSTER_ASSIGNMENTS is fixture-only; set FIXTURE_MODE=1 explicitly" >&2
        exit 1
    fi
    if [[ "$FIXTURE_MODE" != "1" ]]; then
        [[ -f "$FROZEN_INPUT_MANIFEST" ]] || {
            echo "Production run requires FROZEN_INPUT_MANIFEST" >&2
            exit 1
        }
        [[ -n "$EXPECTED_MMSEQS_VERSION" ]] || {
            echo "Production run requires exact EXPECTED_MMSEQS_VERSION" >&2
            exit 1
        }
        [[ "$FRAMEWORK_REVISION" =~ ^[0-9a-f]{40}$ ]] || {
            echo "Production run requires FRAMEWORK_REVISION as exactly 40 lowercase hex characters" >&2
            exit 1
        }
        if [[ "$DIAGNOSTIC_PILOT" != "1" ]]; then
            [[ -f "$ATTRITION_POLICY" ]] || {
                echo "Production run requires reviewed ATTRITION_POLICY" >&2
                exit 1
            }
        fi
        required_hashes=(UNIREF90_FASTA_SHA256 IDMAPPING_SHA256 GOA_SHA256 GO_OBO_SHA256)
        [[ "$UNIPROT_SOURCE_SCOPE" == "trembl-only" ]] || required_hashes+=(UNIPROT_SPROT_SEQUENCES_SHA256)
        [[ "$UNIPROT_SOURCE_SCOPE" == "sprot-only" ]] || required_hashes+=(UNIPROT_TREMBL_SEQUENCES_SHA256)
        for hash_variable in "${required_hashes[@]}"
        do
            [[ -n "${!hash_variable:-}" ]] || {
                echo "Production run requires $hash_variable to pin the frozen input" >&2
                exit 1
            }
        done
    fi
    if [[ -z "${CLUSTER_ASSIGNMENTS:-}" ]]; then
        if [[ "$MMSEQS_BIN" == */* ]]; then
            [[ -x "$MMSEQS_BIN" ]] || { echo "MMseqs2 is not executable: $MMSEQS_BIN" >&2; exit 1; }
        else
            command -v "$MMSEQS_BIN" >/dev/null 2>&1 || {
                echo "MMseqs2 is unavailable; set MMSEQS_BIN to the compute-node executable" >&2
                exit 1
            }
        fi
    fi
fi

COMMAND=(
    "$PYTHON_BIN" -m homology_cluster_benchmark build
    --identity "$IDENTITY"
    --split-policy "$SPLIT_POLICY"
    --training-population "$TRAINING_POPULATION"
    --uniprot-source-scope "$UNIPROT_SOURCE_SCOPE"
    --mmseqs-bin "$MMSEQS_BIN"
    --output-dir "$OUTPUT_ROOT"
    --temp-dir "$TEMP_DIR"
    --threads "$THREADS"
    --run-id "$RUN_ID"
    --seed "$SEED"
    --min-count "$MIN_COUNT"
    --uniprot-release "$UNIPROT_RELEASE"
    --goa-release "$GOA_RELEASE"
    --ontology-release "$ONTOLOGY_RELEASE"
    --scratch-safety-multiplier "$SCRATCH_SAFETY_MULTIPLIER"
    --minimum-free-disk-gb "$MINIMUM_FREE_DISK_GB"
    --mmseqs-work-multiplier "$MMSEQS_WORK_MULTIPLIER"
    --publication-safety-multiplier "$PUBLICATION_SAFETY_MULTIPLIER"
    --excluded-sample-per-reason "$EXCLUDED_SAMPLE_PER_REASON"
)

if [[ -n "$EXPECTED_MMSEQS_VERSION" ]]; then
    COMMAND+=(--expected-mmseqs-version "$EXPECTED_MMSEQS_VERSION")
fi
if [[ -n "$FROZEN_INPUT_MANIFEST" ]]; then
    COMMAND+=(--frozen-input-manifest "$FROZEN_INPUT_MANIFEST")
fi
if [[ -n "$HOMOLOGY_COMMON_PREPROCESSING_CACHE" ]]; then
    COMMAND+=(--common-preprocessing-cache "$HOMOLOGY_COMMON_PREPROCESSING_CACHE")
fi
if [[ -n "$ATTRITION_POLICY" ]]; then
    COMMAND+=(--attrition-policy "$ATTRITION_POLICY")
fi
if [[ -n "$ATTRITION_OVERRIDE" ]]; then
    COMMAND+=(--attrition-override "$ATTRITION_OVERRIDE")
fi
if [[ -n "$FRAMEWORK_REVISION" ]]; then
    COMMAND+=(--framework-revision "$FRAMEWORK_REVISION")
fi
if [[ -n "$REQUESTED_SLOTS" ]]; then
    COMMAND+=(--requested-slots "$REQUESTED_SLOTS")
fi
if [[ -n "$ALLOCATED_SLOTS" ]]; then
    COMMAND+=(--allocated-slots "$ALLOCATED_SLOTS")
fi
if [[ -n "$PERSISTENT_RESULTS_ROOT" ]]; then
    COMMAND+=(--persistent-results-root "$PERSISTENT_RESULTS_ROOT")
fi

append_input uniref90-fasta "${UNIREF90_FASTA:-}" "${UNIREF90_FASTA_URL:-}" "${UNIREF90_FASTA_SHA256:-}"
append_input idmapping "${IDMAPPING:-}" "${IDMAPPING_URL:-}" "${IDMAPPING_SHA256:-}"
if [[ "$UNIPROT_SOURCE_SCOPE" != "trembl-only" ]]; then
    append_input uniprot-sprot-sequences "${UNIPROT_SPROT_SEQUENCES:-}" "${UNIPROT_SPROT_SEQUENCES_URL:-}" "${UNIPROT_SPROT_SEQUENCES_SHA256:-}"
fi
if [[ "$UNIPROT_SOURCE_SCOPE" != "sprot-only" ]]; then
    append_input uniprot-trembl-sequences "${UNIPROT_TREMBL_SEQUENCES:-}" "${UNIPROT_TREMBL_SEQUENCES_URL:-}" "${UNIPROT_TREMBL_SEQUENCES_SHA256:-}"
fi
append_input goa "${GOA:-}" "${GOA_URL:-}" "${GOA_SHA256:-}"
append_input go-obo "${GO_OBO:-}" "$GO_OBO_URL" "${GO_OBO_SHA256:-}"

if [[ -n "${CLUSTER_ASSIGNMENTS:-}" ]]; then
    [[ -f "$CLUSTER_ASSIGNMENTS" ]] || { echo "Cluster assignment fixture does not exist: $CLUSTER_ASSIGNMENTS" >&2; exit 1; }
    COMMAND+=(--cluster-assignments "$CLUSTER_ASSIGNMENTS")
fi
[[ "$FIXTURE_MODE" == "1" ]] && COMMAND+=(--fixture-mode)
[[ "$DIAGNOSTIC_PILOT" == "1" ]] && COMMAND+=(--diagnostic-pilot)
[[ "$NO_DOWNLOADS" == "1" ]] && COMMAND+=(--no-downloads)
[[ "$DRY_RUN" == "1" ]] && COMMAND+=(--dry-run)
[[ "$KEEP_TEMP" == "1" ]] && COMMAND+=(--keep-temp)
COMMAND+=("$@")

echo "Framework root: $FRAMEWORK_ROOT"
echo "Output root   : $OUTPUT_ROOT"
echo "Temporary root: $TEMP_DIR"
printf 'Command        : '
printf '%q ' "${COMMAND[@]}"
printf '\n'

export PYTHONPATH="$BUILDER_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
set +e
if [[ -n "$LOG_FILE" ]]; then
    mkdir -p "$(dirname "$LOG_FILE")"
    "${COMMAND[@]}" > >(tee "$LOG_FILE") 2>&1 &
else
    "${COMMAND[@]}" &
fi
ACTIVE_CHILD_PID=$!
while true; do
    wait "$ACTIVE_CHILD_PID"
    status=$?
    if ! kill -0 "$ACTIVE_CHILD_PID" 2>/dev/null; then
        break
    fi
done
if [[ -n "$SIGNAL_WATCHDOG_PID" ]]; then
    kill "$SIGNAL_WATCHDOG_PID" 2>/dev/null || true
    wait "$SIGNAL_WATCHDOG_PID" 2>/dev/null || true
    SIGNAL_WATCHDOG_PID=""
fi
ACTIVE_CHILD_PID=""
set -e
if [[ "$SIGNAL_STATUS" != "0" ]]; then
    exit "$SIGNAL_STATUS"
fi
exit "$status"
