#!/usr/bin/env bash

set -euo pipefail

LAUNCHER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRAMEWORK_ROOT="$(cd "$LAUNCHER_DIR/../.." && pwd)"
WORKER="$FRAMEWORK_ROOT/hpc_jobs/active/hpc_homology_cluster_benchmark.sh"
PYTHON_BIN="${PYTHON_BIN:-python3}"
DRY_RUN="${DRY_RUN:-0}"
REQUESTED_SLOTS=8
SPLIT_POLICY="${SPLIT_POLICY:-sequence-balanced}"
TRAINING_POPULATION="${TRAINING_POPULATION:-annotated-only}"
SEED="${SEED:-0}"
MIN_COUNT="${MIN_COUNT:-50}"
NO_DOWNLOADS=1
FIXTURE_MODE=0
UNIPROT_RELEASE="${UNIPROT_RELEASE:-2026_02}"
GOA_RELEASE="${GOA_RELEASE:-234}"
ONTOLOGY_RELEASE="${ONTOLOGY_RELEASE:-releases/2026-06-15}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)-$$-${RANDOM}}"

launcher_error() {
    echo "$*" >&2
    exit 2
}

require_value() {
    local name="$1"
    [[ -n "${!name:-}" ]] || launcher_error "Required launcher variable is missing: $name"
}

require_file() {
    local name="$1"
    require_value "$name"
    [[ -f "${!name}" ]] || launcher_error "$name is not a local file: ${!name}"
}

file_sha256() {
    "$PYTHON_BIN" -c \
        'import hashlib,sys; h=hashlib.sha256(); f=open(sys.argv[1], "rb"); [h.update(x) for x in iter(lambda:f.read(1048576), b"")]; print(h.hexdigest())' \
        "$1"
}

require_sha256() {
    local name="$1"
    require_value "$name"
    [[ "${!name}" =~ ^[0-9a-f]{64}$ ]] || launcher_error "$name must be 64 lowercase hexadecimal characters"
}

require_pinned_file() {
    local path_name="$1"
    local hash_name="$2"
    require_file "$path_name"
    require_sha256 "$hash_name"
    local observed
    observed="$(file_sha256 "${!path_name}")"
    [[ "$observed" == "${!hash_name}" ]] || launcher_error \
        "$path_name does not match reviewed $hash_name"
}

validate_common_contract() {
    case "$DRY_RUN" in 0|1) ;; *) launcher_error "DRY_RUN must be 0 or 1" ;; esac
    require_value UNIPROT_SOURCE_SCOPE
    case "$UNIPROT_SOURCE_SCOPE" in
        sprot-only|trembl-only|sprot-and-trembl) ;;
        *) launcher_error "UNIPROT_SOURCE_SCOPE must be explicit" ;;
    esac
    require_value FRAMEWORK_REVISION
    [[ "$FRAMEWORK_REVISION" =~ ^[0-9a-f]{40}$ ]] || launcher_error \
        "FRAMEWORK_REVISION must be exactly 40 lowercase hexadecimal characters"
    require_file FROZEN_INPUT_MANIFEST
    EXPECTED_FROZEN_INPUT_MANIFEST_SHA256="$(file_sha256 "$FROZEN_INPUT_MANIFEST")"
    export EXPECTED_FROZEN_INPUT_MANIFEST_SHA256
    require_pinned_file UNIREF90_FASTA UNIREF90_FASTA_SHA256
    require_pinned_file IDMAPPING IDMAPPING_SHA256
    require_pinned_file GOA GOA_SHA256
    require_pinned_file GO_OBO GO_OBO_SHA256
    require_value RESULTS_ROOT
    [[ "$RESULTS_ROOT" == /* ]] || launcher_error "RESULTS_ROOT must be absolute"
    require_value EXPECTED_MMSEQS_VERSION
    [[ "$EXPECTED_MMSEQS_VERSION" != *'<'* && "$EXPECTED_MMSEQS_VERSION" != *'>'* ]] || \
        launcher_error "EXPECTED_MMSEQS_VERSION still contains a placeholder"
    require_value MMSEQS_BIN
    case "$SPLIT_POLICY" in cluster-count-random|sequence-balanced) ;; *) launcher_error "Unsupported SPLIT_POLICY" ;; esac
    [[ "$TRAINING_POPULATION" == "annotated-only" ]] || launcher_error \
        "TRAINING_POPULATION must remain annotated-only"
    [[ "$SEED" =~ ^[0-9]+$ ]] || launcher_error "SEED must be a non-negative integer"
    [[ "$MIN_COUNT" =~ ^[1-9][0-9]*$ ]] || launcher_error "MIN_COUNT must be a positive integer"
    [[ "$NO_DOWNLOADS" == "1" && "$REQUESTED_SLOTS" == "8" && "$FIXTURE_MODE" == "0" ]] || \
        launcher_error "Production launcher constants were altered"
    require_value UNIPROT_RELEASE
    require_value GOA_RELEASE
    require_value ONTOLOGY_RELEASE
    if [[ "$UNIPROT_SOURCE_SCOPE" != "trembl-only" ]]; then
        require_pinned_file UNIPROT_SPROT_SEQUENCES UNIPROT_SPROT_SEQUENCES_SHA256
    elif [[ -n "${UNIPROT_SPROT_SEQUENCES:-}${UNIPROT_SPROT_SEQUENCES_SHA256:-}" ]]; then
        launcher_error "trembl-only forbids Swiss-Prot path and hash variables"
    fi
    if [[ "$UNIPROT_SOURCE_SCOPE" != "sprot-only" ]]; then
        require_pinned_file UNIPROT_TREMBL_SEQUENCES UNIPROT_TREMBL_SEQUENCES_SHA256
    elif [[ -n "${UNIPROT_TREMBL_SEQUENCES:-}${UNIPROT_TREMBL_SEQUENCES_SHA256:-}" ]]; then
        launcher_error "sprot-only forbids TrEMBL path and hash variables"
    fi
    [[ "$RUN_ID" =~ ^[A-Za-z0-9._-]+$ && "$RUN_ID" =~ [A-Za-z0-9] ]] || \
        launcher_error "RUN_ID must be path-safe and contain an alphanumeric"
}

export_names() {
    local names=(
        UNIPROT_SOURCE_SCOPE FRAMEWORK_REVISION FROZEN_INPUT_MANIFEST
        UNIREF90_FASTA UNIREF90_FASTA_SHA256 IDMAPPING IDMAPPING_SHA256
        GOA GOA_SHA256 GO_OBO GO_OBO_SHA256 RESULTS_ROOT EXPECTED_MMSEQS_VERSION
        MMSEQS_BIN SPLIT_POLICY TRAINING_POPULATION SEED MIN_COUNT NO_DOWNLOADS
        RUN_ID REQUESTED_SLOTS FIXTURE_MODE DIAGNOSTIC_PILOT
        UNIPROT_RELEASE GOA_RELEASE ONTOLOGY_RELEASE
        EXPECTED_FROZEN_INPUT_MANIFEST_SHA256
    )
    [[ "$UNIPROT_SOURCE_SCOPE" == "trembl-only" ]] || names+=(
        UNIPROT_SPROT_SEQUENCES UNIPROT_SPROT_SEQUENCES_SHA256
    )
    [[ "$UNIPROT_SOURCE_SCOPE" == "sprot-only" ]] || names+=(
        UNIPROT_TREMBL_SEQUENCES UNIPROT_TREMBL_SEQUENCES_SHA256
    )
    [[ -z "${ATTRITION_POLICY:-}" ]] || names+=(ATTRITION_POLICY)
    [[ -z "${PILOT_APPROVAL:-}" ]] || names+=(PILOT_APPROVAL)
    [[ -z "${PILOT_RUN_DIR:-}" ]] || names+=(PILOT_RUN_DIR)
    [[ -z "${PILOT_COMPLETION_MARKER:-}" ]] || names+=(PILOT_COMPLETION_MARKER)
    [[ -z "${PILOT_ATTRITION_REPORT:-}" ]] || names+=(PILOT_ATTRITION_REPORT)
    [[ -z "${PILOT_TASK_CONTEXT:-}" ]] || names+=(PILOT_TASK_CONTEXT)
    [[ -z "${PILOT_MEASUREMENT_EVIDENCE:-}" ]] || names+=(PILOT_MEASUREMENT_EVIDENCE)
    [[ -z "${EXPECTED_ATTRITION_POLICY_SHA256:-}" ]] || names+=(EXPECTED_ATTRITION_POLICY_SHA256)
    [[ -z "${EXPECTED_PILOT_APPROVAL_SHA256:-}" ]] || names+=(EXPECTED_PILOT_APPROVAL_SHA256)
    [[ -z "${EXPECTED_PILOT_COMPLETION_MARKER_SHA256:-}" ]] || names+=(EXPECTED_PILOT_COMPLETION_MARKER_SHA256)
    [[ -z "${EXPECTED_PILOT_ATTRITION_REPORT_SHA256:-}" ]] || names+=(EXPECTED_PILOT_ATTRITION_REPORT_SHA256)
    [[ -z "${EXPECTED_PILOT_TASK_CONTEXT_SHA256:-}" ]] || names+=(EXPECTED_PILOT_TASK_CONTEXT_SHA256)
    [[ -z "${EXPECTED_PILOT_MEASUREMENT_EVIDENCE_SHA256:-}" ]] || names+=(EXPECTED_PILOT_MEASUREMENT_EVIDENCE_SHA256)
    local joined
    joined="$(IFS=,; echo "${names[*]}")"
    printf '%s' "$joined"
}

print_preview() {
    local mode="$1"
    local task_range="$2"
    shift 2
    echo "Launcher mode       : $mode"
    echo "Task range          : $task_range"
    echo "UniProt source scope: $UNIPROT_SOURCE_SCOPE"
    echo "Framework commit    : $FRAMEWORK_REVISION"
    echo "Frozen manifest     : $FROZEN_INPUT_MANIFEST"
    echo "Split policy        : $SPLIT_POLICY"
    echo "Training population : $TRAINING_POPULATION"
    echo "Result root         : $RESULTS_ROOT"
    echo "Requested PE        : smp $REQUESTED_SLOTS"
    echo "Attrition policy    : ${ATTRITION_POLICY:-none-diagnostic-pilot}"
    echo "Pilot approval      : ${PILOT_APPROVAL:-not-applicable}"
    echo "Completion marker   : ${PILOT_COMPLETION_MARKER:-not-applicable}"
    echo "Attrition report    : ${PILOT_ATTRITION_REPORT:-not-applicable}"
    echo "Exported variables  : $(export_names)"
    local -a export_list
    local export_name
    IFS=',' read -r -a export_list <<< "$(export_names)"
    echo "Exported values:"
    for export_name in "${export_list[@]}"; do
        printf '  %s=%q\n' "$export_name" "${!export_name}"
    done
    printf 'qsub command        : '
    printf '%q ' "$@"
    printf '\n'
}

launch_array() {
    local mode="$1"
    local task_range="$2"
    local exports
    exports="$(export_names)"
    local command=(qsub -t "$task_range" -pe smp "$REQUESTED_SLOTS" -v "$exports" "$WORKER")
    print_preview "$mode" "$task_range" "${command[@]}"
    if [[ "$DRY_RUN" == "1" ]]; then
        echo "Dry run only: Grid Engine was not contacted."
        return 0
    fi
    "${command[@]}"
}
