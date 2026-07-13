#!/usr/bin/env bash
# Build the contemporary temporal benchmark from frozen UniProt, GOA and GO inputs.
# Existing files are reused; missing inputs are downloaded into the selected DB_ROOT.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRAMEWORK_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILDER_ROOT="$FRAMEWORK_ROOT/benchmark_builders/contemporary_cafa"
FILTER_DAT="$SCRIPT_DIR/filter_uniprot_dat.py"
EXTRACT_MEMBER="$SCRIPT_DIR/extract_tar_member.py"
TARGET_TAXA="$BUILDER_ROOT/src/cafa_benchmark_builder/resources/cafa3_target_taxa.txt"

PROFILE="${PROFILE:-contemporary-cafa3-style}"
DB_ROOT="${DB_ROOT:-$HOME/protein_databases}"
RUN_ROOT="${RUN_ROOT:-$PWD/contemporary_cafa_${PROFILE}_$(date +%Y%m%d_%H%M%S)}"
WORK_DIR="${WORK_DIR:-$RUN_ROOT/work}"
OUTPUT_DIR="$RUN_ROOT/outputs"
REPORT_DIR="$RUN_ROOT/reports"
LOG_DIR="$RUN_ROOT/logs"
PYTHON_BIN="${PYTHON_BIN:-python}"
ALLOW_SPROT_ONLY="${ALLOW_SPROT_ONLY:-0}"
REMOVE_ARCHIVES_AFTER_EXTRACT="${REMOVE_ARCHIVES_AFTER_EXTRACT:-0}"
PIGZ_THREADS="${PIGZ_THREADS:-1}"
MIN_COUNT="${MIN_COUNT:-50}"
SPLIT="${SPLIT:-0.9}"
SEED="${SEED:-0}"
NO_STRICT_QC="${NO_STRICT_QC:-0}"
SKIP_INPUT_CHECKSUMS="${SKIP_INPUT_CHECKSUMS:-0}"
ALLOW_DEPENDENCY_MISMATCH="${ALLOW_DEPENDENCY_MISMATCH:-0}"
T1_ENDPOINT_POLICY="${T1_ENDPOINT_POLICY:-snapshot-membership}"
T1_BACKFILL_POLICY="${T1_BACKFILL_POLICY:-allow}"

T0_UNIPROT_DIR="$DB_ROOT/uniprot/release_2025_01"
T1_UNIPROT_DIR="$DB_ROOT/uniprot/release_2026_02"
T0_GOA_DIR="$DB_ROOT/goa/release_2025_01"
T1_GOA_DIR="$DB_ROOT/goa/release_2026_02"
ONTOLOGY_T0_BENCHMARK_DIR="$DB_ROOT/ontology/release_2025-02-06"
ONTOLOGY_T0_SOURCE_DIR="$DB_ROOT/ontology/release_2025-03-16"
ONTOLOGY_T1_SOURCE_DIR="$DB_ROOT/ontology/release_2026-06-19"

T0_SPROT_DAT="$T0_UNIPROT_DIR/uniprot_sprot.dat.gz"
T0_SPROT_ARCHIVE="$T0_UNIPROT_DIR/uniprot_sprot-only2025_01.tar.gz"
T0_FULL_ARCHIVE="$T0_UNIPROT_DIR/knowledgebase2025_01.tar.gz"
T0_TREMBL_FILTERED="$T0_UNIPROT_DIR/uniprot_trembl_cafa3_targets.dat.gz"
T1_SPROT_DAT="$T1_UNIPROT_DIR/uniprot_sprot.dat.gz"
T1_TREMBL_FULL="$T1_UNIPROT_DIR/uniprot_trembl.dat.gz"
T1_TREMBL_FILTERED="$T1_UNIPROT_DIR/uniprot_trembl_cafa3_targets.dat.gz"
T0_GOA="$T0_GOA_DIR/goa_uniprot_all.gaf.225.gz"
T1_GOA="$T1_GOA_DIR/goa_uniprot_all.gaf.234.gz"
T0_BENCHMARK_OBO="$ONTOLOGY_T0_BENCHMARK_DIR/go-basic.obo"
T0_SOURCE_OBO="$ONTOLOGY_T0_SOURCE_DIR/go-basic.obo"
T1_SOURCE_OBO="$ONTOLOGY_T1_SOURCE_DIR/go-basic.obo"

UNIPROT_2025_BASE="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2025_01/knowledgebase"
UNIPROT_CURRENT_BASE="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"
UNIPROT_2026_ARCHIVE_BASE="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2026_02/knowledgebase"
GOA_OLD_BASE="https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT"
GOA_CURRENT_BASE="https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT"

mkdir -p "$WORK_DIR" "$OUTPUT_DIR" "$REPORT_DIR" "$LOG_DIR"
ACQUISITION_LOG="$LOG_DIR/input_acquisition.tsv"
printf 'input\taction\tsource\tdestination\n' > "$ACQUISITION_LOG"

record_acquisition() {
    printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" >> "$ACQUISITION_LOG"
}

maybe_remove_staged_source() {
    local path="$1"
    if [[ "$REMOVE_ARCHIVES_AFTER_EXTRACT" == "1" && "$path" == "$DB_ROOT/"* ]]; then
        rm -f "$path"
    fi
}

require_file() {
    if [[ ! -f "$1" ]]; then
        echo "Missing required input: $1" >&2
        exit 1
    fi
}

compress_stream() {
    if command -v pigz >/dev/null 2>&1; then
        pigz -p "$PIGZ_THREADS" -c
    else
        gzip -c
    fi
}

decompress_stream() {
    if command -v pigz >/dev/null 2>&1; then
        pigz -p "$PIGZ_THREADS" -dc
    else
        gzip -dc
    fi
}

download_file() {
    local label="$1"
    local destination="$2"
    shift 2
    if [[ -f "$destination" ]]; then
        record_acquisition "$label" reused "$destination" "$destination"
        return 0
    fi
    mkdir -p "$(dirname "$destination")"
    local temporary="${destination}.part"
    local url
    for url in "$@"; do
        echo "Downloading $label from $url"
        rm -f "$temporary"
        if wget --progress=dot:giga -O "$temporary" "$url"; then
            mv "$temporary" "$destination"
            record_acquisition "$label" downloaded "$url" "$destination"
            return 0
        fi
    done
    rm -f "$temporary"
    echo "Could not download $label from any configured URL" >&2
    return 1
}

extract_gz_member_from_archive() {
    local label="$1"
    local archive="$2"
    local suffix="$3"
    local destination="$4"
    mkdir -p "$(dirname "$destination")"
    local temporary="${destination}.part"
    echo "Extracting $suffix from $archive"
    rm -f "$temporary"
    if "$PYTHON_BIN" "$EXTRACT_MEMBER" --archive "$archive" --suffix "$suffix" > "$temporary"; then
        mv "$temporary" "$destination"
        record_acquisition "$label" extracted "$archive#$suffix" "$destination"
        maybe_remove_staged_source "$archive"
    else
        rm -f "$temporary"
        return 1
    fi
}

stream_gz_member_from_url() {
    local label="$1"
    local url="$2"
    local suffix="$3"
    local destination="$4"
    mkdir -p "$(dirname "$destination")"
    local temporary="${destination}.part"
    echo "Streaming $label from $url"
    rm -f "$temporary"
    if wget --progress=dot:giga -O - "$url" \
        | "$PYTHON_BIN" "$EXTRACT_MEMBER" --suffix "$suffix" > "$temporary"; then
        mv "$temporary" "$destination"
        record_acquisition "$label" streamed "$url#$suffix" "$destination"
    else
        rm -f "$temporary"
        return 1
    fi
}

filter_dat_gz() {
    local label="$1"
    local source="$2"
    local destination="$3"
    mkdir -p "$(dirname "$destination")"
    local temporary="${destination}.part"
    echo "Filtering $label to CAFA3 target taxa"
    rm -f "$temporary"
    if decompress_stream < "$source" \
        | "$PYTHON_BIN" "$FILTER_DAT" --taxa-file "$TARGET_TAXA" \
        | compress_stream > "$temporary"; then
        mv "$temporary" "$destination"
        record_acquisition "$label" filtered "$source" "$destination"
        maybe_remove_staged_source "$source"
    else
        rm -f "$temporary"
        return 1
    fi
}

filter_tar_member() {
    local label="$1"
    local archive="$2"
    local suffix="$3"
    local destination="$4"
    mkdir -p "$(dirname "$destination")"
    local temporary="${destination}.part"
    echo "Streaming $suffix from $archive and filtering to CAFA3 target taxa"
    rm -f "$temporary"
    if "$PYTHON_BIN" "$EXTRACT_MEMBER" --archive "$archive" --suffix "$suffix" \
        | decompress_stream \
        | "$PYTHON_BIN" "$FILTER_DAT" --taxa-file "$TARGET_TAXA" \
        | compress_stream > "$temporary"; then
        mv "$temporary" "$destination"
        record_acquisition "$label" filtered "$archive#$suffix" "$destination"
        maybe_remove_staged_source "$archive"
    else
        rm -f "$temporary"
        return 1
    fi
}

filter_tar_member_url() {
    local label="$1"
    local url="$2"
    local suffix="$3"
    local destination="$4"
    mkdir -p "$(dirname "$destination")"
    local temporary="${destination}.part"
    echo "Streaming $label from $url and retaining only CAFA3 target taxa"
    rm -f "$temporary"
    if wget --progress=dot:giga -O - "$url" \
        | "$PYTHON_BIN" "$EXTRACT_MEMBER" --suffix "$suffix" \
        | decompress_stream \
        | "$PYTHON_BIN" "$FILTER_DAT" --taxa-file "$TARGET_TAXA" \
        | compress_stream > "$temporary"; then
        mv "$temporary" "$destination"
        record_acquisition "$label" streamed-filtered "$url#$suffix" "$destination"
    else
        rm -f "$temporary"
        return 1
    fi
}

filter_dat_url() {
    local label="$1"
    local url="$2"
    local destination="$3"
    mkdir -p "$(dirname "$destination")"
    local temporary="${destination}.part"
    echo "Streaming $label from $url and retaining only CAFA3 target taxa"
    rm -f "$temporary"
    if wget --progress=dot:giga -O - "$url" \
        | decompress_stream \
        | "$PYTHON_BIN" "$FILTER_DAT" --taxa-file "$TARGET_TAXA" \
        | compress_stream > "$temporary"; then
        mv "$temporary" "$destination"
        record_acquisition "$label" streamed-filtered "$url" "$destination"
    else
        rm -f "$temporary"
        return 1
    fi
}

current_uniprot_is_2026_02() {
    local marker="$WORK_DIR/uniprot_current_reldate.txt"
    if ! wget -q -O "$marker" "$UNIPROT_CURRENT_BASE/reldate.txt"; then
        return 1
    fi
    grep -q '2026_02' "$marker"
}

ensure_ontologies() {
    download_file go_t0_benchmark "$T0_BENCHMARK_OBO" \
        "https://release.geneontology.org/2025-02-06/ontology/go-basic.obo"
    download_file go_t0_source "$T0_SOURCE_OBO" \
        "https://release.geneontology.org/2025-03-16/ontology/go-basic.obo"
    download_file go_t1_source "$T1_SOURCE_OBO" \
        "https://release.geneontology.org/2026-06-19/ontology/go-basic.obo"
    grep -q '^data-version: releases/2025-02-06' "$T0_BENCHMARK_OBO"
    grep -q '^data-version: releases/2025-03-16' "$T0_SOURCE_OBO"
    grep -q '^data-version: releases/2026-06-15' "$T1_SOURCE_OBO"
}

ensure_goa() {
    download_file goa_225 "$T0_GOA" \
        "$GOA_OLD_BASE/goa_uniprot_all.gaf.225.gz"
    download_file goa_234 "$T1_GOA" \
        "$GOA_OLD_BASE/goa_uniprot_all.gaf.234.gz" \
        "$GOA_CURRENT_BASE/goa_uniprot_all.gaf.gz"
    if ! zgrep -m 1 '^!date-generated: 2025-03-08' "$T0_GOA" >/dev/null \
        || ! zgrep -m 1 '^!go-version:.*2025-03-07' "$T0_GOA" >/dev/null; then
        echo "GOA 225 header does not match the frozen t0 snapshot" >&2
        exit 1
    fi
    if ! zgrep -m 1 '^!date-generated: 2026-06-17' "$T1_GOA" >/dev/null \
        || ! zgrep -m 1 '^!go-version:.*2026-06-15' "$T1_GOA" >/dev/null; then
        echo "GOA 234 header does not match the frozen t1 snapshot" >&2
        exit 1
    fi
}

ensure_default_uniprot_inputs() {
    if [[ ! -f "$T0_SPROT_DAT" ]]; then
        if [[ -f "$T0_SPROT_ARCHIVE" ]]; then
            extract_gz_member_from_archive uniprot_t0_sprot "$T0_SPROT_ARCHIVE" \
                uniprot_sprot.dat.gz "$T0_SPROT_DAT"
        else
            stream_gz_member_from_url uniprot_t0_sprot \
                "$UNIPROT_2025_BASE/uniprot_sprot-only2025_01.tar.gz" \
                uniprot_sprot.dat.gz "$T0_SPROT_DAT"
        fi
    fi

    if [[ "$ALLOW_SPROT_ONLY" != "1" && ! -f "$T0_TREMBL_FILTERED" ]]; then
        if [[ -n "${T0_TREMBL_DAT_SOURCE:-}" ]]; then
            filter_dat_gz uniprot_t0_trembl "$T0_TREMBL_DAT_SOURCE" "$T0_TREMBL_FILTERED"
        elif [[ -n "${T0_TREMBL_ARCHIVE_SOURCE:-}" ]]; then
            filter_tar_member uniprot_t0_trembl "$T0_TREMBL_ARCHIVE_SOURCE" \
                uniprot_trembl.dat.gz "$T0_TREMBL_FILTERED"
        elif [[ -f "$T0_UNIPROT_DIR/uniprot_trembl.dat.gz" ]]; then
            filter_dat_gz uniprot_t0_trembl "$T0_UNIPROT_DIR/uniprot_trembl.dat.gz" \
                "$T0_TREMBL_FILTERED"
        elif [[ -f "$T0_FULL_ARCHIVE" ]]; then
            filter_tar_member uniprot_t0_trembl "$T0_FULL_ARCHIVE" \
                uniprot_trembl.dat.gz "$T0_TREMBL_FILTERED"
        else
            filter_tar_member_url uniprot_t0_trembl \
                "$UNIPROT_2025_BASE/knowledgebase2025_01.tar.gz" \
                uniprot_trembl.dat.gz "$T0_TREMBL_FILTERED"
        fi
    fi

    local current_release=0
    if current_uniprot_is_2026_02; then
        current_release=1
    fi
    if [[ ! -f "$T1_SPROT_DAT" ]]; then
        if [[ "$current_release" == "1" ]]; then
            download_file uniprot_t1_sprot "$T1_SPROT_DAT" \
                "$UNIPROT_CURRENT_BASE/uniprot_sprot.dat.gz"
        else
            stream_gz_member_from_url uniprot_t1_sprot \
                "$UNIPROT_2026_ARCHIVE_BASE/uniprot_sprot-only2026_02.tar.gz" \
                uniprot_sprot.dat.gz "$T1_SPROT_DAT"
        fi
    fi
    if [[ "$ALLOW_SPROT_ONLY" != "1" && ! -f "$T1_TREMBL_FILTERED" ]]; then
        if [[ -n "${T1_TREMBL_DAT_SOURCE:-}" ]]; then
            filter_dat_gz uniprot_t1_trembl "$T1_TREMBL_DAT_SOURCE" "$T1_TREMBL_FILTERED"
        elif [[ -f "$T1_TREMBL_FULL" ]]; then
            filter_dat_gz uniprot_t1_trembl "$T1_TREMBL_FULL" "$T1_TREMBL_FILTERED"
        elif [[ "$current_release" == "1" ]]; then
            filter_dat_url uniprot_t1_trembl "$UNIPROT_CURRENT_BASE/uniprot_trembl.dat.gz" \
                "$T1_TREMBL_FILTERED"
        else
            filter_tar_member_url uniprot_t1_trembl \
                "$UNIPROT_2026_ARCHIVE_BASE/knowledgebase2026_02.tar.gz" \
                uniprot_trembl.dat.gz "$T1_TREMBL_FILTERED"
        fi
    fi
}

case "$PROFILE" in
    contemporary-cafa3-style|supervisor) ;;
    *)
        echo "The 2025->2026 runner supports contemporary-cafa3-style or supervisor, got: $PROFILE" >&2
        exit 1
        ;;
esac

case "$T1_ENDPOINT_POLICY" in
    assigned-date-proxy|snapshot-membership) ;;
    *)
        echo "Unknown T1_ENDPOINT_POLICY: $T1_ENDPOINT_POLICY" >&2
        exit 1
        ;;
esac
case "$T1_BACKFILL_POLICY" in
    exclude-pre-t0|allow) ;;
    *)
        echo "Unknown T1_BACKFILL_POLICY: $T1_BACKFILL_POLICY" >&2
        exit 1
        ;;
esac

for command in "$PYTHON_BIN" gzip grep tee wget zgrep; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "Missing required command: $command" >&2
        exit 1
    fi
done
require_file "$FILTER_DAT"
require_file "$EXTRACT_MEMBER"
require_file "$TARGET_TAXA"

"$PYTHON_BIN" -c 'import numpy, pandas'

ensure_ontologies
ensure_goa

T0_INPUTS=()
T1_INPUTS=()
if [[ -n "${UNIPROT_T0_INPUTS:-}" || -n "${UNIPROT_T1_INPUTS:-}" ]]; then
    if [[ -z "${UNIPROT_T0_INPUTS:-}" || -z "${UNIPROT_T1_INPUTS:-}" ]]; then
        echo "UNIPROT_T0_INPUTS and UNIPROT_T1_INPUTS must be supplied together" >&2
        exit 1
    fi
    IFS=':' read -r -a T0_INPUTS <<< "$UNIPROT_T0_INPUTS"
    IFS=':' read -r -a T1_INPUTS <<< "$UNIPROT_T1_INPUTS"
    for path in "${T0_INPUTS[@]}" "${T1_INPUTS[@]}"; do
        require_file "$path"
        record_acquisition uniprot_custom reused "$path" "$path"
    done
else
    ensure_default_uniprot_inputs
    T0_INPUTS+=("$T0_SPROT_DAT")
    T1_INPUTS+=("$T1_SPROT_DAT")
    if [[ "$ALLOW_SPROT_ONLY" != "1" ]]; then
        T0_INPUTS+=("$T0_TREMBL_FILTERED")
        T1_INPUTS+=("$T1_TREMBL_FILTERED")
    fi
fi

for path in "${T0_INPUTS[@]}" "${T1_INPUTS[@]}"; do
    require_file "$path"
done

COMMAND=(
    "$PYTHON_BIN" -m cafa_benchmark_builder
    --source-mode snapshots
    --profile "$PROFILE"
    --goa-t0 "$T0_GOA"
    --goa-t1 "$T1_GOA"
    --go-obo "$T0_BENCHMARK_OBO"
    --go-obo-t0 "$T0_SOURCE_OBO"
    --go-obo-t1 "$T1_SOURCE_OBO"
    --output-dir "$OUTPUT_DIR"
    --report-dir "$REPORT_DIR"
    --min-count "$MIN_COUNT"
    --split "$SPLIT"
    --seed "$SEED"
    --t1-endpoint-policy "$T1_ENDPOINT_POLICY"
)
if [[ "$T1_BACKFILL_POLICY" == "allow" ]]; then
    COMMAND+=(--allow-t1-backfill)
else
    COMMAND+=(--exclude-t1-backfill)
fi
for path in "${T0_INPUTS[@]}"; do
    COMMAND+=(--uniprot-t0 "$path")
done
for path in "${T1_INPUTS[@]}"; do
    COMMAND+=(--uniprot-t1 "$path")
done
if [[ "$ALLOW_SPROT_ONLY" == "1" ]]; then
    COMMAND+=(--target-reviewed-only)
fi
if [[ "$NO_STRICT_QC" == "1" ]]; then
    COMMAND+=(--no-strict-qc)
fi
if [[ "$SKIP_INPUT_CHECKSUMS" == "1" ]]; then
    COMMAND+=(--skip-input-checksums)
fi

printf '%q ' "${COMMAND[@]}" > "$LOG_DIR/command.txt"
printf '\n' >> "$LOG_DIR/command.txt"
{
    echo "profile=$PROFILE"
    echo "db_root=$DB_ROOT"
    echo "run_root=$RUN_ROOT"
    echo "python=$($PYTHON_BIN --version 2>&1)"
    echo "hostname=$(hostname)"
    echo "min_count=$MIN_COUNT"
    echo "split=$SPLIT"
    echo "seed=$SEED"
    echo "t1_endpoint_policy=$T1_ENDPOINT_POLICY"
    echo "t1_backfill_policy=$T1_BACKFILL_POLICY"
    echo "strict_qc=$((1 - NO_STRICT_QC))"
} > "$LOG_DIR/environment.txt"

echo "Running contemporary CAFA benchmark builder"
echo "Profile : $PROFILE"
echo "Endpoint: $T1_ENDPOINT_POLICY"
echo "Backfill: $T1_BACKFILL_POLICY"
echo "Outputs : $OUTPUT_DIR"
echo "Reports : $REPORT_DIR"
echo

export PYTHONPATH="$BUILDER_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export CAFA_BUILDER_USE_PIGZ="${CAFA_BUILDER_USE_PIGZ:-1}"
"${COMMAND[@]}" 2>&1 | tee "$LOG_DIR/builder.log"

echo
echo "Benchmark build completed successfully: $RUN_ROOT"
