#!/usr/bin/env bash
# Build the 2025 -> 2026 temporal benchmark from explicit database snapshots.
# This script is machine-agnostic; the SGE wrapper only stages data and calls it.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRAMEWORK_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILDER_ROOT="$FRAMEWORK_ROOT/benchmark_builders/contemporary_cafa"

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

T0_UNIPROT_DIR="$DB_ROOT/uniprot/release_2025_01"
T1_UNIPROT_DIR="$DB_ROOT/uniprot/release_2026_02"
T0_GOA="$DB_ROOT/goa/release_2025_01/goa_uniprot_all.gaf.225.gz"
T1_GOA="$DB_ROOT/goa/release_2026_02/goa_uniprot_all.gaf.234.gz"
T0_OBO="$DB_ROOT/ontology/release_2025-03-07/go-basic.obo"
T1_OBO="$DB_ROOT/ontology/release_2026-06-15/go-basic.obo"

require_file() {
    if [[ ! -f "$1" ]]; then
        echo "Missing required input: $1" >&2
        exit 1
    fi
}

extract_dat_archive() {
    local archive="$1"
    local label="$2"
    local database="$3"
    local destination="$WORK_DIR/uniprot_${label}"
    mkdir -p "$destination"
    echo "Extracting $(basename "$archive") into $destination" >&2
    local member_list="$destination/archive_members.txt"
    tar -tzf "$archive" > "$member_list"
    local member
    member="$(awk -v expected="uniprot_${database}.dat.gz" \
        '$0 == expected || $0 ~ ("/" expected "$") {print; exit}' "$member_list")"
    if [[ -z "$member" ]]; then
        echo "No uniprot_${database}.dat.gz member found in $archive" >&2
        exit 1
    fi
    tar -xzf "$archive" -C "$destination" "$member"
    rm -f "$member_list"
    local dat_file
    dat_file="$(find "$destination" -type f -name '*.dat.gz' -print -quit)"
    if [[ -z "$dat_file" ]]; then
        echo "No .dat.gz file found in $archive" >&2
        exit 1
    fi
    if [[ "$REMOVE_ARCHIVES_AFTER_EXTRACT" == "1" ]]; then
        rm -f "$archive"
    fi
    printf '%s\n' "$dat_file"
}

mkdir -p "$WORK_DIR" "$OUTPUT_DIR" "$REPORT_DIR" "$LOG_DIR"

for command in "$PYTHON_BIN" tar find tee zgrep; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "Missing required command: $command" >&2
        exit 1
    fi
done

require_file "$T0_GOA"
require_file "$T1_GOA"
require_file "$T0_OBO"
require_file "$T1_OBO"

if ! zgrep -m 1 '^!go-version:.*2025-03-07' "$T0_GOA" >/dev/null; then
    echo "GOA 225 does not declare GO version 2025-03-07" >&2
    exit 1
fi
if ! zgrep -m 1 '^!go-version:.*2026-06-15' "$T1_GOA" >/dev/null; then
    echo "GOA 234 does not declare GO version 2026-06-15" >&2
    exit 1
fi

"$PYTHON_BIN" -c '
import importlib.metadata
expected = {"numpy": "2.0.2", "pandas": "2.3.3"}
actual = {name: importlib.metadata.version(name) for name in expected}
if actual != expected:
    raise SystemExit(f"Builder dependency mismatch: expected {expected}, found {actual}")
'

T0_INPUTS=()
T1_INPUTS=()

if [[ -n "${UNIPROT_T0_INPUTS:-}" ]]; then
    IFS=':' read -r -a T0_INPUTS <<< "$UNIPROT_T0_INPUTS"
else
    T0_SPROT="$T0_UNIPROT_DIR/uniprot_sprot.dat.gz"
    if [[ ! -f "$T0_SPROT" ]]; then
        T0_SPROT="$(extract_dat_archive \
            "$T0_UNIPROT_DIR/uniprot_sprot-only2025_01.tar.gz" sprot_t0 sprot)"
    fi
    T0_INPUTS+=("$T0_SPROT")

    T0_TREMBL="$T0_UNIPROT_DIR/uniprot_trembl.dat.gz"
    if [[ ! -f "$T0_TREMBL" && -f "$T0_UNIPROT_DIR/uniprot_trembl-only2025_01.tar.gz" ]]; then
        T0_TREMBL="$(extract_dat_archive \
            "$T0_UNIPROT_DIR/uniprot_trembl-only2025_01.tar.gz" trembl_t0 trembl)"
    fi
    if [[ -f "$T0_TREMBL" ]]; then
        T0_INPUTS+=("$T0_TREMBL")
    elif [[ "$ALLOW_SPROT_ONLY" != "1" ]]; then
        echo "Missing t0 TrEMBL input. Set ALLOW_SPROT_ONLY=1 only for a diagnostic Swiss-Prot-only build." >&2
        exit 1
    fi
fi

if [[ -n "${UNIPROT_T1_INPUTS:-}" ]]; then
    IFS=':' read -r -a T1_INPUTS <<< "$UNIPROT_T1_INPUTS"
else
    T1_INPUTS+=("$T1_UNIPROT_DIR/uniprot_sprot.dat.gz")
    if [[ -f "$T1_UNIPROT_DIR/uniprot_trembl.dat.gz" ]]; then
        T1_INPUTS+=("$T1_UNIPROT_DIR/uniprot_trembl.dat.gz")
    elif [[ "$ALLOW_SPROT_ONLY" != "1" ]]; then
        echo "Missing t1 TrEMBL input. Set ALLOW_SPROT_ONLY=1 only for a diagnostic Swiss-Prot-only build." >&2
        exit 1
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
    --go-obo "$T0_OBO"
    --go-obo-t0 "$T0_OBO"
    --go-obo-t1 "$T1_OBO"
    --output-dir "$OUTPUT_DIR"
    --report-dir "$REPORT_DIR"
)
for path in "${T0_INPUTS[@]}"; do
    COMMAND+=(--uniprot-t0 "$path")
done
for path in "${T1_INPUTS[@]}"; do
    COMMAND+=(--uniprot-t1 "$path")
done
if [[ "$ALLOW_SPROT_ONLY" == "1" ]]; then
    COMMAND+=(--target-reviewed-only)
fi

printf '%q ' "${COMMAND[@]}" > "$LOG_DIR/command.txt"
printf '\n' >> "$LOG_DIR/command.txt"
{
    echo "profile=$PROFILE"
    echo "db_root=$DB_ROOT"
    echo "run_root=$RUN_ROOT"
    echo "python=$($PYTHON_BIN --version 2>&1)"
    echo "hostname=$(hostname)"
} > "$LOG_DIR/environment.txt"

echo "Running contemporary CAFA benchmark builder"
echo "Profile : $PROFILE"
echo "Outputs : $OUTPUT_DIR"
echo "Reports : $REPORT_DIR"
echo

export PYTHONPATH="$BUILDER_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export CAFA_BUILDER_USE_PIGZ="${CAFA_BUILDER_USE_PIGZ:-1}"
"${COMMAND[@]}" 2>&1 | tee "$LOG_DIR/builder.log"

echo
echo "Benchmark build completed successfully: $RUN_ROOT"
