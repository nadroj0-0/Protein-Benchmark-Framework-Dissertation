#!/bin/bash
# Populate the dissertation SAN with authenticated, frozen public inputs.
#
# This is deliberately separate from protein_database_download.sh. That older
# script documents the historical home-directory acquisition workflow; this
# script owns the persistent /SAN/bioinf/bmpfp contract.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRAMEWORK_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SPEC_FILE="${SAN_INPUT_SPEC:-${SCRIPT_DIR}/san_frozen_inputs.tsv}"
FILTER_DAT="$FRAMEWORK_ROOT/scripts/benchmark_generation/filter_uniprot_dat.py"
EXTRACT_MEMBER="$FRAMEWORK_ROOT/scripts/benchmark_generation/extract_tar_member.py"
TARGET_TAXA="$FRAMEWORK_ROOT/benchmark_builders/contemporary_cafa/src/cafa_benchmark_builder/resources/cafa3_target_taxa.txt"
# shellcheck source=../reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"
DEFAULT_ROOT="/SAN/bioinf/bmpfp"
ROOT="$DEFAULT_ROOT"
RESERVE_GB=40
DRY_RUN=0
VERIFY_ONLY=0
FULL_VERIFY=0
LIST_ONLY=0
PROFILES=()
DOWNLOADED=0
SKIPPED=0
VERIFIED=0
DERIVED_CREATED=0
DERIVED_SKIPPED=0
FALLBACK_LOCK_DIR=""

DERIVED_T0_ROLE="uniprot_trembl_cafa3_targets_t0"
DERIVED_T0_RELEASE="2025_01"
DERIVED_T0_RELATIVE="derived_inputs/uniprot/cafa3_target_taxa/2025_01/uniprot_trembl_cafa3_targets.dat.gz"
DERIVED_T1_ROLE="uniprot_trembl_cafa3_targets_t1"
DERIVED_T1_RELEASE="2026_02"
DERIVED_T1_RELATIVE="derived_inputs/uniprot/cafa3_target_taxa/2026_02/uniprot_trembl_cafa3_targets.dat.gz"
HOMOLOGY_CACHE_ROLE="homology_common_preprocessing_2026_02"
HOMOLOGY_CACHE_SCOPE="sprot-and-trembl"
HOMOLOGY_CACHE_RELATIVE="derived_inputs/homology/2026_02/goa_234/sprot-and-trembl/common_preprocessing"
HOMOLOGY_CACHE_ALLOWANCE_GB="${HOMOLOGY_CACHE_ALLOWANCE_GB:-50}"
CACHE_MARKER="CACHE_COMPLETE.json"

usage() {
    cat <<'EOF'
Usage:
  bash scripts/data_acquisition/populate_san_frozen_inputs.sh [options]

Options:
  --root PATH           Persistent store root (default: /SAN/bioinf/bmpfp).
  --profile NAME        Select a profile; repeat or comma-separate values.
                        Choices: temporal, homology, embedding-inputs,
                        references, tools, all. Default: all.
  --reserve-gb N        Space that must remain after missing downloads
                        (default: 40).
  --dry-run             Print the acquisition plan without writing or fetching.
  --verify-only         Require every selected file and fully verify it offline.
  --full-verify         Rehash and structurally inspect existing files. Newly
                        downloaded and derived files are always fully verified.
  --list                List selected catalogue entries and exit.
  --help                Show this help.

Examples:
  # Show the complete plan and required bytes.
  bash scripts/data_acquisition/populate_san_frozen_inputs.sh --dry-run

  # Populate only the 2026 homology inputs and pinned MMseqs2 tool archive.
  bash scripts/data_acquisition/populate_san_frozen_inputs.sh \
    --profile homology --profile tools

  # Offline integrity audit of everything already held in SAN.
  bash scripts/data_acquisition/populate_san_frozen_inputs.sh \
    --profile all --verify-only
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --root)
            [[ $# -ge 2 ]] || die "--root requires a path"
            ROOT="$2"
            shift 2
            ;;
        --profile)
            [[ $# -ge 2 ]] || die "--profile requires a value"
            IFS=',' read -r -a requested_profiles <<< "$2"
            for requested_profile in "${requested_profiles[@]}"; do
                [[ -n "$requested_profile" ]] || die "Empty --profile value"
                PROFILES+=("$requested_profile")
            done
            shift 2
            ;;
        --reserve-gb)
            [[ $# -ge 2 ]] || die "--reserve-gb requires an integer"
            RESERVE_GB="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --verify-only)
            VERIFY_ONLY=1
            FULL_VERIFY=1
            shift
            ;;
        --full-verify)
            FULL_VERIFY=1
            shift
            ;;
        --list)
            LIST_ONLY=1
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

[[ -f "$SPEC_FILE" ]] || die "Input specification not found: $SPEC_FILE"
[[ "$RESERVE_GB" =~ ^[0-9]+$ ]] || die "--reserve-gb must be a non-negative integer"
[[ "$HOMOLOGY_CACHE_ALLOWANCE_GB" =~ ^[1-9][0-9]*$ ]] || \
    die "HOMOLOGY_CACHE_ALLOWANCE_GB must be a positive integer"
[[ "$ROOT" = /* ]] || die "--root must be an absolute path"
[[ "$ROOT" != "/" ]] || die "Refusing to use / as the store root"
if [[ ${#PROFILES[@]} -eq 0 ]]; then
    PROFILES=("all")
fi

for profile in "${PROFILES[@]}"; do
    case "$profile" in
        temporal|homology|embedding-inputs|references|tools|all) ;;
        *) die "Unknown profile: $profile" ;;
    esac
done

profile_selected() {
    local row_profiles="$1"
    local selected
    local row_profile
    local row_values=()
    for selected in "${PROFILES[@]}"; do
        [[ "$selected" == "all" ]] && return 0
    done
    IFS=',' read -r -a row_values <<< "$row_profiles"
    for selected in "${PROFILES[@]}"; do
        for row_profile in "${row_values[@]}"; do
            [[ "$selected" == "$row_profile" ]] && return 0
        done
    done
    return 1
}

spec_value_for_role() {
    local wanted_role="$1"
    local field_number="$2"
    awk -F '\t' -v wanted="$wanted_role" -v field="$field_number" '
        $1 !~ /^#/ && $2 == wanted {print $field; found=1; exit}
        END {exit !found}
    ' "$SPEC_FILE"
}

file_size() {
    if stat -c '%s' "$1" >/dev/null 2>&1; then
        stat -c '%s' "$1"
    else
        stat -f '%z' "$1"
    fi
}

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

md5_file() {
    if command -v md5sum >/dev/null 2>&1; then
        md5sum "$1" | awk '{print $1}'
    else
        md5 -q "$1"
    fi
}

digest_file() {
    local algorithm="$1"
    local path="$2"
    case "$algorithm" in
        sha256) sha256_file "$path" ;;
        md5) md5_file "$path" ;;
        *) die "Unsupported checksum algorithm: $algorithm" ;;
    esac
}

human_gib() {
    awk -v bytes="$1" 'BEGIN {printf "%.2f GiB", bytes / 1073741824}'
}

validate_spec() {
    local profiles role release relative_path url expected_bytes checksum_algorithm
    local expected_checksum validator
    local seen_roles=""
    local seen_paths=""
    while IFS=$'\t' read -r profiles role release relative_path url expected_bytes \
        checksum_algorithm expected_checksum validator; do
        [[ -n "$profiles" ]] || continue
        [[ "$profiles" == \#* ]] && continue
        [[ -n "$role" && -n "$release" && -n "$relative_path" && -n "$url" ]] || \
            die "Malformed row in $SPEC_FILE"
        [[ "$role" =~ ^[A-Za-z0-9._-]+$ ]] || \
            die "Unsafe artifact role in $SPEC_FILE: $role"
        [[ "$relative_path" != /* && "$relative_path" != *".."* ]] || \
            die "Unsafe relative path in $SPEC_FILE: $relative_path"
        [[ "$expected_bytes" == "-" || "$expected_bytes" =~ ^[1-9][0-9]*$ ]] || \
            die "Invalid expected byte count for $role: $expected_bytes"
        case "$checksum_algorithm" in
            -)
                [[ "$expected_checksum" == "-" ]] || \
                    die "Checksum value without algorithm for $role"
                ;;
            sha256)
                [[ "$expected_checksum" =~ ^[0-9a-f]{64}$ ]] || \
                    die "Invalid SHA-256 for $role"
                ;;
            md5)
                [[ "$expected_checksum" =~ ^[0-9a-f]{32}$ ]] || \
                    die "Invalid MD5 for $role"
                ;;
            *) die "Unsupported checksum algorithm for $role: $checksum_algorithm" ;;
        esac
        case "$validator" in
            text|csv|obo|gzip|tar-gzip|uniprot-dat-gzip|fasta-gzip|goa-gzip|\
            md5-text|hdf5|deepgo-archive|mmseqs-archive) ;;
            *) die "Unknown validator for $role: $validator" ;;
        esac
        [[ "$seen_roles" != *"|${role}|"* ]] || die "Duplicate role in $SPEC_FILE: $role"
        [[ "$seen_paths" != *"|${relative_path}|"* ]] || \
            die "Duplicate path in $SPEC_FILE: $relative_path"
        seen_roles+="|${role}|"
        seen_paths+="|${relative_path}|"
    done < "$SPEC_FILE"
}

validate_artifact() {
    local path="$1"
    local validator="$2"
    local listing
    [[ -s "$path" ]] || die "Artifact is missing or empty: $path"
    case "$validator" in
        text)
            LC_ALL=C grep -q '[^[:space:]]' "$path" || die "Text file is blank: $path"
            ;;
        csv)
            # Zenodo 7409660's canonical mf-training.csv alone uses singular
            # "protein"; retain the authenticated raw bytes and normalize only
            # in the downstream PFP preparation workspace.
            head -n 1 "$path" | grep -Eq '^(protein|proteins),sequences,' || \
                die "Unexpected PFP CSV header: $path"
            ;;
        obo)
            grep -q '^format-version:' "$path" || die "Missing OBO format header: $path"
            grep -q '^\[Term\]' "$path" || die "OBO contains no terms: $path"
            ;;
        gzip)
            gzip -t "$path"
            ;;
        tar-gzip)
            tar -tzf "$path" >/dev/null
            ;;
        uniprot-dat-gzip)
            gzip -t "$path"
            zgrep -m 1 '^ID   ' "$path" >/dev/null || die "No UniProt DAT records: $path"
            ;;
        fasta-gzip)
            gzip -t "$path"
            zgrep -m 1 '^>' "$path" >/dev/null || die "No FASTA records: $path"
            ;;
        goa-gzip)
            gzip -t "$path"
            zgrep -m 1 '^!gaf-version:' "$path" >/dev/null || \
                die "Missing GAF header: $path"
            ;;
        md5-text)
            grep -Eq '^[0-9a-fA-F]{32}([[:space:]]|$)' "$path" || \
                die "Invalid MD5 sidecar: $path"
            ;;
        hdf5)
            [[ "$(od -An -t x1 -N 8 "$path" | tr -d ' \n')" == "894844460d0a1a0a" ]] || \
                die "Invalid HDF5 signature: $path"
            ;;
        deepgo-archive)
            listing="${path}.listing.$$.tmp"
            tar -tzf "$path" > "$listing"
            grep -Eq '(^|/)train_data\.pkl$' "$listing" || \
                die "DeepGOPlus archive lacks train_data.pkl: $path"
            grep -Eq '(^|/)test_data\.pkl$' "$listing" || \
                die "DeepGOPlus archive lacks test_data.pkl: $path"
            grep -Eq '(^|/)terms\.pkl$' "$listing" || \
                die "DeepGOPlus archive lacks terms.pkl: $path"
            rm -f "$listing"
            ;;
        mmseqs-archive)
            listing="${path}.listing.$$.tmp"
            tar -tzf "$path" > "$listing"
            grep -Eq '(^|/)mmseqs/bin/mmseqs$' "$listing" || \
                die "MMseqs2 archive lacks mmseqs/bin/mmseqs: $path"
            rm -f "$listing"
            ;;
    esac
}

verify_expected_size() {
    local path="$1"
    local expected_bytes="$2"
    local observed_bytes
    [[ "$expected_bytes" == "-" ]] && return 0
    observed_bytes="$(file_size "$path")"
    [[ "$observed_bytes" == "$expected_bytes" ]] || \
        die "Size mismatch for $path: expected=$expected_bytes observed=$observed_bytes"
}

verify_expected_checksum() {
    local path="$1"
    local algorithm="$2"
    local expected_checksum="$3"
    local observed_checksum
    [[ "$algorithm" == "-" ]] && return 0
    observed_checksum="$(digest_file "$algorithm" "$path")"
    [[ "$observed_checksum" == "$expected_checksum" ]] || \
        die "$algorithm mismatch for $path: expected=$expected_checksum observed=$observed_checksum"
}

write_artifact_metadata() {
    local role="$1"
    local release="$2"
    local relative_path="$3"
    local url="$4"
    local expected_bytes="$5"
    local checksum_algorithm="$6"
    local expected_checksum="$7"
    local path="$8"
    local acquisition="$9"
    local supplied_sha="${10:-}"
    local observed_bytes observed_sha timestamp sha_tmp provenance_tmp
    observed_bytes="$(file_size "$path")"
    if [[ -n "$supplied_sha" ]]; then
        observed_sha="$supplied_sha"
    else
        observed_sha="$(sha256_file "$path")"
    fi
    timestamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    sha_tmp="${path}.sha256.partial.$$"
    provenance_tmp="${path}.provenance.tsv.partial.$$"
    printf '%s  %s\n' "$observed_sha" "$(basename "$path")" > "$sha_tmp"
    mv "$sha_tmp" "${path}.sha256"
    printf 'role\trelease\trelative_path\turl\texpected_bytes\tobserved_bytes\texpected_checksum_algorithm\texpected_checksum\tobserved_sha256\tacquisition\tacquired_utc\n' \
        > "$provenance_tmp"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$role" "$release" "$relative_path" "$url" "$expected_bytes" "$observed_bytes" \
        "$checksum_algorithm" "$expected_checksum" "$observed_sha" "$acquisition" "$timestamp" \
        >> "$provenance_tmp"
    mv "$provenance_tmp" "${path}.provenance.tsv"
}

verify_recorded_sha256() {
    local path="$1"
    local sidecar="${path}.sha256"
    local expected observed
    [[ -s "$sidecar" ]] || die "Missing SHA-256 sidecar for unpinned artifact: $path"
    expected="$(awk 'NR == 1 {print $1}' "$sidecar")"
    [[ "$expected" =~ ^[0-9a-f]{64}$ ]] || die "Invalid SHA-256 sidecar: $sidecar"
    observed="$(sha256_file "$path")"
    [[ "$observed" == "$expected" ]] || \
        die "Recorded SHA-256 mismatch for $path: expected=$expected observed=$observed"
}

recorded_sha256_value() {
    local path="$1"
    local sidecar="${path}.sha256"
    local recorded
    [[ -s "$sidecar" ]] || die "Missing SHA-256 sidecar: $sidecar"
    recorded="$(awk 'NR == 1 {print $1}' "$sidecar")"
    [[ "$recorded" =~ ^[0-9a-f]{64}$ ]] || die "Invalid SHA-256 sidecar: $sidecar"
    printf '%s\n' "$recorded"
}

derivation_value() {
    local path="$1"
    local key="$2"
    awk -F '\t' -v wanted="$key" '
        $1 == wanted {print $2; found=1; exit}
        END {exit !found}
    ' "$path"
}

require_derivation_value() {
    local path="$1"
    local key="$2"
    local expected="$3"
    local observed
    observed="$(derivation_value "$path" "$key")" || \
        die "Derived provenance lacks $key: $path"
    [[ "$observed" == "$expected" ]] || \
        die "Derived provenance mismatch for $key in $path: expected=$expected observed=$observed"
}

count_uniprot_records() {
    gzip -dc "$1" | awk '/^ID   / {count++} END {print count + 0}'
}

resolve_python_bin() {
    local candidate
    local candidates=()
    [[ -z "${PYTHON_BIN:-}" ]] || candidates+=("$PYTHON_BIN")
    candidates+=(python3 python)
    for candidate in "${candidates[@]}"; do
        command -v "$candidate" >/dev/null 2>&1 || continue
        if "$candidate" -c 'import sys; raise SystemExit(sys.version_info < (3, 9))' \
            >/dev/null 2>&1; then
            PYTHON_BIN="$(command -v "$candidate")"
            return 0
        fi
    done
    die "Python >=3.9 is required for derived inputs; activate mmfp or set PYTHON_BIN"
}

write_derived_metadata() {
    local role="$1"
    local release="$2"
    local relative_path="$3"
    local source_role="$4"
    local source_relative="$5"
    local source_path="$6"
    local archive_member="$7"
    local output_path="$8"
    local record_count="$9"
    local source_sha taxa_sha filter_sha output_sha output_bytes timestamp
    local derivation_tmp

    source_sha="$(recorded_sha256_value "$source_path")"
    taxa_sha="$(sha256_file "$TARGET_TAXA")"
    filter_sha="$(sha256_file "$FILTER_DAT")"
    output_sha="$(sha256_file "$output_path")"
    output_bytes="$(file_size "$output_path")"
    timestamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

    write_artifact_metadata "$role" "$release" "$relative_path" \
        "derived://${source_role}#${archive_member}" "-" "-" "-" "$output_path" \
        "derived-cafa3-target-taxa-filter" "$output_sha"

    derivation_tmp="${output_path}.derivation.tsv.partial.$$"
    printf 'schema_version\t1\n' > "$derivation_tmp"
    printf 'role\t%s\n' "$role" >> "$derivation_tmp"
    printf 'release\t%s\n' "$release" >> "$derivation_tmp"
    printf 'source_artifact_id\t%s\n' "$source_role" >> "$derivation_tmp"
    printf 'source_relative_path\t%s\n' "$source_relative" >> "$derivation_tmp"
    printf 'source_sha256\t%s\n' "$source_sha" >> "$derivation_tmp"
    printf 'archive_member\t%s\n' "$archive_member" >> "$derivation_tmp"
    printf 'target_taxa_sha256\t%s\n' "$taxa_sha" >> "$derivation_tmp"
    printf 'filter_script_sha256\t%s\n' "$filter_sha" >> "$derivation_tmp"
    printf 'output_sha256\t%s\n' "$output_sha" >> "$derivation_tmp"
    printf 'output_bytes\t%s\n' "$output_bytes" >> "$derivation_tmp"
    printf 'record_count\t%s\n' "$record_count" >> "$derivation_tmp"
    printf 'generated_utc\t%s\n' "$timestamp" >> "$derivation_tmp"
    mv "$derivation_tmp" "${output_path}.derivation.tsv"
}

verify_derived_contract() {
    local role="$1"
    local release="$2"
    local source_role="$3"
    local source_relative="$4"
    local source_path="$5"
    local archive_member="$6"
    local output_path="$7"
    local derivation="${output_path}.derivation.tsv"
    local source_sha taxa_sha filter_sha output_sha output_bytes record_count

    source_sha="$(recorded_sha256_value "$source_path")"
    taxa_sha="$(sha256_file "$TARGET_TAXA")"
    filter_sha="$(sha256_file "$FILTER_DAT")"
    output_sha="$(recorded_sha256_value "$output_path")"
    output_bytes="$(file_size "$output_path")"

    require_derivation_value "$derivation" schema_version 1
    require_derivation_value "$derivation" role "$role"
    require_derivation_value "$derivation" release "$release"
    require_derivation_value "$derivation" source_artifact_id "$source_role"
    require_derivation_value "$derivation" source_relative_path "$source_relative"
    require_derivation_value "$derivation" source_sha256 "$source_sha"
    require_derivation_value "$derivation" archive_member "$archive_member"
    require_derivation_value "$derivation" target_taxa_sha256 "$taxa_sha"
    require_derivation_value "$derivation" filter_script_sha256 "$filter_sha"
    require_derivation_value "$derivation" output_sha256 "$output_sha"
    require_derivation_value "$derivation" output_bytes "$output_bytes"

    if [[ "$FULL_VERIFY" == "1" ]]; then
        verify_recorded_sha256 "$output_path"
        validate_artifact "$output_path" uniprot-dat-gzip
        record_count="$(count_uniprot_records "$output_path")"
        require_derivation_value "$derivation" record_count "$record_count"
        VERIFIED=$((VERIFIED + 1))
    fi
}

build_or_verify_derived_trembl() {
    local role="$1"
    local release="$2"
    local relative_path="$3"
    local source_role="$4"
    local archive_member="$5"
    local source_relative source_path destination partial record_count

    source_relative="$(spec_value_for_role "$source_role" 4)" || \
        die "Derived input $role requires missing specification role: $source_role"
    source_path="$ROOT/$source_relative"
    destination="$ROOT/$relative_path"
    [[ -s "$source_path" ]] || die "Derived input source is missing: $source_path"
    [[ -s "${source_path}.sha256" && -s "${source_path}.provenance.tsv" ]] || \
        die "Derived input source is not authenticated: $source_path"

    echo
    echo "[$role] $relative_path"
    if [[ -s "$destination" && -s "${destination}.sha256" && \
          -s "${destination}.provenance.tsv" && -s "${destination}.derivation.tsv" ]]; then
        verify_derived_contract "$role" "$release" "$source_role" "$source_relative" \
            "$source_path" "$archive_member" "$destination"
        echo "  present with matching derivation contract; skipping generation"
        DERIVED_SKIPPED=$((DERIVED_SKIPPED + 1))
        return 0
    fi

    if [[ "$VERIFY_ONLY" == "1" ]]; then
        die "Selected derived artifact is missing or incomplete: $destination"
    fi

    add_mmfp_singularity_bind "$ROOT"
    resolve_python_bin
    command -v gzip >/dev/null 2>&1 || die "gzip is required to derive filtered TrEMBL caches"
    [[ -f "$FILTER_DAT" ]] || die "UniProt DAT filter is missing: $FILTER_DAT"
    [[ -f "$EXTRACT_MEMBER" ]] || die "Archive member extractor is missing: $EXTRACT_MEMBER"
    [[ -f "$TARGET_TAXA" ]] || die "CAFA3 target taxa file is missing: $TARGET_TAXA"

    mkdir -p "$(dirname "$destination")"
    partial="${destination}.partial.$$"
    rm -f "$partial"
    echo "  deriving CAFA3-target-taxa TrEMBL cache from $source_path"
    if [[ "$archive_member" == "-" ]]; then
        if ! gzip -dc "$source_path" \
            | "$PYTHON_BIN" "$FILTER_DAT" --taxa-file "$TARGET_TAXA" \
            | gzip -n -c > "$partial"; then
            rm -f "$partial"
            die "Failed to derive $role"
        fi
    else
        if ! "$PYTHON_BIN" "$EXTRACT_MEMBER" --archive "$source_path" --suffix "$archive_member" \
            | gzip -dc \
            | "$PYTHON_BIN" "$FILTER_DAT" --taxa-file "$TARGET_TAXA" \
            | gzip -n -c > "$partial"; then
            rm -f "$partial"
            die "Failed to derive $role"
        fi
    fi

    validate_artifact "$partial" uniprot-dat-gzip
    record_count="$(count_uniprot_records "$partial")"
    (( record_count > 0 )) || die "Derived TrEMBL cache contains zero records: $partial"
    mv "$partial" "$destination"
    write_derived_metadata "$role" "$release" "$relative_path" "$source_role" \
        "$source_relative" "$source_path" "$archive_member" "$destination" "$record_count"
    echo "  derived and published atomically: records=$record_count"
    DERIVED_CREATED=$((DERIVED_CREATED + 1))
    VERIFIED=$((VERIFIED + 1))
}

process_temporal_derived_inputs() {
    profile_selected temporal || return 0
    build_or_verify_derived_trembl "$DERIVED_T0_ROLE" "$DERIVED_T0_RELEASE" \
        "$DERIVED_T0_RELATIVE" uniprot_knowledgebase_t0 uniprot_trembl.dat.gz
    build_or_verify_derived_trembl "$DERIVED_T1_ROLE" "$DERIVED_T1_RELEASE" \
        "$DERIVED_T1_RELATIVE" uniprot_trembl_t1 -
}

homology_input_path() {
    local role="$1"
    local relative
    relative="$(spec_value_for_role "$role" 4)" || \
        die "Homology common cache requires missing specification role: $role"
    printf '%s/%s\n' "$ROOT" "$relative"
}

homology_input_url() {
    spec_value_for_role "$1" 5 || \
        die "Homology common cache requires a URL for role: $1"
}

homology_expected_sha_arguments() {
    local spec_role logical_role source_path
    for binding in \
        "uniref90_t1|uniref90_fasta" \
        "idmapping_t1|idmapping" \
        "uniprot_sprot_t1|uniprot_sprot_sequences" \
        "uniprot_trembl_t1|uniprot_trembl_sequences" \
        "goa_t1|goa" \
        "go_basic_t1|go_obo"; do
        IFS='|' read -r spec_role logical_role <<< "$binding"
        source_path="$(homology_input_path "$spec_role")"
        printf '%s=%s\n' "$logical_role" "$(recorded_sha256_value "$source_path")"
    done
}

process_homology_derived_inputs() {
    profile_selected homology || return 0
    local destination="$ROOT/$HOMOLOGY_CACHE_RELATIVE"
    local marker="$destination/CACHE_COMPLETE.json"
    local work_root manifest policy revision binding
    local verify_args=()
    local build_args=()
    local expected_bindings=()

    echo
    echo "[$HOMOLOGY_CACHE_ROLE] $HOMOLOGY_CACHE_RELATIVE"
    resolve_python_bin
    add_mmfp_singularity_bind "$ROOT"
    export PYTHONPATH="$FRAMEWORK_ROOT/benchmark_builders/homology_cluster/src${PYTHONPATH:+:$PYTHONPATH}"
    export MMFP_PYTHONPATH="$PYTHONPATH"
    while IFS= read -r binding; do
        expected_bindings+=("$binding")
        verify_args+=(--expected-input-sha256 "$binding")
    done < <(homology_expected_sha_arguments)
    verify_args+=(
        --cache-dir "$destination"
        --source-scope "$HOMOLOGY_CACHE_SCOPE"
    )
    [[ "$FULL_VERIFY" == "1" ]] && verify_args+=(--full-hashes)

    if [[ -s "$marker" ]] && "$PYTHON_BIN" -m homology_cluster_benchmark.common_cache \
        verify "${verify_args[@]}"; then
        echo "  present with matching input, policy, code, and file contracts; skipping generation"
        DERIVED_SKIPPED=$((DERIVED_SKIPPED + 1))
        [[ "$FULL_VERIFY" == "1" ]] && VERIFIED=$((VERIFIED + 1))
        return 0
    fi
    if [[ "$VERIFY_ONLY" == "1" ]]; then
        die "Selected homology common preprocessing cache is missing or invalid: $destination"
    fi

    revision="${FRAMEWORK_REVISION:-}"
    if [[ -z "$revision" ]]; then
        revision="$(cd "$FRAMEWORK_ROOT" && git rev-parse HEAD)" || \
            die "FRAMEWORK_REVISION is unset and the framework is not a Git checkout"
    fi
    [[ "$revision" =~ ^[0-9a-f]{40}$ ]] || \
        die "Homology cache generation requires a 40-character FRAMEWORK_REVISION"
    work_root="${HOMOLOGY_CACHE_WORK_DIR:-${TMPDIR:-/tmp}/homology-common-cache-${USER:-user}-$$}"
    [[ "$work_root" = /* && "$work_root" != "/" ]] || \
        die "HOMOLOGY_CACHE_WORK_DIR must be an absolute non-root path"
    mkdir -p "$work_root/contracts" "$work_root/build"
    manifest="$work_root/contracts/frozen_input_manifest.json"
    policy="$work_root/contracts/unused_runtime_policy.json"

    echo "  preparing authoritative frozen-input contract"
    "$PYTHON_BIN" -m homology_cluster_benchmark.runtime_contract prepare \
        --manifest-out "$manifest" \
        --policy-out "$policy" \
        --source-scope "$HOMOLOGY_CACHE_SCOPE" \
        --framework-revision "$revision" \
        --uniref90-fasta "$(homology_input_path uniref90_t1)" \
        --uniref90-fasta-url "$(homology_input_url uniref90_t1)" \
        --uniref90-fasta-acquisition provided-persistent-store \
        --idmapping "$(homology_input_path idmapping_t1)" \
        --idmapping-url "$(homology_input_url idmapping_t1)" \
        --idmapping-acquisition provided-persistent-store \
        --uniprot-sprot-sequences "$(homology_input_path uniprot_sprot_t1)" \
        --uniprot-sprot-sequences-url "$(homology_input_url uniprot_sprot_t1)" \
        --uniprot-sprot-sequences-acquisition provided-persistent-store \
        --uniprot-trembl-sequences "$(homology_input_path uniprot_trembl_t1)" \
        --uniprot-trembl-sequences-url "$(homology_input_url uniprot_trembl_t1)" \
        --uniprot-trembl-sequences-acquisition provided-persistent-store \
        --goa "$(homology_input_path goa_t1)" \
        --goa-url "$(homology_input_url goa_t1)" \
        --goa-acquisition provided-persistent-store \
        --go-obo "$(homology_input_path go_basic_t1)" \
        --go-obo-url "$(homology_input_url go_basic_t1)" \
        --go-obo-acquisition provided-persistent-store

    build_args=(
        --output-dir "$destination"
        --work-dir "$work_root/build"
        --frozen-input-manifest "$manifest"
        --source-scope "$HOMOLOGY_CACHE_SCOPE"
        --uniref90-fasta "$(homology_input_path uniref90_t1)"
        --idmapping "$(homology_input_path idmapping_t1)"
        --uniprot-sprot-sequences "$(homology_input_path uniprot_sprot_t1)"
        --uniprot-trembl-sequences "$(homology_input_path uniprot_trembl_t1)"
        --goa "$(homology_input_path goa_t1)"
        --go-obo "$(homology_input_path go_basic_t1)"
        --replace-existing
    )
    echo "  building threshold-independent preprocessing state once"
    "$PYTHON_BIN" -m homology_cluster_benchmark.common_cache build "${build_args[@]}"
    local final_verify_args=(
        --cache-dir "$destination"
        --source-scope "$HOMOLOGY_CACHE_SCOPE"
        --full-hashes
    )
    for binding in "${expected_bindings[@]}"; do
        final_verify_args+=(--expected-input-sha256 "$binding")
    done
    "$PYTHON_BIN" -m homology_cluster_benchmark.common_cache \
        verify "${final_verify_args[@]}"
    echo "  derived and published atomically: $destination"
    DERIVED_CREATED=$((DERIVED_CREATED + 1))
    VERIFIED=$((VERIFIED + 1))
}

fetch_text() {
    local url="$1"
    if command -v curl >/dev/null 2>&1; then
        if curl --fail --location --silent --show-error "$url"; then
            return 0
        fi
        echo "curl failed while fetching release metadata; trying wget" >&2
    fi
    if command -v wget >/dev/null 2>&1; then
        wget --quiet --output-document=- "$url"
    else
        die "curl failed and wget is not available"
    fi
}

download_to_partial() {
    local url="$1"
    local partial="$2"
    if [[ -s "$partial" ]]; then
        echo "  resuming partial: $partial"
    fi
    if command -v wget >/dev/null 2>&1; then
        if wget --continue --tries=5 --timeout=60 --progress=dot:giga \
            --output-document="$partial" "$url"; then
            [[ -s "$partial" ]] || die "Download is empty: $url"
            return 0
        fi
        echo "  wget failed; preserving the partial and trying curl" >&2
    fi
    if command -v curl >/dev/null 2>&1; then
        curl --fail --location --retry 5 --continue-at - --output "$partial" "$url"
    else
        die "wget failed and curl is not available"
    fi
    [[ -s "$partial" ]] || die "Download is empty: $url"
}

catalog_metadata() {
    local catalog="$ROOT/manifests/frozen_input_catalog.tsv"
    local temporary="${catalog}.partial.$$"
    local sidecar
    printf 'role\trelease\trelative_path\turl\texpected_bytes\tobserved_bytes\texpected_checksum_algorithm\texpected_checksum\tobserved_sha256\tacquisition\tacquired_utc\n' \
        > "$temporary"
    while IFS= read -r sidecar; do
        tail -n +2 "$sidecar"
    done < <(find "$ROOT" -type f -name '*.provenance.tsv' -print | sort) | sort -t $'\t' -k3,3 \
        >> "$temporary"
    if [[ -f "$catalog" ]] && cmp -s "$temporary" "$catalog"; then
        rm -f "$temporary"
    else
        mv "$temporary" "$catalog"
    fi
}

write_artifact_path_catalog() {
    local catalog="$ROOT/manifests/artifact_paths.tsv"
    local temporary="${catalog}.partial.$$"
    local profiles role release relative_path url expected_bytes
    local checksum_algorithm expected_checksum validator destination
    printf 'artifact_id\tpath\n' > "$temporary"
    while IFS=$'\t' read -r profiles role release relative_path url expected_bytes \
        checksum_algorithm expected_checksum validator; do
        [[ -n "$profiles" ]] || continue
        [[ "$profiles" == \#* ]] && continue
        destination="$ROOT/$relative_path"
        [[ -s "$destination" && -s "${destination}.sha256" && \
           -s "${destination}.provenance.tsv" ]] || continue
        printf '%s\t%s\n' "$role" "$destination" >> "$temporary"
    done < "$SPEC_FILE"
    for derived_entry in \
        "$DERIVED_T0_ROLE|$DERIVED_T0_RELATIVE" \
        "$DERIVED_T1_ROLE|$DERIVED_T1_RELATIVE"; do
        IFS='|' read -r role relative_path <<< "$derived_entry"
        destination="$ROOT/$relative_path"
        [[ -s "$destination" && -s "${destination}.sha256" && \
           -s "${destination}.provenance.tsv" && -s "${destination}.derivation.tsv" ]] || continue
        printf '%s\t%s\n' "$role" "$destination" >> "$temporary"
    done
    destination="$ROOT/$HOMOLOGY_CACHE_RELATIVE/$CACHE_MARKER"
    if [[ -s "$destination" ]]; then
        printf '%s\t%s\n' "$HOMOLOGY_CACHE_ROLE" "$destination" >> "$temporary"
    fi
    if [[ -f "$catalog" ]] && cmp -s "$temporary" "$catalog"; then
        rm -f "$temporary"
    else
        mv "$temporary" "$catalog"
    fi
}

release_lock() {
    if [[ -n "$FALLBACK_LOCK_DIR" ]]; then
        rmdir "$FALLBACK_LOCK_DIR" 2>/dev/null || true
    fi
}

validate_spec

selected_count=0
missing_count=0
missing_known_bytes=0
missing_unknown_count=0
needs_uniprot_guard=0
needs_goa_guard=0

while IFS=$'\t' read -r profiles role release relative_path url expected_bytes \
    checksum_algorithm expected_checksum validator; do
    [[ -n "$profiles" ]] || continue
    [[ "$profiles" == \#* ]] && continue
    profile_selected "$profiles" || continue
    selected_count=$((selected_count + 1))
    destination="$ROOT/$relative_path"
    if [[ ! -s "$destination" ]]; then
        missing_count=$((missing_count + 1))
        if [[ "$expected_bytes" == "-" ]]; then
            missing_unknown_count=$((missing_unknown_count + 1))
        else
            remaining_bytes="$expected_bytes"
            if [[ -s "${destination}.partial" ]]; then
                partial_bytes="$(file_size "${destination}.partial")"
                (( partial_bytes <= expected_bytes )) || \
                    die "Partial file is larger than expected: ${destination}.partial"
                remaining_bytes=$((expected_bytes - partial_bytes))
            fi
            missing_known_bytes=$((missing_known_bytes + remaining_bytes))
        fi
        [[ "$url" == *"ftp.uniprot.org/pub/databases/uniprot/current_release"* ]] && \
            needs_uniprot_guard=1
        [[ "$url" == *"ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT"* || \
           "$url" == *"ftp.ebi.ac.uk/pub/databases/GO/goa/current_release_numbers.txt"* ]] && \
            needs_goa_guard=1
    fi
    if [[ "$LIST_ONLY" == "1" || "$DRY_RUN" == "1" ]]; then
        printf '%-18s %-32s %-16s %s\n' "$profiles" "$role" "$release" "$destination"
    fi
done < "$SPEC_FILE"

if profile_selected temporal; then
    selected_count=$((selected_count + 2))
    for derived_entry in \
        "$DERIVED_T0_ROLE|$DERIVED_T0_RELATIVE" \
        "$DERIVED_T1_ROLE|$DERIVED_T1_RELATIVE"; do
        IFS='|' read -r role relative_path <<< "$derived_entry"
        destination="$ROOT/$relative_path"
        if [[ ! -s "$destination" || ! -s "${destination}.sha256" || \
              ! -s "${destination}.provenance.tsv" || ! -s "${destination}.derivation.tsv" ]]; then
            missing_count=$((missing_count + 1))
            missing_unknown_count=$((missing_unknown_count + 1))
        fi
        if [[ "$LIST_ONLY" == "1" || "$DRY_RUN" == "1" ]]; then
            printf '%-18s %-32s %-16s %s\n' temporal "$role" derived "$destination"
        fi
    done
fi

if profile_selected homology; then
    selected_count=$((selected_count + 1))
    destination="$ROOT/$HOMOLOGY_CACHE_RELATIVE/$CACHE_MARKER"
    if [[ ! -s "$destination" ]]; then
        missing_count=$((missing_count + 1))
        missing_known_bytes=$((
            missing_known_bytes + HOMOLOGY_CACHE_ALLOWANCE_GB * 1024 * 1024 * 1024
        ))
    fi
    if [[ "$LIST_ONLY" == "1" || "$DRY_RUN" == "1" ]]; then
        printf '%-18s %-32s %-16s %s\n' \
            homology "$HOMOLOGY_CACHE_ROLE" derived "$destination"
    fi
fi

(( selected_count > 0 )) || die "No catalogue entries selected"
unknown_allowance_bytes=$((missing_unknown_count * 1024 * 1024 * 1024))
planned_bytes=$((missing_known_bytes + unknown_allowance_bytes))

echo
echo "SAN frozen-input acquisition"
echo "  root:              $ROOT"
echo "  profiles:          ${PROFILES[*]}"
echo "  selected files:    $selected_count"
echo "  missing files:     $missing_count"
echo "  known bytes left:  $missing_known_bytes ($(human_gib "$missing_known_bytes"))"
echo "  unknown allowance: $unknown_allowance_bytes ($(human_gib "$unknown_allowance_bytes"))"
echo "  retained reserve:  ${RESERVE_GB} GiB"

if [[ "$LIST_ONLY" == "1" || "$DRY_RUN" == "1" ]]; then
    echo "No files were changed."
    exit 0
fi

if [[ "$ROOT" == "$DEFAULT_ROOT" ]]; then
    [[ -d "$ROOT" ]] || die "SAN root is not mounted or visible: $ROOT"
else
    mkdir -p "$ROOT"
fi
[[ -w "$ROOT" ]] || die "Store root is not writable: $ROOT"
mkdir -p "$ROOT/manifests"

if command -v flock >/dev/null 2>&1; then
    exec 9> "$ROOT/manifests/.populate_san_frozen_inputs.lock"
    flock -n 9 || die "Another SAN acquisition process is already running"
else
    FALLBACK_LOCK_DIR="$ROOT/manifests/.populate_san_frozen_inputs.lock.d"
    mkdir "$FALLBACK_LOCK_DIR" 2>/dev/null || \
        die "Another SAN acquisition process may be running: $FALLBACK_LOCK_DIR"
    trap release_lock EXIT INT TERM
fi

if [[ "$VERIFY_ONLY" == "0" && "$missing_count" -gt 0 ]]; then
    available_kib="$(df -Pk "$ROOT" | awk 'NR == 2 {print $4}')"
    available_bytes=$((available_kib * 1024))
    required_bytes=$((planned_bytes + RESERVE_GB * 1024 * 1024 * 1024))
    (( available_bytes >= required_bytes )) || die \
        "Insufficient SAN space: available=$(human_gib "$available_bytes") required_with_reserve=$(human_gib "$required_bytes")"
fi

UNIPROT_METADATA_BEFORE=""
GOA_METADATA_BEFORE=""
if [[ "$VERIFY_ONLY" == "0" && "$needs_uniprot_guard" == "1" ]]; then
    echo "Checking the mutable UniProt endpoint before acquisition"
    UNIPROT_METADATA_BEFORE="$(fetch_text 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt')"
    grep -Eq 'UniProt Release[[:space:]]+2026_02([^0-9]|$)' <<< "$UNIPROT_METADATA_BEFORE" || \
        die "UniProt current_release is no longer 2026_02; refusing to mislabel new data"
fi
if [[ "$VERIFY_ONLY" == "0" && "$needs_goa_guard" == "1" ]]; then
    echo "Checking the mutable GOA endpoint before acquisition"
    GOA_METADATA_BEFORE="$(fetch_text 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/current_release_numbers.txt')"
    awk '$1 == "uniprot" && $2 == "234" && $3 == "2026-06-17" {found=1} END {exit !found}' \
        <<< "$GOA_METADATA_BEFORE" || \
        die "GOA current release is no longer UniProt-GOA 234 from 2026-06-17"
fi

while IFS=$'\t' read -r profiles role release relative_path url expected_bytes \
    checksum_algorithm expected_checksum validator; do
    [[ -n "$profiles" ]] || continue
    [[ "$profiles" == \#* ]] && continue
    profile_selected "$profiles" || continue
    destination="$ROOT/$relative_path"
    provenance="${destination}.provenance.tsv"
    sidecar="${destination}.sha256"
    echo
    echo "[$role] $relative_path"

    if [[ -s "$destination" ]]; then
        verify_expected_size "$destination" "$expected_bytes"
        if [[ "$FULL_VERIFY" == "1" || ! -s "$provenance" || ! -s "$sidecar" ]]; then
            if [[ "$checksum_algorithm" == "-" && ! -s "$sidecar" ]]; then
                die "Existing unpinned file has no trusted SHA-256 sidecar: $destination"
            fi
            if [[ "$checksum_algorithm" == "-" ]]; then
                verify_recorded_sha256 "$destination"
            else
                verify_expected_checksum "$destination" "$checksum_algorithm" "$expected_checksum"
            fi
            validate_artifact "$destination" "$validator"
            VERIFIED=$((VERIFIED + 1))
        fi
        if [[ ! -s "$provenance" || ! -s "$sidecar" ]]; then
            write_artifact_metadata "$role" "$release" "$relative_path" "$url" \
                "$expected_bytes" "$checksum_algorithm" "$expected_checksum" "$destination" \
                "verified-existing-pinned"
        fi
        echo "  present and valid; skipping download"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    if [[ "$VERIFY_ONLY" == "1" ]]; then
        die "Selected artifact is missing: $destination"
    fi

    mkdir -p "$(dirname "$destination")"
    partial="${destination}.partial"
    echo "  downloading: $url"
    download_to_partial "$url" "$partial"
    verify_expected_size "$partial" "$expected_bytes"
    observed_sha="$(sha256_file "$partial")"
    if [[ "$checksum_algorithm" == "sha256" ]]; then
        [[ "$observed_sha" == "$expected_checksum" ]] || \
            die "sha256 mismatch for $partial: expected=$expected_checksum observed=$observed_sha"
    else
        verify_expected_checksum "$partial" "$checksum_algorithm" "$expected_checksum"
    fi
    validate_artifact "$partial" "$validator"
    mv "$partial" "$destination"
    write_artifact_metadata "$role" "$release" "$relative_path" "$url" \
        "$expected_bytes" "$checksum_algorithm" "$expected_checksum" "$destination" "downloaded" \
        "$observed_sha"
    echo "  authenticated and published atomically"
    DOWNLOADED=$((DOWNLOADED + 1))
    VERIFIED=$((VERIFIED + 1))
done < "$SPEC_FILE"

if [[ "$needs_uniprot_guard" == "1" ]]; then
    UNIPROT_METADATA_AFTER="$(fetch_text 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt')"
    [[ "$UNIPROT_METADATA_AFTER" == "$UNIPROT_METADATA_BEFORE" ]] || \
        die "UniProt current_release metadata changed during acquisition"
fi
if [[ "$needs_goa_guard" == "1" ]]; then
    GOA_METADATA_AFTER="$(fetch_text 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/current_release_numbers.txt')"
    [[ "$GOA_METADATA_AFTER" == "$GOA_METADATA_BEFORE" ]] || \
        die "GOA current release metadata changed during acquisition"
fi

process_temporal_derived_inputs
process_homology_derived_inputs

catalog_metadata
write_artifact_path_catalog

echo
echo "Completed successfully"
echo "  downloaded: $DOWNLOADED"
echo "  skipped:    $SKIPPED"
echo "  derived:    $DERIVED_CREATED created, $DERIVED_SKIPPED skipped"
echo "  full checks: $VERIFIED"
echo "  catalogue:  $ROOT/manifests/frozen_input_catalog.tsv"
echo "  path map:   $ROOT/manifests/artifact_paths.tsv"
du -sh "$ROOT"
