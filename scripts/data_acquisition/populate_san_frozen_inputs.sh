#!/bin/bash
# Populate the dissertation SAN with authenticated, frozen public inputs.
#
# This is deliberately separate from protein_database_download.sh. That older
# script documents the historical home-directory acquisition workflow; this
# script owns the persistent /SAN/bioinf/bmpfp contract.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPEC_FILE="${SAN_INPUT_SPEC:-${SCRIPT_DIR}/san_frozen_inputs.tsv}"
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
FALLBACK_LOCK_DIR=""

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
                        downloaded files are always fully verified.
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

catalog_metadata

echo
echo "Completed successfully"
echo "  downloaded: $DOWNLOADED"
echo "  skipped:    $SKIPPED"
echo "  full checks: $VERIFIED"
echo "  catalogue:  $ROOT/manifests/frozen_input_catalog.tsv"
du -sh "$ROOT"
