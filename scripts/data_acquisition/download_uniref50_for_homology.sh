#!/usr/bin/env bash
# Stage the current UniProt 2026_02 UniRef50 FASTA without touching UniRef90.

set -euo pipefail

DESTINATION="${1:-}"
[[ -n "$DESTINATION" && "$DESTINATION" == /* && "$DESTINATION" != "/" ]] || {
    echo "Usage: $0 /absolute/path/uniref50.fasta.gz" >&2
    exit 2
}

URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta.gz"
RELNOTES_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt"
MINIMUM_FREE_GB="${MINIMUM_FREE_GB:-60}"
[[ "$MINIMUM_FREE_GB" =~ ^[1-9][0-9]*$ ]] || {
    echo "MINIMUM_FREE_GB must be a positive integer" >&2
    exit 2
}

mkdir -p "$(dirname "$DESTINATION")"
if [[ -s "$DESTINATION" && -s "$DESTINATION.sha256" ]]; then
    (cd "$(dirname "$DESTINATION")" && sha256sum -c "$(basename "$DESTINATION").sha256")
    gzip -t "$DESTINATION"
    echo "Validated existing UniRef50 input: $DESTINATION"
    exit 0
fi
[[ ! -e "$DESTINATION" ]] || {
    echo "Refusing an existing UniRef50 file without a valid SHA-256 sidecar: $DESTINATION" >&2
    exit 1
}

free_kb="$(df -Pk "$(dirname "$DESTINATION")" | awk 'NR==2 {print $4}')"
required_kb="$((MINIMUM_FREE_GB * 1024 * 1024))"
(( free_kb >= required_kb )) || {
    echo "Destination has ${free_kb} KiB free; ${required_kb} KiB is required" >&2
    exit 1
}

partial="$DESTINATION.partial"
relnotes_before="$DESTINATION.relnotes.before.partial"
relnotes_after="$DESTINATION.relnotes.after.partial"
checksum_partial="$DESTINATION.sha256.partial"
provenance_partial="$DESTINATION.provenance.tsv.partial"
cleanup() {
    local status=$?
    if [[ "$status" != "0" ]]; then
        rm -f "$relnotes_before" "$relnotes_after" "$checksum_partial" "$provenance_partial"
    fi
    exit "$status"
}
trap cleanup EXIT

curl --fail --location --retry 5 --output "$relnotes_before" "$RELNOTES_URL"
grep -Eq 'UniProt Release[[:space:]]+2026_02([^0-9]|$)' "$relnotes_before" || {
    echo "The mutable current_release endpoint is no longer UniProt release 2026_02" >&2
    exit 1
}

curl --fail --location --retry 5 --continue-at - --output "$partial" "$URL"
[[ -s "$partial" ]] || { echo "Downloaded UniRef50 file is empty" >&2; exit 1; }
gzip -t "$partial"
python3 - "$partial" <<'PY'
import gzip
import sys

with gzip.open(sys.argv[1], "rt", encoding="utf-8") as handle:
    first = handle.readline().strip()
if not first.startswith(">UniRef50_"):
    raise SystemExit(f"Unexpected UniRef50 FASTA header: {first!r}")
PY

curl --fail --location --retry 5 --output "$relnotes_after" "$RELNOTES_URL"
cmp -s "$relnotes_before" "$relnotes_after" || {
    echo "UniProt current_release changed during the download" >&2
    exit 1
}

digest="$(sha256sum "$partial" | awk '{print $1}')"
printf '%s  %s\n' "$digest" "$(basename "$DESTINATION")" > "$checksum_partial"
printf 'role\trelease\turl\tsha256\tsize_bytes\n' > "$provenance_partial"
printf 'uniref50_t1\t2026_02\t%s\t%s\t%s\n' \
    "$URL" "$digest" "$(stat -c '%s' "$partial")" >> "$provenance_partial"

mv "$partial" "$DESTINATION"
mv "$checksum_partial" "$DESTINATION.sha256"
mv "$provenance_partial" "$DESTINATION.provenance.tsv"
mv "$relnotes_after" "$DESTINATION.relnotes.txt"
rm -f "$relnotes_before"
echo "Published validated UniRef50 input: $DESTINATION"
echo "SHA-256: $digest"
