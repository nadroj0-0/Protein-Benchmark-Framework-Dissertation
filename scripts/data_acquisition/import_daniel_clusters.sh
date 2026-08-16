#!/usr/bin/env bash
set -euo pipefail

BASE="/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02/supervisor_daniel_buchan"
INCOMING="$BASE/incoming/2026-08-06"
ARCHIVE="$INCOMING/cluster_assignments.tar.gz"
OUTPUT_ROOT="$BASE/mmseqs_cluster_assignments/uniref50_sensitivity_4"
EXPECTED_ARCHIVE_SHA256="014056109ed306c0d7125deffcdc2cdfd46a85a6a72cb1f901a8432e9dbdeaf2"
EXPECTED_ROWS=38794121

actual_archive_sha256="$(sha256sum "$ARCHIVE" | awk '{print $1}')"
[[ "$actual_archive_sha256" == "$EXPECTED_ARCHIVE_SHA256" ]] || {
  echo "Archive SHA-256 mismatch" >&2
  exit 1
}

for identity in 25 20 15; do
  member="threshold_${identity}/cluster_assignments_threshold_${identity}.tsv"
  destination="$OUTPUT_ROOT/identity_${identity}/raw"
  raw_partial="$destination/.cluster_assignments.tsv.partial"
  gzip_partial="$destination/.cluster_assignments.tsv.gz.partial"
  final="$destination/cluster_assignments.tsv.gz"

  mkdir -p "$destination"
  rm -f "$raw_partial" "$gzip_partial"
  tar -xOzf "$ARCHIVE" "$member" > "$raw_partial"

  awk -F '\t' -v expected="$EXPECTED_ROWS" '
    NF != 2 || $1 == "" || $2 == "" { bad += 1 }
    END {
      printf "rows=%d bad_rows=%d\n", NR, bad + 0
      if (NR != expected || bad != 0) exit 1
    }
  ' "$raw_partial" > "$destination/validation.txt.partial"

  raw_sha256="$(sha256sum "$raw_partial" | awk '{print $1}')"
  raw_size="$(stat -c '%s' "$raw_partial")"
  gzip -n -c "$raw_partial" > "$gzip_partial"
  gzip -t "$gzip_partial"
  mv "$gzip_partial" "$final"
  mv "$destination/validation.txt.partial" "$destination/validation.txt"

  compressed_sha256="$(sha256sum "$final" | awk '{print $1}')"
  compressed_size="$(stat -c '%s' "$final")"
  printf '%s  %s\n' "$compressed_sha256" "cluster_assignments.tsv.gz" \
    > "$destination/cluster_assignments.tsv.gz.sha256"
  printf '%s  %s\n' "$raw_sha256" "cluster_assignments_threshold_${identity}.tsv" \
    > "$destination/source_raw.tsv.sha256"

  {
    printf 'field\tvalue\n'
    printf 'producer\tDaniel Buchan\n'
    printf 'received_date\t2026-08-06\n'
    printf 'identity_percent\t%s\n' "$identity"
    printf 'source_archive\t%s\n' "$ARCHIVE"
    printf 'source_archive_sha256\t%s\n' "$actual_archive_sha256"
    printf 'source_archive_member\t%s\n' "$member"
    printf 'source_raw_size_bytes\t%s\n' "$raw_size"
    printf 'source_raw_sha256\t%s\n' "$raw_sha256"
    printf 'published_gzip_size_bytes\t%s\n' "$compressed_size"
    printf 'published_gzip_sha256\t%s\n' "$compressed_sha256"
    printf 'expected_members\t%s\n' "$EXPECTED_ROWS"
    printf 'validation\texact row count; two non-empty tab-separated columns; gzip integrity\n'
  } > "$destination/IMPORT_PROVENANCE.tsv"

  rm -f "$raw_partial"
  echo "Published identity ${identity}%: $final"
done
