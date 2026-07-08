#!/bin/bash
set -euo pipefail

ROOT="${1:-$HOME/protein_databases}"
FULL_STATS="${FULL_STATS:-false}"
VERIFY_GZIP="${VERIFY_GZIP:-false}"
SAMPLE_N="${SAMPLE_N:-100000}"

line() { printf '%*s\n' 80 '' | tr ' ' '='; }
section() { echo; line; echo "$1"; line; }
exists() { [[ -e "$1" ]]; }

inspect_file() {
    local f="$1"
    if exists "$f"; then
        echo "✓ $f"
        ls -lh "$f"
    else
        echo "✗ missing: $f"
    fi
}

preview_gzip() {
    local f="$1"
    [[ ! -f "$f" ]] && return
    echo "quick gzip readability check:"
    if zcat "$f" | head -5 >/dev/null 2>&1; then
        echo "  OK"
    else
        echo "  FAILED"
    fi
    echo "first non-comment lines:"
    zgrep -m 5 -v '^!' "$f" 2>/dev/null || true
}

inspect_gzip() {
    local f="$1"
    inspect_file "$f"
    preview_gzip "$f"
}

verify_gzip() {
    local f="$1"
    inspect_file "$f"
    [[ ! -f "$f" ]] && return
    if [[ "$VERIFY_GZIP" == "true" ]]; then
        echo "full gzip integrity check:"
        if gzip -t "$f" 2>/tmp/gzip_check.err; then
            echo "  OK"
        else
            echo "  FAILED"
            cat /tmp/gzip_check.err
        fi
    else
        echo "full gzip integrity check: skipped (set VERIFY_GZIP=true to enable)"
    fi
    preview_gzip "$f"
}

inspect_tar() {
    local f="$1"
    inspect_file "$f"
    if [[ -f "$f" ]]; then
        echo "tar contents preview:"
        tar -tzf "$f" | head -20 || true
    fi
}

count_fasta_sequences() {
    local f="$1"
    [[ ! -f "$f" ]] && return
    echo
    echo "FASTA sequence count:"
    zgrep -c '^>' "$f" || true
    echo "First FASTA header:"
    zgrep -m 1 '^>' "$f" || true
    echo "Last FASTA header: skipped in quick mode"
}

count_uniprot_dat_entries() {
    local f="$1"
    [[ ! -f "$f" ]] && return
    echo
    echo "DAT entry count:"
    zgrep -c '^ID   ' "$f" || true
    echo "First ID lines:"
    zgrep -m 5 '^ID   ' "$f" || true
}

count_uniprot_xml_entries() {
    local f="$1"
    [[ ! -f "$f" ]] && return
    echo
    echo "XML entry count:"
    zgrep -c '<entry ' "$f" || true
    echo "First entry line:"
    zgrep -m 1 '<entry ' "$f" || true
}

inspect_goa_sample() {
    local f="$1"
    [[ ! -f "$f" ]] && return
    echo
    echo "GAF header/comment preview:"
    zcat "$f" | head -20 || true
    echo
    echo "First data rows:"
    zgrep -m 5 -v '^!' "$f" || true
    echo
    echo "Column count check on first 1000 data rows:"
    zgrep -m 1000 -v '^!' "$f" | awk -F'\t' '{c[NF]++} END{for (n in c) print n, c[n]}' | sort -n || true
    echo
    echo "Evidence-code counts in first ${SAMPLE_N} data rows:"
    zgrep -m "$SAMPLE_N" -v '^!' "$f" | awk -F'\t' '{c[$7]++} END{for (e in c) print e, c[e]}' | sort || true
    echo
    echo "Aspect counts in first ${SAMPLE_N} data rows:"
    zgrep -m "$SAMPLE_N" -v '^!' "$f" | awk -F'\t' '{c[$9]++} END{for (a in c) print a, c[a]}' | sort || true
    echo
    echo "Date range in first ${SAMPLE_N} data rows:"
    zgrep -m "$SAMPLE_N" -v '^!' "$f" | awk -F'\t' 'NR==1{min=$14;max=$14} {if($14<min)min=$14;if($14>max)max=$14} END{print min, max}' || true
    echo
    echo "Example CAFA target taxa rows, first hits:"
    for tax in 9606 10090 10116 3702 7955 8355 44689 559292 284812 83333 7227; do
        echo "taxon:$tax"
        zgrep -m 2 "taxon:$tax" "$f" || true
    done
}

inspect_goa_full_stats() {
    local f="$1"
    [[ ! -f "$f" ]] && return
    if [[ "$FULL_STATS" != "true" ]]; then
        echo
        echo "Full GOA statistics skipped (run FULL_STATS=true $0 to enable)."
        return
    fi
    echo
    echo "FULL GOA STATISTICS"
    echo "Total annotations:"
    zgrep -vc '^!' "$f" || true
    echo
    echo "Unique proteins:"
    zgrep -v '^!' "$f" | cut -f2 | sort -u | wc -l || true
    echo
    echo "Unique GO terms:"
    zgrep -v '^!' "$f" | cut -f5 | sort -u | wc -l || true
    echo
    echo "Unique taxa:"
    zgrep -v '^!' "$f" | cut -f13 | sort -u | wc -l || true
    echo
    echo "Evidence-code counts (whole file):"
    zgrep -v '^!' "$f" | cut -f7 | sort | uniq -c | sort -nr || true
    echo
    echo "Aspect counts (whole file):"
    zgrep -v '^!' "$f" | cut -f9 | sort | uniq -c | sort -nr || true
    echo
    echo "Most common GO terms (top 20):"
    zgrep -v '^!' "$f" | cut -f5 | sort | uniq -c | sort -nr | head -20 || true
    echo
    echo "Most common taxa (top 20):"
    zgrep -v '^!' "$f" | cut -f13 | sort | uniq -c | sort -nr | head -20 || true
    echo
    echo "Annotation date range (whole file):"
    zgrep -v '^!' "$f" | cut -f14 | sort | awk 'NR==1{min=$1} {max=$1} END{print min, max}' || true
}

inspect_ontology() {
    local f="$1"
    inspect_file "$f"
    [[ ! -f "$f" ]] && return
    echo
    echo "Preview: $f"
    head -40 "$f" || true
    echo
    echo "Term count:"
    grep -c '^\[Term\]' "$f" || true
    echo "is_a count:"
    grep -c '^is_a:' "$f" || true
    echo "part_of count:"
    grep -c 'relationship: part_of' "$f" || true
}

section "ROOT"
echo "$ROOT"
if [[ ! -d "$ROOT" ]]; then
    echo "Root directory does not exist."
    exit 0
fi

section "DIRECTORY TREE"
tree -L 4 "$ROOT" 2>/dev/null || find "$ROOT" -maxdepth 4 -print

section "DISK USAGE"
du -sh "$ROOT" || true
du -sh "$ROOT"/* 2>/dev/null || true

section "UNIPROT 2025_01"
U25="$ROOT/uniprot/release_2025_01"
inspect_tar "$U25/uniprot_sprot-only2025_01.tar.gz"
inspect_file "$U25/relnotes.txt"
inspect_file "$U25/changes.html"
inspect_file "$U25/news.html"
if [[ -d "$U25/extracted" ]]; then
    echo
    echo "Extracted files:"
    find "$U25/extracted" -maxdepth 2 -type f -exec ls -lh {} \;
    echo
    inspect_gzip "$U25/extracted/uniprot_sprot.fasta.gz"
    count_fasta_sequences "$U25/extracted/uniprot_sprot.fasta.gz"
    inspect_gzip "$U25/extracted/uniprot_sprot.dat.gz"
    count_uniprot_dat_entries "$U25/extracted/uniprot_sprot.dat.gz"
    inspect_gzip "$U25/extracted/uniprot_sprot.xml.gz"
    count_uniprot_xml_entries "$U25/extracted/uniprot_sprot.xml.gz"
    inspect_gzip "$U25/extracted/uniprot_sprot_varsplic.fasta.gz"
    count_fasta_sequences "$U25/extracted/uniprot_sprot_varsplic.fasta.gz"
else
    echo "No extracted directory found."
fi

section "UNIPROT 2026_02"
U26="$ROOT/uniprot/release_2026_02"
inspect_gzip "$U26/uniprot_sprot.fasta.gz"
count_fasta_sequences "$U26/uniprot_sprot.fasta.gz"
inspect_gzip "$U26/uniprot_sprot.dat.gz"
count_uniprot_dat_entries "$U26/uniprot_sprot.dat.gz"
inspect_gzip "$U26/uniprot_sprot.xml.gz"
count_uniprot_xml_entries "$U26/uniprot_sprot.xml.gz"
inspect_gzip "$U26/uniprot_sprot_varsplic.fasta.gz"
count_fasta_sequences "$U26/uniprot_sprot_varsplic.fasta.gz"
inspect_file "$U26/relnotes.txt"
inspect_file "$U26/changes.html"

section "UNIPROT 2026_02 CONSISTENCY CHECK"
if [[ -f "$U26/uniprot_sprot.fasta.gz" && -f "$U26/uniprot_sprot.dat.gz" && -f "$U26/uniprot_sprot.xml.gz" ]]; then
    FASTA_COUNT=$(zgrep -c '^>' "$U26/uniprot_sprot.fasta.gz" || echo 0)
    DAT_COUNT=$(zgrep -c '^ID   ' "$U26/uniprot_sprot.dat.gz" || echo 0)
    XML_COUNT=$(zgrep -c '<entry ' "$U26/uniprot_sprot.xml.gz" || echo 0)
    echo "FASTA entries: $FASTA_COUNT"
    echo "DAT entries:   $DAT_COUNT"
    echo "XML entries:   $XML_COUNT"
    if [[ "$FASTA_COUNT" == "$DAT_COUNT" && "$DAT_COUNT" == "$XML_COUNT" ]]; then
        echo "✓ All core UniProt 2026_02 formats have matching entry counts."
    else
        echo "⚠ UniProt 2026_02 format counts differ."
    fi
else
    echo "Skipped: missing one or more core UniProt 2026_02 files."
fi

section "GOA 2025_01 / RELEASE 225"
G25="$ROOT/goa/release_2025_01"
inspect_gzip "$G25/goa_uniprot_all.gaf.225.gz"
inspect_file "$G25/README"
inspect_goa_sample "$G25/goa_uniprot_all.gaf.225.gz"
inspect_goa_full_stats "$G25/goa_uniprot_all.gaf.225.gz"

section "GOA 2026_02 / RELEASE 234"
G26="$ROOT/goa/release_2026_02"
inspect_gzip "$G26/goa_uniprot_all.gaf.234.gz"
inspect_file "$G26/goa_uniprot_all.gaf.234.gz.md5"
inspect_file "$G26/README"
inspect_goa_sample "$G26/goa_uniprot_all.gaf.234.gz"
inspect_goa_full_stats "$G26/goa_uniprot_all.gaf.234.gz"

section "ONTOLOGY"
O="$ROOT/ontology"
inspect_ontology "$O/go-basic.obo"
inspect_ontology "$O/go.obo"

section "MANIFEST"
inspect_file "$ROOT/MANIFEST.md"
if [[ -f "$ROOT/MANIFEST.md" ]]; then
    cat "$ROOT/MANIFEST.md"
fi

section "SUMMARY"
echo "Completed inspection."
echo "Quick mode used by default. For expensive whole-file GOA statistics run:"
echo "FULL_STATS=true ./inspect_protein_databases.sh"
echo "For full gzip verification run:"
echo "VERIFY_GZIP=true ./inspect_protein_databases.sh"

