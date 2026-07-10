#!/bin/bash
#
# protein_database_download.sh
#
# Downloads and organises the reference biological databases used throughout
# the Protein Benchmark Framework dissertation.
#
# Current downloads:
#   - UniProt Swiss-Prot release 2025_01 (t0)
#   - GOA release 225 (corresponding to UniProt 2025_01)
#
# The directory layout is:
#
#   ~/protein_databases/
#       uniprot/
#           release_2025_01/
#           release_2026_02/
#       goa/
#           release_2025_01/
#           release_2026_02/
#       ontology/
#       string/
#       alphafold/
#
# The script is safe to rerun:
#   - creates directories if absent
#   - resumes interrupted downloads using wget -c
#
# Full target-organism coverage also needs TrEMBL. Because those archives are
# very large, opt in explicitly with DOWNLOAD_TREMBL=1.
#

set -euo pipefail

ROOT="$HOME/protein_databases"
TODAY="$(date +%F)"
DOWNLOAD_TREMBL="${DOWNLOAD_TREMBL:-0}"

UNIPROT_T0_RELEASE="2025_01"
UNIPROT_T1_RELEASE="2026_02"

GOA_T0_RELEASE="225"
GOA_T1_RELEASE="234"

# Public GO product releases used by the benchmark builder. GOA 225 declares a
# 2025-03-07 ontology that is not retained as a standalone GO release, so the
# prediction graph is frozen to the last preceding release (2025-02-06) and the
# first following release (2025-03-16) is used only for source-ID resolution.
# The 2026-06-19 product contains data-version releases/2026-06-15, exactly
# matching GOA 234.
ONTOLOGY_T0_BENCHMARK_RELEASE="2025-02-06"
ONTOLOGY_T0_SOURCE_RELEASE="2025-03-16"
ONTOLOGY_T1_SOURCE_RELEASE="2026-06-19"

UNIPROT_T0_DIR="$ROOT/uniprot/release_${UNIPROT_T0_RELEASE}"
UNIPROT_T1_DIR="$ROOT/uniprot/release_${UNIPROT_T1_RELEASE}"
GOA_T0_DIR="$ROOT/goa/release_${UNIPROT_T0_RELEASE}"
GOA_T1_DIR="$ROOT/goa/release_${UNIPROT_T1_RELEASE}"
ONTOLOGY_T0_DIR="$ROOT/ontology/release_${ONTOLOGY_T0_BENCHMARK_RELEASE}"
ONTOLOGY_T0_SOURCE_DIR="$ROOT/ontology/release_${ONTOLOGY_T0_SOURCE_RELEASE}"
ONTOLOGY_T1_DIR="$ROOT/ontology/release_${ONTOLOGY_T1_SOURCE_RELEASE}"

download_if_missing() {
    local url="$1"
    local out="$2"

    if [[ -f "$out" ]]; then
        echo "✓ $(basename "$out") already exists - skipping"
    else
        echo "Downloading $(basename "$out")..."
        wget -c "$url" -O "$out"
    fi
}

assert_file_contains() {
    local file="$1"
    local expected="$2"
    if ! grep -q "$expected" "$file"; then
        echo "Release check failed: $file does not contain $expected" >&2
        echo "A mutable current endpoint may have advanced; use the archived release instead." >&2
        exit 1
    fi
}

echo "=============================================================="
echo "Protein Database Download Utility"
echo "=============================================================="

mkdir -p \
    "$UNIPROT_T0_DIR" \
    "$UNIPROT_T1_DIR" \
    "$GOA_T0_DIR" \
    "$GOA_T1_DIR" \
    "$ONTOLOGY_T0_DIR" \
    "$ONTOLOGY_T0_SOURCE_DIR" \
    "$ONTOLOGY_T1_DIR" \
    "$ROOT/string" \
    "$ROOT/alphafold"

echo
echo "[1/8] UniProt ${UNIPROT_T0_RELEASE}"

download_if_missing \
"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-${UNIPROT_T0_RELEASE}/knowledgebase/uniprot_sprot-only${UNIPROT_T0_RELEASE}.tar.gz" \
"$UNIPROT_T0_DIR/uniprot_sprot-only${UNIPROT_T0_RELEASE}.tar.gz"

if [[ "$DOWNLOAD_TREMBL" == "1" ]]; then
    download_if_missing \
    "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-${UNIPROT_T0_RELEASE}/knowledgebase/knowledgebase${UNIPROT_T0_RELEASE}.tar.gz" \
    "$UNIPROT_T0_DIR/knowledgebase${UNIPROT_T0_RELEASE}.tar.gz"
fi

download_if_missing \
"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-${UNIPROT_T0_RELEASE}/relnotes.txt" \
"$UNIPROT_T0_DIR/relnotes.txt"

download_if_missing \
"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-${UNIPROT_T0_RELEASE}/changes.html" \
"$UNIPROT_T0_DIR/changes.html"

download_if_missing \
"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-${UNIPROT_T0_RELEASE}/news.html" \
"$UNIPROT_T0_DIR/news.html"

echo
echo "[2/8] UniProt ${UNIPROT_T1_RELEASE}"

UNIPROT_CURRENT_BASE="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"

download_if_missing \
"https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt" \
"$UNIPROT_T1_DIR/relnotes.txt"

assert_file_contains "$UNIPROT_T1_DIR/relnotes.txt" "$UNIPROT_T1_RELEASE"

for f in \
    uniprot_sprot.dat.gz \
    uniprot_sprot.fasta.gz \
    uniprot_sprot.xml.gz \
    uniprot_sprot_varsplic.fasta.gz
do
    download_if_missing \
    "$UNIPROT_CURRENT_BASE/$f" \
    "$UNIPROT_T1_DIR/$f"
done

if [[ "$DOWNLOAD_TREMBL" == "1" ]]; then
    download_if_missing \
    "$UNIPROT_CURRENT_BASE/uniprot_trembl.dat.gz" \
    "$UNIPROT_T1_DIR/uniprot_trembl.dat.gz"
fi

download_if_missing \
"https://ftp.uniprot.org/pub/databases/uniprot/current_release/changes.html" \
"$UNIPROT_T1_DIR/changes.html"

echo
echo "[3/8] GOA ${GOA_T0_RELEASE}"

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.${GOA_T0_RELEASE}.gz" \
"$GOA_T0_DIR/goa_uniprot_all.gaf.${GOA_T0_RELEASE}.gz"

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/README" \
"$GOA_T0_DIR/README"

echo
echo "[4/8] GOA ${GOA_T1_RELEASE}"

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/current_release_numbers.txt" \
"$GOA_T1_DIR/current_release_numbers.txt"

assert_file_contains "$GOA_T1_DIR/current_release_numbers.txt" "$GOA_T1_RELEASE"

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz" \
"$GOA_T1_DIR/goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz"

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz.md5" \
"$GOA_T1_DIR/goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz.md5"

EXPECTED_GOA_T1_MD5="$(awk '{print $1}' "$GOA_T1_DIR/goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz.md5")"
ACTUAL_GOA_T1_MD5="$(md5sum "$GOA_T1_DIR/goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz" | awk '{print $1}')"
if [[ "$EXPECTED_GOA_T1_MD5" != "$ACTUAL_GOA_T1_MD5" ]]; then
    echo "GOA ${GOA_T1_RELEASE} checksum mismatch" >&2
    exit 1
fi
if ! zgrep -m 1 '^!date-generated:.*2026-06-17' \
    "$GOA_T1_DIR/goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz" >/dev/null; then
    echo "Downloaded GOA file is not the frozen release 234 snapshot from 2026-06-17" >&2
    exit 1
fi
if ! zgrep -m 1 '^!go-version:.*2026-06-15' \
    "$GOA_T1_DIR/goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz" >/dev/null; then
    echo "Downloaded GOA file does not declare GO 2026-06-15" >&2
    exit 1
fi

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README" \
"$GOA_T1_DIR/README"

echo
echo "[5/8] GO ontologies ${ONTOLOGY_T0_BENCHMARK_RELEASE} and ${ONTOLOGY_T0_SOURCE_RELEASE}"

download_if_missing \
"https://release.geneontology.org/${ONTOLOGY_T0_BENCHMARK_RELEASE}/ontology/go-basic.obo" \
"$ONTOLOGY_T0_DIR/go-basic.obo"

download_if_missing \
"https://release.geneontology.org/${ONTOLOGY_T0_BENCHMARK_RELEASE}/ontology/go.obo" \
"$ONTOLOGY_T0_DIR/go.obo"

download_if_missing \
"https://release.geneontology.org/${ONTOLOGY_T0_BENCHMARK_RELEASE}/summary.txt" \
"$ONTOLOGY_T0_DIR/summary.txt"

download_if_missing \
"https://release.geneontology.org/${ONTOLOGY_T0_SOURCE_RELEASE}/ontology/go-basic.obo" \
"$ONTOLOGY_T0_SOURCE_DIR/go-basic.obo"

echo
echo "[6/8] GO ontology ${ONTOLOGY_T1_SOURCE_RELEASE}"

download_if_missing \
"https://release.geneontology.org/${ONTOLOGY_T1_SOURCE_RELEASE}/ontology/go-basic.obo" \
"$ONTOLOGY_T1_DIR/go-basic.obo"

download_if_missing \
"https://release.geneontology.org/${ONTOLOGY_T1_SOURCE_RELEASE}/ontology/go.obo" \
"$ONTOLOGY_T1_DIR/go.obo"

download_if_missing \
"https://release.geneontology.org/${ONTOLOGY_T1_SOURCE_RELEASE}/summary.txt" \
"$ONTOLOGY_T1_DIR/summary.txt"

echo
echo "[7/8] Writing manifest"

cat > "$ROOT/MANIFEST.md" <<EOF
# Protein Database Manifest

Generated: ${TODAY}

## Purpose

Reference biological databases for the contemporary CAFA-style temporal benchmark.

---

# Training snapshot (t0)

## UniProt

Release: ${UNIPROT_T0_RELEASE}

Directory:

\`\`\`
${UNIPROT_T0_DIR}
\`\`\`

Primary archive:

\`\`\`
uniprot_sprot-only${UNIPROT_T0_RELEASE}.tar.gz
\`\`\`

Source:

\`\`\`
https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-${UNIPROT_T0_RELEASE}/
\`\`\`

## GOA

Release: ${GOA_T0_RELEASE}

Directory:

\`\`\`
${GOA_T0_DIR}
\`\`\`

Primary annotation file:

\`\`\`
goa_uniprot_all.gaf.${GOA_T0_RELEASE}.gz
\`\`\`

Source:

\`\`\`
https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/
\`\`\`

## Gene Ontology

Release: ${ONTOLOGY_T0_BENCHMARK_RELEASE}

Directory:

\`\`\`
${ONTOLOGY_T0_DIR}
\`\`\`

Files:

\`\`\`
go-basic.obo
go.obo
\`\`\`

Source:

\`\`\`
https://release.geneontology.org/${ONTOLOGY_T0_BENCHMARK_RELEASE}/ontology/
\`\`\`

---

# Evaluation snapshot (t1)

## UniProt

Release: ${UNIPROT_T1_RELEASE}

Directory:

\`\`\`
${UNIPROT_T1_DIR}
\`\`\`

Files:

\`\`\`
uniprot_sprot.dat.gz
uniprot_sprot.fasta.gz
uniprot_sprot.xml.gz
uniprot_sprot_varsplic.fasta.gz
relnotes.txt
changes.html
\`\`\`

Source:

\`\`\`
https://ftp.uniprot.org/pub/databases/uniprot/current_release/
\`\`\`

## GOA

Current release: ${GOA_T1_RELEASE}

Directory:

\`\`\`
${GOA_T1_DIR}
\`\`\`

Files:

\`\`\`
current_release_numbers.txt
goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz
goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz.md5
README
\`\`\`

Source:

\`\`\`
https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/
\`\`\`

## Gene Ontology

Current release

Directory:

\`\`\`
${ONTOLOGY_T1_DIR}
\`\`\`

Files:

\`\`\`
go-basic.obo
go.obo
\`\`\`

Source:

\`\`\`
https://current.geneontology.org/ontology/
\`\`\`

---

## Notes

- Training (t0) uses archived, fixed snapshots to ensure reproducibility.
- Evaluation (t1) uses the current releases available at download time.
- GOA current release metadata is recorded in \`current_release_numbers.txt\`.
- Large GOA GAF files should be streamed directly from the compressed archive rather than decompressed to disk where possible.
EOF

echo
echo "[8/8] Summary"

find "$ROOT" -type f | sort

echo
echo "Disk usage:"
du -sh "$ROOT"

echo
echo "Finished successfully."
