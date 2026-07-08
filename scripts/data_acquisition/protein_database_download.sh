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

set -euo pipefail

ROOT="$HOME/protein_databases"
TODAY="$(date +%F)"

UNIPROT_T0_RELEASE="2025_01"
UNIPROT_T1_RELEASE="2026_02"

GOA_T0_RELEASE="225"
GOA_T1_RELEASE="234"

UNIPROT_T0_DIR="$ROOT/uniprot/release_${UNIPROT_T0_RELEASE}"
UNIPROT_T1_DIR="$ROOT/uniprot/release_${UNIPROT_T1_RELEASE}"
GOA_T0_DIR="$ROOT/goa/release_${UNIPROT_T0_RELEASE}"
GOA_T1_DIR="$ROOT/goa/release_${UNIPROT_T1_RELEASE}"
ONTOLOGY_DIR="$ROOT/ontology"

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

echo "=============================================================="
echo "Protein Database Download Utility"
echo "=============================================================="

mkdir -p \
    "$UNIPROT_T0_DIR" \
    "$UNIPROT_T1_DIR" \
    "$GOA_T0_DIR" \
    "$GOA_T1_DIR" \
    "$ONTOLOGY_DIR" \
    "$ROOT/string" \
    "$ROOT/alphafold"

echo
echo "[1/7] UniProt ${UNIPROT_T0_RELEASE}"

download_if_missing \
"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-${UNIPROT_T0_RELEASE}/knowledgebase/uniprot_sprot-only${UNIPROT_T0_RELEASE}.tar.gz" \
"$UNIPROT_T0_DIR/uniprot_sprot-only${UNIPROT_T0_RELEASE}.tar.gz"

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
echo "[2/7] UniProt ${UNIPROT_T1_RELEASE}"

UNIPROT_CURRENT_BASE="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"

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

download_if_missing \
"https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt" \
"$UNIPROT_T1_DIR/relnotes.txt"

download_if_missing \
"https://ftp.uniprot.org/pub/databases/uniprot/current_release/changes.html" \
"$UNIPROT_T1_DIR/changes.html"

echo
echo "[3/7] GOA ${GOA_T0_RELEASE}"

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.${GOA_T0_RELEASE}.gz" \
"$GOA_T0_DIR/goa_uniprot_all.gaf.${GOA_T0_RELEASE}.gz"

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/README" \
"$GOA_T0_DIR/README"

echo
echo "[4/7] GOA ${GOA_T1_RELEASE}"

#download_if_missing \
#"https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz" \
#"$GOA_T1_DIR/goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz"

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz.md5" \
"$GOA_T1_DIR/goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz.md5"

download_if_missing \
"https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README" \
"$GOA_T1_DIR/README"

echo
echo "[5/7] GO ontology"

download_if_missing \
"http://current.geneontology.org/ontology/go-basic.obo" \
"$ONTOLOGY_DIR/go-basic.obo"

download_if_missing \
"http://current.geneontology.org/ontology/go.obo" \
"$ONTOLOGY_DIR/go.obo"

echo
echo "[6/7] Writing manifest"

cat > "$ROOT/MANIFEST.md" <<EOF
# Protein Database Manifest

Generated: ${TODAY}

## Purpose

Reference biological databases for the contemporary CAFA-style temporal benchmark.

## Training snapshot / t0

### UniProt

Release: ${UNIPROT_T0_RELEASE}

Directory:

\`\`\`
$UNIPROT_T0_DIR
\`\`\`

Primary archive:

\`\`\`
uniprot_sprot-only${UNIPROT_T0_RELEASE}.tar.gz
\`\`\`

Source:

\`\`\`
https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-${UNIPROT_T0_RELEASE}/
\`\`\`

### GOA

GOA release: ${GOA_T0_RELEASE}

Directory:

\`\`\`
$GOA_T0_DIR
\`\`\`

Primary annotation file:

\`\`\`
goa_uniprot_all.gaf.${GOA_T0_RELEASE}.gz
\`\`\`

Source:

\`\`\`
https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/
\`\`\`

## Evaluation snapshot / t1

### UniProt

Release: ${UNIPROT_T1_RELEASE}

Directory:

\`\`\`
$UNIPROT_T1_DIR
\`\`\`

Files:

\`\`\`
uniprot_sprot.dat.gz
uniprot_sprot.fasta.gz
uniprot_sprot.xml.gz
uniprot_sprot_varsplic.fasta.gz
\`\`\`

Source:

\`\`\`
https://ftp.uniprot.org/pub/databases/uniprot/current_release/
\`\`\`

### GOA

GOA release: ${GOA_T1_RELEASE}

Directory:

\`\`\`
$GOA_T1_DIR
\`\`\`

Primary annotation file:

\`\`\`
goa_uniprot_all.gaf.${GOA_T1_RELEASE}.gz
\`\`\`

Source:

\`\`\`
https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/
\`\`\`

## Ontology

Directory:

\`\`\`
$ONTOLOGY_DIR
\`\`\`

Files:

\`\`\`
go-basic.obo
go.obo
\`\`\`

Source:

\`\`\`
http://current.geneontology.org/ontology/
\`\`\`

## Notes

The GOA 2026 file was downloaded from the current GOA directory and saved locally with release number ${GOA_T1_RELEASE} to freeze the snapshot used in this benchmark.

Large GOA GAF files should normally be streamed directly from gzip rather than decompressed, to avoid disk quota problems.
EOF

echo
echo "[7/7] Summary"

find "$ROOT" -type f | sort

echo
echo "Disk usage:"
du -sh "$ROOT"

echo
echo "Finished successfully."
