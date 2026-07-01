#!/bin/bash
# generate_embeddings_dependencies.sh
# Downloads ALL external resources for from-scratch embedding generation into
# PFP/external/, and writes external/dependency_env.sh with the env-var exports.
# Run from the PFP repo root.
set -euo pipefail

EXT="$(pwd)/external"
mkdir -p "${EXT}"
echo "==> External dependencies will live in: ${EXT}"

# --- 1. CAFA Assessment Tool (PPI, Text, Structure ID mapping) ---------
if [ ! -d "${EXT}/CAFA_assessment_tool" ]; then
  echo "==> Cloning CAFA_assessment_tool"
  git clone https://github.com/ashleyzhou972/CAFA_assessment_tool.git "${EXT}/CAFA_assessment_tool"
else
  echo "==> CAFA_assessment_tool already present, skipping"
fi

# --- 2. Raw CAFA3 CSVs — Zenodo 7409660 (verified bit-for-bit vs Zijan's splits).
#        Download 9 split CSVs, authenticate PRISTINE files against Zenodo md5s,
#        THEN normalise MF column 'protein' -> 'proteins' so prepare runs unmodified.
RAW="${EXT}/cafa3_raw"
mkdir -p "${RAW}"
CAFA3_BASE="https://zenodo.org/records/7409660/files"

cat > "${RAW}/.zenodo_md5.txt" <<'EOF'
e9a4b239cd47a7ac80975f63e259581e  bp-test.csv
85c19594547a503956226b9c225efc5d  bp-training.csv
c2674223770d6a8cf680dd9335d51ebe  bp-validation.csv
0e5dc8528ca95e8897b10cddaa12a775  cc-test.csv
074b13dd50fad4a6a4f13e4d8d4105d6  cc-training.csv
cdc8ceefcab4fb8c9278dd07c184327f  cc-validation.csv
2735e408dd57f6de29b1538f6b150d68  mf-test.csv
b31a8f22b5934aef61b76ec3b89296da  mf-training.csv
897921ce5df8174672200320926ccc87  mf-validation.csv
EOF

for aspect in bp cc mf; do
  for split in training validation test; do
    f="${RAW}/${aspect}-${split}.csv"
    if [ ! -f "$f" ]; then
      echo "==> Downloading ${aspect}-${split}.csv"
      wget -c "${CAFA3_BASE}/${aspect}-${split}.csv?download=1" -O "$f"
    fi
  done
done

if [ ! -f "${RAW}/.normalised" ]; then
  echo "==> Authenticating CAFA3 CSV md5s against Zenodo..."
  ( cd "${RAW}"
    if command -v md5sum >/dev/null 2>&1; then
      md5sum -c .zenodo_md5.txt
    else
      while read -r want name; do
        got=$(md5 -q "$name")
        [ "$got" = "$want" ] && echo "  OK  $name" || { echo "  BAD $name ($got != $want)"; exit 1; }
      done < .zenodo_md5.txt
    fi
  )
  echo "==> CAFA3 CSVs authenticated."
  for split in training validation test; do
    f="${RAW}/mf-${split}.csv"
    if head -1 "$f" | grep -q '^protein,'; then
      echo "==> Normalising header: mf-${split}.csv  protein -> proteins"
      if sed --version >/dev/null 2>&1; then
        sed -i '1 s/^protein,/proteins,/' "$f"
      else
        sed -i '' '1 s/^protein,/proteins,/' "$f"
      fi
    fi
  done
  touch "${RAW}/.normalised"
  echo "==> CAFA3 CSVs normalised."
else
  echo "==> CAFA3 CSVs already authenticated + normalised, skipping."
fi

# --- 3. STRING files (PPI): alias (confirmed URL) + network embeddings .h5 (manual)
mkdir -p "${EXT}/string"
STRING_ALIAS_GZ="${EXT}/string/protein.aliases.v12.0.txt.gz"
STRING_ALIAS="${EXT}/string/protein.aliases.v12.0.txt"
if [ ! -f "${STRING_ALIAS}" ]; then
  echo "==> Downloading STRING aliases v12.0 (~3.2 GB)"
  wget -c "https://stringdb-downloads.org/download/protein.aliases.v12.0.txt.gz" -O "${STRING_ALIAS_GZ}"
  gunzip "${STRING_ALIAS_GZ}"
else
  echo "==> STRING aliases already present, skipping"
fi

STRING_H5="${EXT}/string/protein.network.embeddings.v12.0.h5"
if [ ! -f "${STRING_H5}" ]; then
  echo "==> Downloading STRING network embeddings v12.0 (.h5) (~17.9 GB)"
  wget -c "https://stringdb-downloads.org/download/protein.network.embeddings.v12.0.h5" -O "${STRING_H5}"
else
  echo "==> STRING network embeddings already present, skipping"
fi


# --- 4. AlphaFold + UniProt are runtime API downloads, nothing to pre-fetch.
echo "==> AlphaFold & UniProt are runtime API downloads (handled in modality scripts)."

# --- 5. Emit env-var exports.
cat > "${EXT}/dependency_env.sh" <<EOF
# Source before running modality scripts:  source external/dependency_env.sh
export CAFA_ASSESSMENT_DIR="${EXT}/CAFA_assessment_tool"
export CAFA3_RAW_DIR="${RAW}"
export STRING_H5_FILE="${STRING_H5}"
export STRING_ALIAS_FILE="${STRING_ALIAS}"
EOF

echo ""
echo "==> Wrote env exports to ${EXT}/dependency_env.sh"
