#!/bin/bash
# generate_embeddings_dependencies.sh
# Downloads ALL external resources for from-scratch embedding generation into
# PFP/external/, and writes external/dependency_env.sh with the env-var exports.
# Run from the PFP repo root.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
if [ -f "${REPO_ROOT}/configs/paths.local.sh" ]; then
  # Machine-specific paths are intentionally not committed.
  # shellcheck disable=SC1091
  source "${REPO_ROOT}/configs/paths.local.sh"
fi
# shellcheck source=../reproduction_common.sh
source "${REPO_ROOT}/scripts/reproduction_common.sh"
artifact_catalog_configure "${REPO_ROOT}" "${ARTIFACT_CATALOG:-}"

PFP_ROOT="${PFP_ROOT:-$(pwd)}"
EXT="${PFP_EXTERNAL_DIR:-${PFP_ROOT}/external}"
DATA_DIR="${PFP_DATA_DIR:-${PFP_ROOT}/data}"
DEPENDENCY_ENV="${DEPENDENCY_ENV:-${EXT}/dependency_env.sh}"
CAFA_ASSESSMENT_REPO_URL="${CAFA_ASSESSMENT_REPO_URL:-https://github.com/ashleyzhou972/CAFA_assessment_tool.git}"
CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-${EXT}/CAFA_assessment_tool}"
CAFA_ASSESSMENT_COMMIT="${CAFA_ASSESSMENT_COMMIT:-}"
PFP_CAFA3_RAW_DIR="${PFP_CAFA3_RAW_DIR:-${EXT}/cafa3_raw}"
PFP_STRING_DIR="${PFP_STRING_DIR:-${EXT}/string}"
CAFA3_BASE="${CAFA3_BASE:-https://zenodo.org/records/7409660/files}"
STRING_DOWNLOAD_BASE="${STRING_DOWNLOAD_BASE:-https://stringdb-downloads.org/download}"
CAFA3_SOURCE_DIR="${CAFA3_SOURCE_DIR:-${SAN_CAFA3_RAW_DIR:-}}"
EMBEDDING_DEPENDENCY_PROFILE="${EMBEDDING_DEPENDENCY_PROFILE:-all}"
case "${EMBEDDING_DEPENDENCY_PROFILE}" in
  all|sequence|text|structure|ppi) ;;
  *)
    echo "Unknown EMBEDDING_DEPENDENCY_PROFILE: ${EMBEDDING_DEPENDENCY_PROFILE}" >&2
    exit 2
    ;;
esac

mkdir -p "${EXT}"
echo "==> External dependencies will live in: ${EXT}"
echo "==> PFP data directory: ${DATA_DIR}"

git_in_dir() {
  local directory="$1"
  shift
  (cd "$directory" && git "$@")
}

# --- 1. CAFA Assessment Tool (PPI, Text, Structure ID mapping) ---------
if [ ! -d "${CAFA_ASSESSMENT_DIR}" ]; then
  echo "==> Cloning CAFA_assessment_tool"
  mkdir -p "$(dirname "${CAFA_ASSESSMENT_DIR}")"
  git clone "${CAFA_ASSESSMENT_REPO_URL}" "${CAFA_ASSESSMENT_DIR}"
else
  echo "==> CAFA_assessment_tool already present, skipping"
fi
if [ -n "${CAFA_ASSESSMENT_COMMIT}" ]; then
  [ -d "${CAFA_ASSESSMENT_DIR}/.git" ] || {
    echo "Cannot pin non-Git CAFA assessment directory: ${CAFA_ASSESSMENT_DIR}" >&2
    exit 1
  }
  git_in_dir "${CAFA_ASSESSMENT_DIR}" checkout --detach "${CAFA_ASSESSMENT_COMMIT}"
  observed_cafa_commit="$(git_in_dir "${CAFA_ASSESSMENT_DIR}" rev-parse HEAD)"
  case "${observed_cafa_commit}" in
    "${CAFA_ASSESSMENT_COMMIT}"*) ;;
    *)
      echo "CAFA assessment commit mismatch: ${observed_cafa_commit}" >&2
      exit 1
      ;;
  esac
  echo "==> Pinned CAFA_assessment_tool: ${observed_cafa_commit}"
fi
# --- 1b. Stage the CAFA3-era GO ontology expected at data/go.obo -------
# reproduce_embeddings_retrain_eval.sh checks `data/go.obo`, and train.py/eval
# parse it for GO DAG propagation. It ships in the CAFA tool cloned above as
# precrec/go_cafa3.obo (GO release 2016-05-31 — the correct ontology for CAFA3;
# do NOT substitute a current .obo or you contaminate the benchmark).
mkdir -p "${DATA_DIR}"
if [ ! -f "${DATA_DIR}/go.obo" ]; then
  cafa3_go_source="$(resolve_artifact_path cafa3_go_obo "${CAFA3_GO_OBO:-}" || true)"
  if [[ -n "$cafa3_go_source" ]]; then
    cp -p "$cafa3_go_source" "${DATA_DIR}/go.obo"
    echo "==> Staged data/go.obo from existing CAFA3 artifact: $cafa3_go_source"
  else
    cp "${CAFA_ASSESSMENT_DIR}/precrec/go_cafa3.obo" "${DATA_DIR}/go.obo"
    echo "==> Staged data/go.obo from CAFA_assessment_tool (go_cafa3.obo, 2016-05-31)"
  fi
fi

# --- 2. Raw CAFA3 CSVs — Zenodo 7409660 (verified bit-for-bit vs Zijan's splits).
#        Download 9 split CSVs, authenticate PRISTINE files against Zenodo md5s,
#        THEN normalise MF column 'protein' -> 'proteins' so prepare runs unmodified.
RAW="${PFP_CAFA3_RAW_DIR}"

if [ "${EMBEDDING_DEPENDENCY_PROFILE}" = "all" ]; then
  mkdir -p "${RAW}"
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
        role="cafa3_${aspect}_${split}"
        explicit_source=""
        if [[ -n "$CAFA3_SOURCE_DIR" ]]; then
          explicit_source="${CAFA3_SOURCE_DIR}/${aspect}-${split}.csv"
        fi
        artifact_source="$(resolve_artifact_path "$role" "$explicit_source" || true)"
        if [ -n "$artifact_source" ]; then
          echo "==> Staging ${aspect}-${split}.csv from existing artifact: $artifact_source"
          cp -p "$artifact_source" "$f"
        else
          echo "==> Downloading ${aspect}-${split}.csv"
          wget -c "${CAFA3_BASE}/${aspect}-${split}.csv?download=1" -O "$f"
        fi
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
else
  echo "==> Skipping CAFA3 CSV staging for ${EMBEDDING_DEPENDENCY_PROFILE} dependency profile"
fi

# --- 3. STRING files (PPI): alias (confirmed URL) + network embeddings .h5 (manual)
STRING_ALIAS="${STRING_ALIAS_FILE:-${PFP_STRING_DIR}/protein.aliases.v12.0.txt}"
STRING_ALIAS_GZ="${STRING_ALIAS}.gz"
if [ "${EMBEDDING_DEPENDENCY_PROFILE}" = "all" ] || \
   [ "${EMBEDDING_DEPENDENCY_PROFILE}" = "ppi" ]; then
  mkdir -p "${PFP_STRING_DIR}"
  if [ ! -f "${STRING_ALIAS}" ]; then
    alias_source="$(resolve_artifact_path string_aliases "${STRING_ALIAS_GZ_FILE:-${SAN_STRING_ALIAS_GZ:-}}" || true)"
    if [ -n "$alias_source" ]; then
      echo "==> Expanding STRING aliases v12.0 from existing artifact: $alias_source"
      gzip -dc "$alias_source" > "${STRING_ALIAS}.partial"
      mv "${STRING_ALIAS}.partial" "${STRING_ALIAS}"
    else
      echo "==> Downloading STRING aliases v12.0 (~3.2 GB)"
      wget -c "${STRING_DOWNLOAD_BASE}/protein.aliases.v12.0.txt.gz" -O "${STRING_ALIAS_GZ}"
      gunzip "${STRING_ALIAS_GZ}"
    fi
  else
    echo "==> STRING aliases already present, skipping"
  fi

  STRING_H5="$(resolve_artifact_path string_embeddings "${STRING_H5_FILE:-${SAN_STRING_H5:-}}" || true)"
  if [ -n "$STRING_H5" ]; then
    echo "==> Using existing STRING network embeddings: $STRING_H5"
  else
    STRING_H5="${PFP_STRING_DIR}/protein.network.embeddings.v12.0.h5"
  fi
  if [ ! -f "${STRING_H5}" ]; then
    echo "==> Downloading STRING network embeddings v12.0 (.h5) (~17.9 GB)"
    wget -c "${STRING_DOWNLOAD_BASE}/protein.network.embeddings.v12.0.h5" -O "${STRING_H5}"
  else
    echo "==> STRING network embeddings already present, skipping"
  fi
  add_mmfp_singularity_bind "$(dirname "$STRING_H5")"
else
  STRING_H5="${PFP_STRING_DIR}/protein.network.embeddings.v12.0.h5"
  echo "==> Skipping STRING staging for ${EMBEDDING_DEPENDENCY_PROFILE} dependency profile"
fi


# --- 4. AlphaFold + UniProt are runtime API downloads, nothing to pre-fetch.
echo "==> AlphaFold & UniProt are runtime API downloads (handled in modality scripts)."

# --- 5. Emit env-var exports.
mkdir -p "$(dirname "${DEPENDENCY_ENV}")"
cat > "${DEPENDENCY_ENV}" <<EOF
# Source before running modality scripts:  source ${DEPENDENCY_ENV}
export CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR}"
export CAFA3_RAW_DIR="${RAW}"
export STRING_H5_FILE="${STRING_H5}"
export STRING_ALIAS_FILE="${STRING_ALIAS}"
export SINGULARITY_BINDPATH="${SINGULARITY_BINDPATH:-}"
EOF

echo ""
echo "==> Wrote env exports to ${DEPENDENCY_ENV}"
