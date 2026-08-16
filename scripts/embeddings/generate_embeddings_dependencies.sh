#!/bin/bash
# generate_embeddings_dependencies.sh
# Downloads external resources for from-scratch embedding generation and writes
# PFP/external/dependency_env.sh with the env-var exports. The separately pinned
# CAFA Assessment Tool checkout defaults outside PFP so it cannot contaminate
# PFP source authentication.
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
PFP_ROOT="$(cd "$PFP_ROOT" && pwd -P)"
verify_clean_pinned_git_checkout "$PFP_ROOT" "$MMFP_PFP_COMMIT" "PFP" || exit 1
EXT="${PFP_EXTERNAL_DIR:-${PFP_ROOT}/external}"
DATA_DIR="${PFP_DATA_DIR:-${PFP_ROOT}/data}"
DEPENDENCY_ENV="${DEPENDENCY_ENV:-${EXT}/dependency_env.sh}"
CAFA_ASSESSMENT_REPO_URL="${CAFA_ASSESSMENT_REPO_URL:-https://github.com/ashleyzhou972/CAFA_assessment_tool.git}"
CAFA_ASSESSMENT_DIR="${CAFA_ASSESSMENT_DIR:-$(dirname "$PFP_ROOT")/CAFA_assessment_tool}"
CAFA_ASSESSMENT_COMMIT="$MMFP_CAFA_ASSESSMENT_COMMIT"
PFP_CAFA3_RAW_DIR="${PFP_CAFA3_RAW_DIR:-${EXT}/cafa3_raw}"
PFP_STRING_DIR="${PFP_STRING_DIR:-${EXT}/string}"
CAFA3_BASE="${CAFA3_BASE:-https://zenodo.org/records/7409660/files}"
STRING_DOWNLOAD_BASE="${STRING_DOWNLOAD_BASE:-https://stringdb-downloads.org/download}"
CAFA3_SOURCE_DIR="${CAFA3_SOURCE_DIR:-}"
EMBEDDING_DEPENDENCY_PROFILE="${EMBEDDING_DEPENDENCY_PROFILE:-all}"
CAFA_ASSESSMENT_NAME="$(basename "$CAFA_ASSESSMENT_DIR")"
[[ -n "$CAFA_ASSESSMENT_NAME" && "$CAFA_ASSESSMENT_NAME" != "." && \
   "$CAFA_ASSESSMENT_NAME" != ".." ]] || {
  echo "Invalid CAFA Assessment Tool destination: ${CAFA_ASSESSMENT_DIR}" >&2
  exit 1
}
PYTHON_BIN="${PYTHON_BIN:-python3}"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || {
  echo "Python is required to resolve the CAFA Assessment Tool destination" >&2
  exit 1
}
CAFA_ASSESSMENT_DIR="$(
  "$PYTHON_BIN" -c \
    'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve(strict=False))' \
    "$CAFA_ASSESSMENT_DIR"
)"
case "$CAFA_ASSESSMENT_DIR/" in
  "$PFP_ROOT/"|"$PFP_ROOT/"*)
    echo "CAFA Assessment Tool must be outside the authenticated PFP checkout" >&2
    exit 1
    ;;
esac
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
  CAFA_ASSESSMENT_DIR="$(cd "${CAFA_ASSESSMENT_DIR}" && pwd -P)"
  cafa_top_level="$(git -C "${CAFA_ASSESSMENT_DIR}" rev-parse --show-toplevel 2>/dev/null)" || {
    echo "CAFA Assessment Tool must be a Git checkout: ${CAFA_ASSESSMENT_DIR}" >&2
    exit 1
  }
  [[ "$(cd "$cafa_top_level" && pwd -P)" == "$CAFA_ASSESSMENT_DIR" ]] || {
    echo "CAFA Assessment Tool root is not the checkout root: ${CAFA_ASSESSMENT_DIR}" >&2
    exit 1
  }
  git_in_dir "${CAFA_ASSESSMENT_DIR}" checkout --detach "${CAFA_ASSESSMENT_COMMIT}"
  verify_clean_pinned_git_checkout \
    "${CAFA_ASSESSMENT_DIR}" "${CAFA_ASSESSMENT_COMMIT}" "CAFA Assessment Tool" || exit 1
  observed_cafa_commit="$(git_in_dir "${CAFA_ASSESSMENT_DIR}" rev-parse HEAD)"
  echo "==> Pinned CAFA_assessment_tool: ${observed_cafa_commit}"
fi
# --- 1b. Stage the CAFA3-era GO ontology expected at data/go.obo -------
# CAFA3 dependency setup stages `data/go.obo` where upstream train.py/eval
# expect it for GO DAG propagation. It ships in the CAFA tool cloned above as
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
    alias_source="$(resolve_artifact_path string_aliases "${STRING_ALIAS_GZ_FILE:-}" || true)"
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

  STRING_H5="$(resolve_artifact_path string_embeddings "${STRING_H5_FILE:-}" || true)"
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
EOF

echo ""
echo "==> Wrote env exports to ${DEPENDENCY_ENV}"
