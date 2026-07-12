#!/usr/bin/env bash
# Shared shell helpers for the root-level PFP reproduction wrappers.

load_framework_paths() {
  local framework_root="$1"
  if [ -f "${framework_root}/configs/paths.local.sh" ]; then
    # Machine-specific paths are intentionally not committed.
    # shellcheck disable=SC1091
    source "${framework_root}/configs/paths.local.sh"
  fi

  export PFP_GIT_URL="${PFP_GIT_URL:-https://github.com/psipred/PFP.git}"
  export PFP_CLONE_DIR="${PFP_CLONE_DIR:-PFP}"
  export CONDA_EXE="${CONDA_EXE:-/share/apps/miniforge3_mamba/bin/conda}"
  export MMFP_ENV="${MMFP_ENV:-mmfp}"
  export MMFP_ENV_DIR="${MMFP_ENV_DIR:-${HOME}/.conda/envs/${MMFP_ENV}}"
  export MMFP_PYTHON="${MMFP_PYTHON:-3.11}"
}

clone_or_reuse_pfp() {
  if [ -e "${PFP_CLONE_DIR}" ] && [ ! -d "${PFP_CLONE_DIR}/.git" ]; then
    echo "PFP_CLONE_DIR exists but is not a git checkout: ${PFP_CLONE_DIR}" >&2
    exit 1
  fi

  if [ -d "${PFP_CLONE_DIR}/.git" ]; then
    echo "==> Using existing PFP checkout: ${PFP_CLONE_DIR}"
  else
    echo "==> Cloning PFP from ${PFP_GIT_URL} into ${PFP_CLONE_DIR}"
    git clone "${PFP_GIT_URL}" "${PFP_CLONE_DIR}"
  fi

  cd "${PFP_CLONE_DIR}"
}

activate_or_create_mmfp_env() {
  if [ ! -x "${CONDA_EXE}" ]; then
    echo "Missing conda executable: ${CONDA_EXE}" >&2
    echo "Set CONDA_EXE in configs/paths.local.sh or the environment." >&2
    exit 1
  fi

  eval "$("${CONDA_EXE}" shell.bash hook)"

  if [ ! -d "${MMFP_ENV_DIR}" ]; then
    echo "==> Creating Conda environment: ${MMFP_ENV}"

    conda create -y -n "${MMFP_ENV}" "python=${MMFP_PYTHON}"
    conda activate "${MMFP_ENV}"

    python -m pip install --upgrade pip setuptools wheel
    pip install -r requirements.txt --prefer-binary
    pip install requests fair-esm biopython protobuf sentencepiece torch-geometric "biotite==0.41.2"
    pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-2.6.0+cu124.html
    pip install --only-binary=:all: h5py
    pip install --only-binary=:all: tiktoken

    echo "==> Environment created."
  else
    echo "==> Using existing Conda environment: ${MMFP_ENV}"
  fi

  conda activate "${MMFP_ENV}"
}
