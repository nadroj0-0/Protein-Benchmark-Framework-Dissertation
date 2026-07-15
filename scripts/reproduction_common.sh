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
  export MMFP_PYTHON="${MMFP_PYTHON:-3.9.23}"
  export MMFP_TORCH_INDEX_URL="${MMFP_TORCH_INDEX_URL:-https://download.pytorch.org/whl/cu126}"
  export MMFP_PYG_WHEEL_BASE="${MMFP_PYG_WHEEL_BASE:-https://data.pyg.org/whl}"
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

validate_mmfp_env() {
  local python_bin="${1:-$(command -v python)}"
  [[ -n "$python_bin" && -x "$python_bin" ]] || {
    echo "Cannot validate mmfp: Python executable is unavailable: $python_bin" >&2
    return 1
  }

  "$python_bin" - <<'PY'
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
import platform
import sys

EXPECTED_PYTHON = "3.9.23"
EXPECTED = {
    "torch": "2.8.0",
    "numpy": "2.0.2",
    "pandas": "2.3.3",
    "scipy": "1.13.1",
    "tqdm": "4.67.1",
    "scikit-learn": "1.6.1",
    "cafaeval": "1.2.1",
    "obonet": "1.1.1",
    "networkx": "3.2.1",
    "transformers": "4.56.2",
    "sentencepiece": "0.2.1",
    "biopython": "1.85",
    "h5py": "3.14.0",
    "requests": "2.32.5",
    "biotite": "0.38.0",
}
REQUIRED_UNPINNED = (
    "fair-esm",
    "protobuf",
    "torch-geometric",
    "torch-scatter",
    "torch-sparse",
    "tiktoken",
)
IMPORTS = (
    "torch",
    "numpy",
    "pandas",
    "scipy",
    "tqdm",
    "sklearn",
    "cafaeval",
    "obonet",
    "networkx",
    "transformers",
    "sentencepiece",
    "Bio",
    "h5py",
    "requests",
    "biotite",
    "esm",
    "torch_geometric",
    "torch_scatter",
    "torch_sparse",
    "tiktoken",
    "google.protobuf",
)

errors = []
observed_python = platform.python_version()
if observed_python != EXPECTED_PYTHON:
    errors.append(f"Python: expected {EXPECTED_PYTHON}, found {observed_python}")

observed = {}
for distribution, expected in EXPECTED.items():
    try:
        value = version(distribution)
    except PackageNotFoundError:
        errors.append(f"{distribution}: missing (expected {expected})")
        continue
    observed[distribution] = value
    comparable = value.split("+", 1)[0] if distribution == "torch" else value
    if comparable != expected:
        errors.append(f"{distribution}: expected {expected}, found {value}")

for distribution in REQUIRED_UNPINNED:
    try:
        observed[distribution] = version(distribution)
    except PackageNotFoundError:
        errors.append(f"{distribution}: missing")

for module in IMPORTS:
    try:
        import_module(module)
    except Exception as exc:
        errors.append(f"import {module}: {type(exc).__name__}: {exc}")

if "torch" in sys.modules:
    import torch

    torch_base = observed.get("torch", "").split("+", 1)[0]
    parts = torch_base.split(".")
    pytorch_tag = "pt" + "".join(parts[:2]) if len(parts) >= 2 else ""
    cuda_tag = "cpu" if torch.version.cuda is None else "cu" + torch.version.cuda.replace(".", "")
    for distribution in ("torch-scatter", "torch-sparse"):
        extension_version = observed.get(distribution, "")
        wanted_tag = f"+{pytorch_tag}{cuda_tag}"
        if extension_version and wanted_tag not in extension_version:
            errors.append(
                f"{distribution}: build {extension_version} does not match "
                f"PyTorch {torch_base}/{cuda_tag}; expected tag {wanted_tag}"
            )

if errors:
    print("mmfp environment validation FAILED:", file=sys.stderr)
    for error in errors:
        print(f"  - {error}", file=sys.stderr)
    raise SystemExit(1)

print(f"==> mmfp environment validated: Python {observed_python}")
for distribution in (*EXPECTED, *REQUIRED_UNPINNED):
    print(f"    {distribution}=={observed[distribution]}")
PY

  "$python_bin" -m pip check
}

validate_mmfp_creation_host() {
  if [[ "$(uname -s)" != "Linux" || "$MMFP_TORCH_INDEX_URL" != *download.pytorch.org* ]]; then
    return 0
  fi

  local glibc_version
  glibc_version="$(getconf GNU_LIBC_VERSION 2>/dev/null | awk '{print $2}')"
  if [[ -n "$glibc_version" && "$(printf '%s\n' "$glibc_version" 2.28 | sort -V | head -n 1)" != "2.28" ]]; then
    cat >&2 <<EOF
Cannot create the exact mmfp environment natively on this host.
PyTorch 2.8 official Linux wheels require glibc 2.28 or newer; this host has glibc ${glibc_version}.
Use a compatible container runtime rather than weakening Zijian's Python/PyTorch pins.
EOF
    return 1
  fi
}

activate_or_create_mmfp_env() {
  local pyg_wheel_url

  if [ ! -x "${CONDA_EXE}" ]; then
    echo "Missing conda executable: ${CONDA_EXE}" >&2
    echo "Set CONDA_EXE in configs/paths.local.sh or the environment." >&2
    exit 1
  fi

  eval "$("${CONDA_EXE}" shell.bash hook)"

  if [ ! -d "${MMFP_ENV_DIR}" ]; then
    validate_mmfp_creation_host
    echo "==> Creating Conda environment: ${MMFP_ENV}"

    conda create -y -n "${MMFP_ENV}" "python=${MMFP_PYTHON}"
    conda activate "${MMFP_ENV}"

    python -m pip install --upgrade pip setuptools wheel
    python -m pip install \
      --index-url "${MMFP_TORCH_INDEX_URL}" \
      "torch==2.8.0"
    python -m pip install --prefer-binary \
      "numpy==2.0.2" \
      "pandas==2.3.3" \
      "scipy==1.13.1" \
      "tqdm==4.67.1" \
      "scikit-learn==1.6.1" \
      "cafaeval==1.2.1" \
      "obonet==1.1.1" \
      "networkx==3.2.1" \
      "transformers==4.56.2" \
      "sentencepiece==0.2.1" \
      "biopython==1.85" \
      "requests==2.32.5" \
      "biotite==0.38.0" \
      fair-esm protobuf torch-geometric tiktoken
    python -m pip install --only-binary=:all: "h5py==3.14.0"

    pyg_wheel_url="$(python - "${MMFP_PYG_WHEEL_BASE}" <<'PY'
import sys
import torch

base = sys.argv[1].rstrip("/")
cuda = "cpu" if torch.version.cuda is None else "cu" + torch.version.cuda.replace(".", "")
print(f"{base}/torch-2.8.0+{cuda}.html")
PY
)"
    echo "==> Installing PyTorch-Geometric extensions from: ${pyg_wheel_url}"
    python -m pip install --only-binary=:all: \
      torch-scatter torch-sparse -f "${pyg_wheel_url}"

    echo "==> Environment created."
  else
    echo "==> Using existing Conda environment: ${MMFP_ENV}"
  fi

  conda activate "${MMFP_ENV}"
  validate_mmfp_env "$(command -v python)"
}
