#!/usr/bin/env bash
# Shared shell helpers for the root-level PFP reproduction wrappers.

# shellcheck disable=SC1091
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/artifact_catalog.sh"

readonly MMFP_PFP_COMMIT="1e04fd6d6d3c40458fd41ec1a881ed6e24de768e"
readonly MMFP_CAFA_ASSESSMENT_COMMIT="d72f0a5abb66d3224bd808e2015b55f1c9d18340"

load_framework_paths() {
  local framework_root="$1"
  if [ -f "${framework_root}/configs/paths.local.sh" ]; then
    # Machine-specific paths are intentionally not committed.
    # shellcheck disable=SC1091
    source "${framework_root}/configs/paths.local.sh"
  fi

  export PFP_GIT_URL="${PFP_GIT_URL:-https://github.com/psipred/PFP.git}"
  export PFP_CLONE_DIR="${PFP_CLONE_DIR:-PFP}"

  artifact_catalog_configure "$framework_root" "${ARTIFACT_CATALOG:-}"
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
    "fair-esm": "2.0.0",
    "protobuf": "6.33.6",
    "torch-geometric": "2.6.1",
    "torch-scatter": "2.1.2",
    "torch-sparse": "0.6.18",
    "tiktoken": "0.13.0",
}
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
    comparable = (
        value.split("+", 1)[0]
        if distribution in {"torch", "torch-scatter", "torch-sparse"}
        else value
    )
    if comparable != expected:
        errors.append(f"{distribution}: expected {expected}, found {value}")

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
for distribution in EXPECTED:
    print(f"    {distribution}=={observed[distribution]}")
PY

  "$python_bin" -m pip check
}

validate_mmfp_if1_env() {
  local python_bin="$1"
  local numpy_overlay="$2"

  [[ -n "$python_bin" && -x "$python_bin" ]] || {
    echo "Cannot validate IF1 runtime: Python executable is unavailable: $python_bin" >&2
    return 1
  }
  [[ -d "$numpy_overlay/numpy" ]] || {
    echo "Cannot validate IF1 runtime: NumPy overlay is unavailable: $numpy_overlay" >&2
    return 1
  }

  PYTHONPATH="$numpy_overlay${PYTHONPATH:+:$PYTHONPATH}" \
    "$python_bin" - <<'PY'
from importlib.metadata import PackageNotFoundError, version
import json
import platform

EXPECTED = {
    "numpy": "1.26.4",
    "biotite": "0.38.0",
}

observed = {}
errors = []
for distribution, expected in EXPECTED.items():
    try:
        observed[distribution] = version(distribution)
    except PackageNotFoundError:
        errors.append(f"{distribution}: missing (expected {expected})")
        continue
    if observed[distribution] != expected:
        errors.append(
            f"{distribution}: expected {expected}, found {observed[distribution]}"
        )

try:
    observed["fair-esm"] = version("fair-esm")
except PackageNotFoundError:
    errors.append("fair-esm: missing")

for module in ("biotite.structure", "esm.inverse_folding"):
    try:
        __import__(module)
    except Exception as exc:
        errors.append(f"import {module}: {type(exc).__name__}: {exc}")

if errors:
    for error in errors:
        print(f"IF1 runtime error: {error}")
    raise SystemExit(1)

print(json.dumps({
    "biotite": observed["biotite"],
    "fair_esm": observed["fair-esm"],
    "numpy": observed["numpy"],
    "numpy_overlay": True,
    "python": platform.python_version(),
}, indent=2, sort_keys=True))
PY
}

install_mmfp_if1_numpy_overlay() {
  local python_bin="$1"
  local numpy_overlay="$2"

  rm -rf -- "$numpy_overlay"
  mkdir -p "$numpy_overlay"
  "$python_bin" -m pip install \
    --disable-pip-version-check \
    --no-deps \
    --only-binary=:all: \
    --target "$numpy_overlay" \
    "numpy==1.26.4"
}
