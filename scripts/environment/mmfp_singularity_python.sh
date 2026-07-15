#!/usr/bin/env bash
set -euo pipefail

image="${MMFP_SINGULARITY_IMAGE:-$HOME/.mmfp_singularity/python-3.9.23.sif}"
venv="${MMFP_SINGULARITY_VENV:-$HOME/.mmfp_singularity/venv}"
singularity_bin="${SINGULARITY_BIN:-/usr/bin/singularity}"

[[ -x "$singularity_bin" ]] || {
  echo "Missing Singularity executable: $singularity_bin" >&2
  exit 1
}
[[ -s "$image" ]] || {
  echo "Missing MMFP Singularity image: $image" >&2
  exit 1
}
[[ -f "$venv/pyvenv.cfg" ]] || {
  echo "Missing MMFP container virtual environment: $venv" >&2
  exit 1
}

exec "$singularity_bin" exec --nv "$image" "$venv/bin/python" "$@"
