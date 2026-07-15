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

singularity_args=(exec --nv)

# Cluster jobs run from node-local scratch. Bind the resolved filesystem at
# both its public alias and canonical path so Python can open either form.
scratch_alias="/scratch0"
scratch_target="$(readlink -f "$scratch_alias" 2>/dev/null || true)"
if [[ -n "$scratch_target" && -d "$scratch_target" ]]; then
  singularity_args+=(--bind "$scratch_target:$scratch_alias")
  if [[ "$scratch_target" != "$scratch_alias" ]]; then
    singularity_args+=(--bind "$scratch_target:$scratch_target")
  fi
fi

container_command=("$venv/bin/python")
if [[ -n "${MMFP_PYTHONPATH:-}" ]]; then
  container_command=(env "PYTHONPATH=${MMFP_PYTHONPATH}" "$venv/bin/python")
fi

exec "$singularity_bin" "${singularity_args[@]}" "$image" \
  "${container_command[@]}" "$@"
