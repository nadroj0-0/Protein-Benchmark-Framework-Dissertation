#!/usr/bin/env bash
# Rebuild the shared mmfp entrypoint using a newer userspace inside Singularity.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FRAMEWORK_ROOT="$(cd "${HERE}/../.." && pwd)"
# shellcheck source=../reproduction_common.sh
source "${FRAMEWORK_ROOT}/scripts/reproduction_common.sh"
load_framework_paths "$FRAMEWORK_ROOT"

[[ "${REBUILD_MMFP:-}" == "YES" ]] || {
  echo "Refusing destructive rebuild. Set REBUILD_MMFP=YES after stopping jobs that use mmfp." >&2
  exit 2
}
[[ -x "$CONDA_EXE" ]] || { echo "Missing Conda executable: $CONDA_EXE" >&2; exit 1; }
SINGULARITY_BIN="${SINGULARITY_BIN:-/usr/bin/singularity}"
[[ -x "$SINGULARITY_BIN" ]] || { echo "Missing Singularity executable: $SINGULARITY_BIN" >&2; exit 1; }

eval "$("$CONDA_EXE" shell.bash hook)"
conda env remove -y -n "$MMFP_ENV" >/dev/null 2>&1 || true
rm -rf -- "$MMFP_ENV_DIR"
mkdir -p "$MMFP_SINGULARITY_DIR" "${MMFP_SINGULARITY_DIR}/cache"

conda create -y -n "$MMFP_ENV" \
  "python=${MMFP_PYTHON}" pip setuptools wheel squashfs-tools

if [[ ! -s "$MMFP_SINGULARITY_IMAGE" || "${MMFP_REFRESH_IMAGE:-0}" == "1" ]]; then
  PATH="${MMFP_ENV_DIR}/bin:${PATH}" \
    SINGULARITY_CACHEDIR="${MMFP_SINGULARITY_DIR}/cache" \
    "$SINGULARITY_BIN" pull --force \
      "$MMFP_SINGULARITY_IMAGE" "$MMFP_SINGULARITY_IMAGE_URI"
fi

rm -rf -- "$MMFP_SINGULARITY_VENV"
"$SINGULARITY_BIN" exec "$MMFP_SINGULARITY_IMAGE" \
  python -m venv "$MMFP_SINGULARITY_VENV"

bin_dir="${MMFP_ENV_DIR}/bin"
mv "${bin_dir}/python3.9" "${bin_dir}/python3.9.host"
mv "${bin_dir}/pip" "${bin_dir}/pip.host"
mv "${bin_dir}/pip3" "${bin_dir}/pip3.host"
install -m 755 "${HERE}/mmfp_singularity_python.sh" "${bin_dir}/python3.9"
install -m 755 "${HERE}/mmfp_singularity_pip.sh" "${bin_dir}/pip"
ln -sfn pip "${bin_dir}/pip3"
ln -sfn pip "${bin_dir}/pip3.9"

install_mmfp_packages "${bin_dir}/python"
conda activate "$MMFP_ENV"
validate_mmfp_env "$(command -v python)"

echo "==> MMFP Singularity image"
sha256sum "$MMFP_SINGULARITY_IMAGE"
du -sh "$MMFP_SINGULARITY_DIR" "$MMFP_ENV_DIR"
