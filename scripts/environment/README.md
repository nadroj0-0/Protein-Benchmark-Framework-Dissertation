# MMFP environment

`rebuild_mmfp_singularity.sh` recreates the shared `mmfp` Conda entrypoint while
running its scientific Python environment inside Singularity. This is required
on CentOS 7/glibc 2.17 hosts because the official PyTorch 2.8 wheels require
glibc 2.28 or newer.

Stop every job that may use `mmfp`, then run:

```bash
REBUILD_MMFP=YES bash scripts/environment/rebuild_mmfp_singularity.sh
```

The Conda environment name and `$HOME/.conda/envs/mmfp/bin/python` contract stay
unchanged. The Python and pip entrypoints are thin launchers for a Python 3.9.23
venv under `$HOME/.mmfp_singularity`. Package versions remain centralized in
`scripts/reproduction_common.sh`; this directory does not define a second set of
constraints.

The launcher requests `--nv` for GPU passthrough. A login node without a GPU may
print `Could not find any nv files`; this is expected. CUDA availability must be
confirmed in a scheduled GPU smoke test.
