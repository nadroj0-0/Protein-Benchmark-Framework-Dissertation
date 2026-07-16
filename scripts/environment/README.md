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

`mmfp_singularity_python.sh` also accepts `MMFP_PYTHONPATH` for an explicitly
scoped container-side overlay. The contemporary embedding workflow uses this
only for its scratch-local NumPy 1.26.4 IF1 compatibility layer; the default
launcher and the primary MMFP environment remain on NumPy 2.0.2.

The base image intentionally remains minimal and does not install Git. Homology
HPC wrappers verify their scratch checkout with host Git before invoking Python,
then pass the verified commit, clean state, and exact repository path into the
container. The builder validates and records that state without requiring a
second Git installation. Direct/local homology runs continue to inspect Git
normally.
