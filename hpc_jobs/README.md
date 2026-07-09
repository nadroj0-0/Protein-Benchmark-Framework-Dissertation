# HPC Jobs

This directory contains cluster submission wrappers for running the
framework on UCL/SGE machines.

Use `scripts/` for reusable implementation logic. Use `hpc_jobs/` for
`qsub` entrypoints that request resources, prepare scratch space, clone
the framework, run a workflow, and copy results home.

## Layout

```text
hpc_jobs/
├── active/    # Current qsub wrappers used for reproduction jobs
├── examples/  # Scheduler examples/templates
└── archive/   # Historical scripts kept for provenance
```

## Active Jobs

Submit active wrappers from the repository root or by giving `qsub` the
full path:

```bash
qsub hpc_jobs/active/hpc_reproduce_eval_only.sh
qsub hpc_jobs/active/hpc_reproduce_retrain_eval.sh
qsub hpc_jobs/active/hpc_reproduce_embeddings_retrain_eval.sh
qsub hpc_jobs/active/hpc_cafa3_deepgoplus_pickle_generation_validation.sh
qsub hpc_jobs/active/hpc_cafa3_deepgoplus_validation.sh
qsub hpc_jobs/active/hpc_cafa3_historical_validation.sh
```

The active wrappers clone the full framework into node-local scratch and
then call the normal entrypoints under `scripts/`.

For benchmark validation:

- `hpc_cafa3_deepgoplus_pickle_generation_validation.sh` validates the
  official CAFA3/DeepGOPlus file-to-pickle layer.
- `hpc_cafa3_deepgoplus_validation.sh` validates the released DeepGOPlus/TEMPROT
  intermediate path and is the preferred lightweight historical validation.
- `hpc_cafa3_historical_validation.sh` runs the heavier raw-snapshot audit.

Both validation wrappers copy only reports/logs back to home and remove scratch
data at the end of the job.
