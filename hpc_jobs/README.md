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
qsub hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh
```

The active wrappers clone the full framework into node-local scratch and
then call the normal entrypoints under `scripts/`.

The historical and contemporary benchmark-generation wrappers activate and use
the shared `mmfp` environment directly. They do not create another virtual
environment or replace its NumPy and pandas installations.

`hpc_contemporary_temporal_benchmark.sh` stages any locally available frozen
2025/2026 UniProt, GOA and GO inputs, downloads missing inputs into scratch,
stream-filters full TrEMBL sources to the CAFA3 taxa, and invokes
`scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh`, copies
the complete benchmark run to durable storage and clears scratch. It builds the
benchmark only; it does not launch PFP.

The production contemporary defaults compare the two frozen GOA files by direct
snapshot membership and do not remove t1-only rows using their assigned date:

```text
T1_ENDPOINT_POLICY=snapshot-membership
T1_BACKFILL_POLICY=allow
```

Previously filtered TrEMBL products can be supplied with
`T0_TREMBL_FILTERED_INPUT` and `T1_TREMBL_FILTERED_INPUT`. The wrapper copies
them into scratch under the expected release layout, avoiding another download
and stream-filter of the complete TrEMBL releases.

The wrapper is a `.sh` file because that is the repository convention. Grid
Engine does not require `.qsub`; the `qsub` command uses the embedded `#$`
directives regardless of the filename suffix.

`PROTEIN_DATABASE_ROOT` can point at a persistent database tree. Local files are
an optional cache, not a requirement. If no files are available, the job uses
the official frozen URLs and keeps only the filtered/required products in
scratch. Colon-separated `UNIPROT_T0_INPUTS` and `UNIPROT_T1_INPUTS` are copied
to scratch when explicitly supplied.

For benchmark validation:

- `hpc_cafa3_deepgoplus_pickle_generation_validation.sh` validates the
  official CAFA3/DeepGOPlus file-to-pickle layer.
- `hpc_cafa3_deepgoplus_validation.sh` validates the released DeepGOPlus/TEMPROT
  intermediate path and is the preferred lightweight historical validation.
- `hpc_cafa3_historical_validation.sh` defaults to the released-groundtruth
  historical artifact audit, the closest public pre-freeze training source
  (UniProtKB 2016_08), and the released CAFA3 target FASTA, aggregate ground
  truth, and DeepGOPlus ontology. It retains its regenerated nine CSVs and five
  pickle intermediates with comparison reports. `raw-goa` enables the heavier
  archived-GOA and TrEMBL mapping audit when that distinct forensic question is
  required.

The primary historical submission is:

```bash
qsub -v HISTORICAL_TRAINING_SNAPSHOT=september-2016,TARGET_UNIVERSE_POLICY=official-cafa3-targets,HISTORICAL_TEST_SOURCE=official-groundtruth \
  hpc_jobs/active/hpc_cafa3_historical_validation.sh
```

Set `HISTORICAL_TEST_SOURCE=raw-goa` to rerun the heavier public-snapshot
forensic reconstruction instead. The released-groundtruth default bypasses the
large GOA/TrEMBL downloads because the curated CAFA test labels and target FASTA
are authoritative for that validation claim.

The bounded raw-GOA policy matrix is submitted with:

```bash
bash hpc_jobs/submit_cafa3_raw_experiment_matrix.sh
```

It submits four independent jobs. Together with the completed date-proxy /
backfill / February-ontology baseline, they isolate t1 snapshot membership,
pre-t0 backfill removal, and use of the packaged DeepGOPlus ontology. The
launcher deliberately does not resubmit the baseline. SGE schedules each job
normally; the script does not force simultaneous execution.

Missing historical sources are downloaded into scratch. Optional local archive
overrides are `HISTORICAL_TRAINING_UNIPROT_ARCHIVE` and
`OFFICIAL_CAFA3_ARCHIVE_INPUT`.

Validation wrappers remove scratch after completion. The lightweight validation
jobs copy reports/logs; the raw-snapshot audit additionally retains its generated
CSV and pickle artefacts.
