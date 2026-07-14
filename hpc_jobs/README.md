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
├── launchers/ # Reviewed dry-run/pilot/array launchers for guarded workflows
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
qsub -v BENCHMARK_DIR=/path/to/contemporary/run/outputs \
  hpc_jobs/active/hpc_contemporary_embedding_inventory.sh
```

The active wrappers clone the full framework into node-local scratch and
then call the normal entrypoints under `scripts/`.

The historical and contemporary benchmark-generation wrappers activate and use
the shared `mmfp` environment directly. They do not create another virtual
environment or replace its NumPy and pandas installations.

`hpc_contemporary_embedding_inventory.sh` is a CPU-only integration test for
the completed nine contemporary CSVs against Zijian's published embeddings.
The benchmark path is supplied at submission time, never embedded in the
wrapper. The job downloads Zenodo records 7409660 and 19498341 into scratch,
authenticates the nine source CSVs and three embedding archives, validates the
published cache, and runs the existing provenance-aware planner. Only compact
reports, manifests, ID lists, and `generate_prott5.fasta` are copied to
`$HOME/contemporary_embedding_inventory_results`; archive/cache bytes are not.
Scratch is removed after success, failure, or termination.

The common submission form is:

```bash
qsub -v BENCHMARK_DIR="$HOME/contemporary_cafa_benchmark_results/<run>/outputs" \
  hpc_jobs/active/hpc_contemporary_embedding_inventory.sh
```

Any role can be overridden with `BP_TRAINING_CSV`, `BP_VALIDATION_CSV`,
`BP_TEST_CSV`, and corresponding `CC_*`/`MF_*` variables. With no
`BENCHMARK_DIR`, all nine variables are required. Optional overrides are
`SOURCE_BENCHMARK_DIR`, `PUBLISHED_EMBEDDING_ARCHIVE_DIR`, `ALIASES_FILE`,
`INVENTORY_CONFIG`, `INVENTORY_POLICY`, `REPORT_LEVEL`, and `RESULTS_ROOT`.
The default is to download the canonical source CSVs and published embedding
archives during each job, as requested.

## Homology-cluster identity array

The homology benchmark has three components:

```text
hpc_jobs/active/hpc_homology_cluster_benchmark.sh
hpc_jobs/launchers/submit_homology_cluster_pilot.sh
hpc_jobs/launchers/submit_homology_cluster_array.sh
```

The active file is the array-aware worker; scientific work remains in
`scripts/benchmark_generation/run_homology_cluster_benchmark.sh` and the isolated Python package.
Both launchers source `_homology_cluster_common.sh` for strict shared-input, hash, scope, revision,
methodology, and preview validation.

The deterministic mapping is:

```text
1 -> 30%    2 -> 25%    3 -> 20%    4 -> 15%    5 -> 10%    6 -> 5%
```

The pilot launcher constructs `qsub -t 1 -pe smp 8`; the full launcher constructs
`qsub -t 1-6 -pe smp 8`. `NSLOTS` is authoritative within each task and must equal the MMseqs2
thread count. Six runnable tasks can request up to 48 CPU slots, but scheduler concurrency remains
a Grid Engine decision.

Production arrays require:

- an exact full lowercase framework commit and detached clean checkout;
- one explicit `sprot-only`, `trembl-only`, or `sprot-and-trembl` scope;
- a reviewed scope-aware 5/5/6 frozen manifest;
- shared local checksum-verified inputs and `NO_DOWNLOADS=1`;
- fixed split, annotated-only population, seed, evidence/ontology policy, and exact MMseqs2 version;
- a successful validated task-1/30% diagnostic pilot;
- reviewed attrition policy, task context, separately sourced runtime/memory/scratch/output
  measurements, and human approval.

The full launcher rechecks approval evidence before constructing the command. The queued worker
rechecks launcher-time evidence hashes and reruns approval from its detached checkout before input
staging. Approval binds the reviewed attrition-policy hash and a common pilot run ID across the
marker, task context, and measurement evidence. Authorization reconstructs the pilot ratios and
requires the reviewed policy to accept them; diagnostic reports remain explicitly non-authorizing.
No bypass exists. One array-wide attrition override is rejected because observed failures can
differ across identities.

Task paths include source scope, framework revision, run ID, Grid Engine job ID, task ID, identity,
and split policy. Scratch and persistent paths are atomically claimed. Success, failures, INT, TERM,
and copy failure attempt marker-free diagnostics and then remove only task-owned scratch. Unsafe,
pre-existing, symlinked, or ownership-mismatched cleanup targets are refused; a collision cannot
mutate an earlier publication.

Preview commands are safe and print the exact command plus every exported value:

```bash
DRY_RUN=1 bash hpc_jobs/launchers/submit_homology_cluster_pilot.sh
DRY_RUN=1 bash hpc_jobs/launchers/submit_homology_cluster_array.sh
```

The full preview still validates the approval; only `qsub` is suppressed. Automated tests place a
recording fake `qsub` first in `PATH`, test controlled scheduler failure, and prove dry run does not
call it. Neither launcher was executed against Grid Engine during implementation. No real `qsub`,
`qstat`, `qdel`, SSH, full download, or full MMseqs2 run occurred.

The exact future nine-step user workflow—shared input declaration, pilot preview/submission,
validation/review, approval creation, full preview/submission, monitoring, and six-run
aggregation—is documented in
[`benchmark_builders/homology_cluster/README.md`](../benchmark_builders/homology_cluster/README.md).

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
