# HPC Jobs

## SAN frozen-input acquisition

`active/hpc_populate_san_frozen_inputs.sh` is the Grid Engine wrapper for the
idempotent persistent-input loader. It clones the committed framework into a
small job-owned scratch directory, writes large inputs directly to
`/SAN/bioinf/bmpfp`, and always removes the scratch checkout. Interrupted SAN
downloads remain as resumable `.partial` files.

Submit the full catalogue with:

```bash
qsub hpc_jobs/active/hpc_populate_san_frozen_inputs.sh
```

Select profiles without editing the wrapper:

```bash
qsub -v SAN_INPUT_PROFILES=homology,tools \
  hpc_jobs/active/hpc_populate_san_frozen_inputs.sh
```

The scheduler-neutral implementation and complete inventory are documented in
[`scripts/data_acquisition/README.md`](../scripts/data_acquisition/README.md).

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
qsub hpc_jobs/active/hpc_cafa3_full_from_scratch_reproduction.sh
qsub hpc_jobs/active/hpc_cafa3_deepgoplus_pickle_generation_validation.sh
qsub hpc_jobs/active/hpc_cafa3_deepgoplus_validation.sh
qsub hpc_jobs/active/hpc_cafa3_historical_validation.sh
qsub hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh
qsub -v BENCHMARK_DIR=/path/to/contemporary/run/outputs \
  hpc_jobs/active/hpc_contemporary_embedding_inventory.sh
qsub -v TARGET_BENCHMARK_DIR=/path/to/contemporary/run/outputs \
  hpc_jobs/active/hpc_contemporary_benchmark_reuse_plan.sh
```

### Full CAFA3 From-Scratch Reproduction

`hpc_cafa3_full_from_scratch_reproduction.sh` is the hardened replacement for
the older embedding/retrain/evaluate wrapper when the purpose is to audit the
whole chain. Submit it from a clean committed framework checkout:

```bash
qsub hpc_jobs/active/hpc_cafa3_full_from_scratch_reproduction.sh
```

It clones pinned framework and PFP commits into `/scratch0`, validates the
author-supplied MMFP package pins, downloads and authenticates the nine Zenodo
7409660 CSVs, builds `proteins.fasta`, and runs all four embedding modalities in
parallel after a bounded preflight. The CAFA assessment tool is pinned to the
official supplementary-code tag `v1.0-beta` at commit
`d72f0a5abb66d3224bd808e2015b55f1c9d18340`. The IF1 NumPy/device
compatibility path and PPI zero-CAFA-ID guard are isolated scratch copies;
upstream PFP is unchanged.

Only after regeneration finishes does the job download the three authenticated
embedding archives from Zenodo 19498341. It records byte hashes and numeric
comparisons for every union protein/modality, compresses the row-level table,
and deletes both the extracted published cache and its archives. Fresh training
therefore uses regenerated arrays only. It requires all three checkpoints and
evaluation summaries, while metric or embedding differences are retained as
results rather than disguised as infrastructure failures.

Successful compact results are published under
`$HOME/cafa3_full_reproduction_results/<job>_<timestamp>/`:

```text
cafa3_full_reproduction_report.md
cafa3_full_reproduction_report.json
WORKFLOW_COMPLETE.json
reports/embedding_comparison.csv.gz
reports/embedding_comparison_summary.json
reports/evaluation/
reports/training/
reports/input_acquisition.tsv
reports/modality_status.tsv
logs/
```

Transient generated arrays, published embedding archives/cache, model caches,
STRING inputs and checkpoints stay in scratch and are removed unconditionally
after compact publication, including when publication fails. Arrays that pass
the independent validator are the sole exception: one copy is atomically
published into the provenance-bound SAN embedding state before scratch cleanup.
Optional CLI overrides are
`--results-root`, `--embedding-state-root`, `--embedding-mode`, and
`--text-cutoff-date`; the historical default cutoff is `2016-02-17`.

The initial run now preserves validated arrays in one provenance-bound SAN
state. If any modality remains below the published historical coverage floor,
the result is published with an `.incomplete` suffix and
`GENERATION_INCOMPLETE.json`; training and evaluation do not start. This is a
completed acquisition attempt, not an infrastructure failure.

Retry exactly one missing modality at a time:

```bash
qsub hpc_jobs/active/hpc_cafa3_embedding_retry.sh --modality structure
qsub hpc_jobs/active/hpc_cafa3_embedding_retry.sh --modality text
qsub hpc_jobs/active/hpc_cafa3_embedding_retry.sh --modality ppi
qsub hpc_jobs/active/hpc_cafa3_embedding_retry.sh --modality sequence
```

Only modalities listed in the current SAN `needs_retry.tsv` need a job. The
retry wrapper requests one GPU, recreates canonical inputs in scratch, runs a
20-protein subset-equivalence control, merges only validated successes, copies
compact reports to `$HOME/cafa3_embedding_retry_results`, and always removes
scratch. It does not resubmit itself and does not trigger training.

After the SAN marker passes, continue the original audit without regenerating
accepted embeddings:

```bash
qsub hpc_jobs/active/hpc_cafa3_full_from_scratch_reproduction.sh \
  --embedding-mode resume
```

The default persistent state is
`/SAN/bioinf/bmpfp/embedding_states/cafa3_full_reproduction`; use
`--embedding-state-root` on both wrappers to select another explicit location.

The active wrappers clone the full framework into node-local scratch and
then call the normal entrypoints under `scripts/`.

The historical and contemporary benchmark-generation wrappers activate and use
the shared `mmfp` environment directly. They do not create another virtual
environment or replace its installed packages. Active workflows validate the
environment against Zijian's supplied Python/package versions and required
unpinned dependencies before doing substantive work; a stale or ABI-incompatible
environment fails loudly.

`hpc_contemporary_embedding_inventory.sh` is a CPU-only integration test for
the completed nine contemporary CSVs against Zijian's published embeddings.
The benchmark path is supplied at submission time, never embedded in the
wrapper. The job downloads Zenodo records 7409660 and 19498341 into scratch,
authenticates the nine source CSVs and three embedding archives, validates the
published cache, and runs the existing provenance-aware planner. Only compact
reports, the binary `reuse`/`regenerate` manifests and lists, and
`regenerate/prott5.fasta` are copied to
`$HOME/contemporary_embedding_inventory_results`; archive/cache bytes are not.
Scratch is removed after success, failure, or termination.
The wrapper must be submitted from a clean framework Git checkout unless a
complete `FRAMEWORK_COMMIT` is supplied. It checks out that exact commit in
scratch, so queue delay cannot silently change the reviewed code. Persistent
results are copied to a sibling staging directory and atomically renamed only
after success markers are present. Failed runs are published under a `.failed`
suffix without a success marker. Scratch cleanup remains unconditional even if
home-directory publication fails.

The common submission form is:

```bash
qsub -v BENCHMARK_DIR="$HOME/contemporary_cafa_benchmark_results/<run>/outputs" \
  hpc_jobs/active/hpc_contemporary_embedding_inventory.sh
```

Any role can be overridden with `BP_TRAINING_CSV`, `BP_VALIDATION_CSV`,
`BP_TEST_CSV`, and corresponding `CC_*`/`MF_*` variables. With no
`BENCHMARK_DIR`, all nine variables are required. Optional overrides are
`SOURCE_BENCHMARK_DIR`, `PUBLISHED_EMBEDDING_ARCHIVE_DIR`, `PFP_REFERENCE_DIR`, `ALIASES_FILE`,
`INVENTORY_CONFIG`, `INVENTORY_POLICY`, `REPORT_LEVEL`, `RESULTS_ROOT`, and the
complete 40-character `FRAMEWORK_COMMIT`.
The default is to download the canonical source CSVs and published embedding
archives during each job, as requested.

`hpc_contemporary_benchmark_reuse_plan.sh` runs the newer CSV-only binary reuse
planner. It stages the completed contemporary nine-CSV benchmark in scratch and,
by default, downloads and authenticates Zijian's canonical CAFA3 CSVs from
Zenodo record 7409660. It compares exact case-sensitive protein IDs and complete
sequences, producing only `reuse` and `regenerate` buckets. It does not download
or inspect embedding arrays. Results are published under
`$HOME/contemporary_benchmark_reuse_results`; scratch cleanup is unconditional.
Use `EMBEDDED_BENCHMARK_DIR` to replace the default download with a local set of
nine previously embedded CSVs. Optional names are `TARGET_BENCHMARK_NAME` and
`EMBEDDED_BENCHMARK_NAME`; `RESULTS_ROOT` and `FRAMEWORK_COMMIT` follow the same
rules as the inventory wrapper.

For the completed 2025-2026 run currently retained on Morecambe, submit:

```bash
qsub -v TARGET_BENCHMARK_DIR="$HOME/contemporary_cafa_benchmark_results/7065592_20260714_090900/outputs" \
  hpc_jobs/active/hpc_contemporary_benchmark_reuse_plan.sh
```

### Contemporary embedding generation and assembly

`hpc_contemporary_embedding_generation.sh` is the overnight continuation of
the CSV-only reuse plan. The exact benchmark and reuse-plan directories are
required command-line arguments; the wrapper does not guess which run to use:

```bash
qsub hpc_jobs/active/hpc_contemporary_embedding_generation.sh \
  --target-benchmark-dir "$HOME/contemporary_cafa_benchmark_results/7065592_20260714_090900/outputs" \
  --reuse-plan-dir "$HOME/contemporary_benchmark_reuse_results/7069671_20260715_054253/plan"
```

The job requests three GPU slots and preserves the established PFP parallel
layout: ProtT5 on GPU 0, PubMedBERT text on GPU 1, ESM-IF1 on GPU 2, and STRING
PPI extraction on CPU. It first runs a bounded preflight, then expands the same
scratch workspace to the complete `regenerate` partition. The PFP checkout is
pinned and disposable. A UniProt-only PPI denominator guard is made in a
separate, checksummed scratch copy. IF1 receives a scratch-local NumPy 1.26.4
overlay for compatibility with Biotite 0.38.0 and a checksummed extractor copy
that places fair-esm encoder inputs on the CUDA model device. Upstream PFP and
the primary author-pinned MMFP environment are not edited.

The wrapper downloads and authenticates Zijian's three published embedding
archives in scratch. Assembly uses published arrays only for planner-approved
`reuse` proteins and newly generated arrays only for `regenerate` proteins.
ProtT5 regeneration is mandatory for every regenerate protein. Missing text,
structure, or PPI source coverage is recorded and left absent so PFP retains
its existing zero-vector/mask behavior. PubMedBERT hidden states are reduced
incrementally to the exact `(768,)` CLS vector consumed by PFP, avoiding a
roughly 151 GiB raw text cache.

Successful results are published under
`$HOME/contemporary_embedding_generation_results` and contain:

```text
archives/contemporary_embedding_cache.tar.gz
archives/generated_prott5.tar.gz
archives/generated_text.tar.gz
archives/generated_structure.tar.gz
archives/generated_ppi.tar.gz
reports/assembly/assembly_summary.json
reports/assembly/embedding_assembly.tsv.gz
reports/reuse_plan/{RUN_COMPLETE,output_manifest,run_manifest,summary}.json
reports/archive_manifest.tsv
reports/if1_environment.json
reports/pfp_if1_compatibility.json
logs/{preflight,full}/
WORKFLOW_COMPLETE.json
```

Unpacked `.npy` files and large external inputs are never copied home. Failed
runs retain compact logs/reports and any generated modality archives completed
before the failure under a `.failed` result directory. Scratch cleanup remains
unconditional even if home publication fails. Optional overrides include
`TARGET_BENCHMARK_DIR`, `REUSE_PLAN_DIR`, `TEXT_CUTOFF_DATE`, `RESULTS_ROOT`,
`PREFLIGHT_PER_SPLIT`, the complete `FRAMEWORK_COMMIT`, and the pinned
`PFP_COMMIT`.

The completion JSON uses commit values already resolved by the host wrapper; it
does not invoke `git` inside the MMFP Python container.

### Contemporary embedding state and retries

The completed cache and exact benchmark inputs are preserved at:

```text
/SAN/bioinf/bmpfp/embeddings/contemporary/2025_01_to_2026_02_supervisor/
/SAN/bioinf/bmpfp/benchmarks/contemporary/2025_01_to_2026_02_supervisor/
```

Initialize the archive-backed state once:

```bash
qsub hpc_jobs/active/hpc_contemporary_embedding_state_initialize.sh
```

This CPU-only job verifies the archive, indexes its assembly report, binds the
state to the exact CSV/sequence/manifests/runtime contract, and writes the
state beneath the embedding release as `retry_state/`. It does not extract the
baseline arrays persistently and does not generate embeddings.

After inspecting the initialization report, retry one modality per job:

```bash
qsub hpc_jobs/active/hpc_contemporary_embedding_retry.sh --modality structure
qsub hpc_jobs/active/hpc_contemporary_embedding_retry.sh --modality text
qsub hpc_jobs/active/hpc_contemporary_embedding_retry.sh --modality ppi
```

The retry wrapper prefers frozen inputs under `/SAN/bioinf/bmpfp`: canonical
CAFA3 CSVs are copied into scratch before PFP's header normalization, the STRING
H5 is read directly from SAN, and the compressed STRING alias file is expanded
into scratch. Missing SAN inputs retain the original download fallback. Model
weights and AlphaFold structures acquired during retries are retained once in
the retry state's `source_cache/` for later attempts.

ProtT5 currently has complete coverage and does not need a retry. Each retry
requests only the pending protein/modality pairs, verifies regenerated control
arrays, merges valid successes into one SAN delta, publishes compact reports
under `$HOME/contemporary_embedding_retry_results`, and always deletes its
scratch directory. These commands are intentionally manual; the initializer
does not submit retries automatically.

## Homology-cluster identity array

### Final scratch-first runtime entrypoints

The normal no-persistent-input workflow now has two thin direct `qsub` entrypoints:

```text
hpc_jobs/active/hpc_homology_cluster_runtime_pilot.sh
hpc_jobs/active/hpc_homology_cluster_runtime_array.sh
```

Submit the recommended one-item 30% pilot with:

```bash
qsub hpc_jobs/active/hpc_homology_cluster_runtime_pilot.sh
```

Submit all six identities without a pilot prerequisite with:

```bash
qsub hpc_jobs/active/hpc_homology_cluster_runtime_array.sh
```

The array maps tasks `1-6` to `30, 25, 20, 15, 10, 5%` and uses `-tc 2` to avoid six
simultaneous copies of the roughly 170 GB compressed source collection. Node-local scratch is not
shared between array tasks, so every task with missing path overrides downloads its own inputs.
The 30% diagnostic pilot requests approximately `300G` scratch in total. UCL Grid Engine accounts
the consumable `tmem` and `tscratch` requests per SMP slot, so both wrappers request `tmem=8G` and
`tscratch=38G` across eight slots: 64 GB memory and 304 GB scratch per task. The non-consumable
`scratch0free=300G` remains a host-free-space threshold. The unsupported parser/MMseqs/publication
defaults are reduced globally to neutral `1x` bookkeeping; the pilot records rather than enforces
the resulting estimate. The six-task wrapper uses the same provisional per-task request, but do
not submit that array before
reviewing the pilot's disk report.

Both wrappers delegate to
`scripts/benchmark_generation/run_homology_cluster_runtime_hpc.sh`. That driver:

- uses the existing `mmfp` environment and clones the exact submitted framework commit;
- verifies the detached scratch checkout with host Git and passes that exact clean state into the
  minimal Singularity runtime, avoiding a redundant in-container Git dependency;
- defaults to `sprot-and-trembl`, while accepting `sprot-only` and `trembl-only` explicitly;
- prefers authenticated frozen inputs under `/SAN/bioinf/bmpfp`, stages any supplied source paths
  into task-owned scratch, and downloads only sources still missing;
- checks that mutable UniProt endpoints still mean release `2026_02`, while downloading GOA `234`
  from EBI's immutable historical URL and validating its pinned SHA-256 and embedded release
  metadata;
- downloads pinned MMseqs2 release `18-8cc5c` into scratch when `MMSEQS_BIN` is
  not supplied, and separately validates and records its binary-reported full
  Git commit `8cc5ce367b5638c4306c2d7cfc652dd099a4643f`;
- creates the checksum-bound frozen-input manifest and a clearly labelled automatic, non-blocking
  runtime attrition policy;
- calls the normal homology builder and strict publication validator, then writes an automatic
  review beside every task's five pickles and nine CSVs;
- samples allocated bytes in the job-owned scratch tree every 120 seconds and at named stage
  boundaries, writing `logs/disk_usage.tsv`, `logs/disk_usage_by_path.tsv`, and
  `logs/disk_usage_summary.tsv` with the measured peak;
- copies important results and logs atomically under
  `$HOME/homology_cluster_benchmark_results` by default; and
- always removes the task-owned scratch directory, including when home copy-back fails. A copy
  failure returns non-zero and leaves the Grid Engine `.o` log as the diagnostic.

Optional local inputs use the existing variables `UNIREF90_FASTA`, `IDMAPPING`,
`UNIPROT_SPROT_SEQUENCES`, `UNIPROT_TREMBL_SEQUENCES`, `GOA`, and `GO_OBO`. Each may be paired with
its existing `*_SHA256` variable. Missing paths are downloaded into scratch. `RESULTS_ROOT`,
`FRAMEWORK_REVISION`, `MMSEQS_BIN`, `UNIPROT_SOURCE_SCOPE`, `SPLIT_POLICY`, `SEED`, and `MIN_COUNT`
are also overridable. Large source bytes and MMseqs2 temporary data are never copied home.

The automatic runtime review proves the software/output contract; its permissive attrition bounds
are not a substitute for biological interpretation. Running the pilot first remains good practice,
but no pilot file, approval, or completion marker is required by the full runtime array.

### Pre-staged reviewed-input launchers

The older, more restrictive reviewed-input workflow has three components:

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

The pilot launcher constructs `qsub -t 1 -pe smp 2`; the full launcher constructs
`qsub -t 1-6 -pe smp 2`. `NSLOTS` is authoritative within each task and must equal the MMseqs2
thread count. Six runnable tasks can request up to 12 CPU slots, but scheduler concurrency remains
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

The full launcher in this older workflow rechecks approval evidence before constructing the command. The queued worker
rechecks launcher-time evidence hashes and reruns approval from its detached checkout before input
staging. Approval binds the reviewed attrition-policy hash and a common pilot run ID across the
marker, task context, and measurement evidence. Authorization reconstructs the pilot ratios and
requires the reviewed policy to accept them; diagnostic reports remain explicitly non-authorizing.
That strict reviewed-input path deliberately has no bypass. It is separate from the final
scratch-first runtime array above, whose pilot is optional. One array-wide attrition override is rejected because observed failures can
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
