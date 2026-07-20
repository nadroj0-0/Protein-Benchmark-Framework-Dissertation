# Benchmark-Agnostic PFP Execution

This directory is the additive compatibility layer between a completed
nine-CSV benchmark, its embedding cache, and the immutable upstream PFP
checkout. It does not fork PFP and it does not generate embeddings.

The central finding behind this layer is that upstream PFP is already mostly
benchmark-agnostic:

- `scripts/prepare_cafa3_data.py` accepts any directory containing the nine
  expected CSV filenames;
- `train.py` accepts explicit data, ontology, output, seed and text-cache
  paths; and
- PFP's dataset and CAFA evaluator already handle missing non-sequence
  modalities through zero vectors and binary masks.

The unsafe parts were orchestration and validation. Upstream only warns about
split overlap, silently coerces malformed arrays, can skip stale output, and
can finish without CAFA metrics. The framework validates every boundary before
calling PFP and checks every required output afterwards.

## Supported operations

`run_pfp_benchmark.sh` has three execution modes:

- `prepare-only`: validate and materialize a nine-CSV benchmark;
- `eval-only`: evaluate an existing standard PFP checkpoint tree; and
- `train-eval`: train fresh checkpoints and evaluate them.

It has two modality modes:

- `full`: use all four modalities. Sequence coverage must be 100%. Missing
  text, structure and PPI arrays remain absent so PFP applies its original
  zero-vector/mask behavior. Every present array must be valid.
- `sequence-only`: require 100% ProtT5 coverage and deliberately mask text,
  structure and PPI. This preserves PFP's own four-branch masked control; it is
  not a newly designed standalone sequence architecture.

The full and sequence-only results should both be reported. Full mode holds
Zijian's behavior constant, while sequence-only provides a clean comparison
when modality coverage differs substantially between historical and current
proteins.

## Validation contract

Before preparation, the wrapper requires:

- exactly the nine `bp/cc/mf-{training,validation,test}.csv` files;
- safe unique protein IDs and nonempty valid sequences;
- binary GO labels and configurable all-zero-row policy;
- exact GO column identity and order across all three splits per ontology;
- configurable global or per-ontology protein-ID and exact-sequence overlap policy;
- every GO term to exist in the supplied OBO and belong to the correct
  namespace; and
- no duplicate rows unless a policy explicitly permits them.

The CAFA3 config alone accepts the published singular `protein` header. It
normalizes that header in disposable working storage and records the alias;
the source CSV is never edited.

After unmodified upstream preparation, the wrapper checks raw-to-materialized
protein order, sequences, sparse labels and ordered GO terms. Terms with zero
training positives are retained and reported, which is required when a frozen
homology term universe moves all positive examples into another split.

Before training or evaluation, every present target array is loaded with
pickle disabled and checked for its exact one-dimensional shape, numeric dtype
and finite values. Missing non-sequence arrays are listed in the issues TSV and
reported per aspect/split. Missing or invalid sequence arrays are fatal.

Each aspect is trained separately with upstream `train.py --single`, so a
failure cannot be swallowed by PFP's multi-experiment exception handler. Fresh
checkpoints are then evaluated again by the framework through PFP's normal
model/dataset code but a strict cafaeval call: `ia`, `norm=cafa`, `prop=max`,
and `no_orphans=false` may not be silently discarded. A run is complete only when every
selected aspect has a checkpoint and finite, physically valid `cafa_fmax`,
`cafa_wfmax`, and `cafa_smin` values.

`--capture-predictions` is opt-in. It observes the arrays from that same CAFA
inference call and atomically publishes a standalone sensitivity bundle. The
bundle includes exact truth and score arrays, IA bytes, benchmark and prepared
data fingerprints, validation reports, config/OBO hashes, seed, and code
revisions. Capture requires those provenance reports; ordinary evaluation does
not.

## Frozen model behavior

All configs deliberately enforce Zijian's model settings:

- ProtT5 sequence input;
- gated bilinear fusion;
- late fusion in hybrid mode;
- hidden dimension 512;
- model dropout 0.4;
- modality dropout 0.1; and
- auxiliary loss weight 0.8.

Changing these values is a model experiment, not benchmark compatibility, and
is rejected by this layer.

Every benchmark and ontology requires fresh checkpoints because its ordered GO
term output space can differ. Evaluation-only additionally requires
`--reference-data-dir`: freshly prepared names, labels, sequences and term
order must match the data paired with the checkpoint. This closes the otherwise
silent same-dimension/different-term-order failure mode.
For the CAFA3 v1.5 control, the CAFA3 config additionally requires
`--reference-source-archive`; the supplied directory must match the 30 relevant
members of the catalogued author archive byte for byte.

The direct runner rejects a dirty framework checkout by default. The
`--allow-dirty-framework` escape hatch exists only for tiny development
fixtures and marks the reported framework revision as dirty; it must not be
used for dissertation model runs. The HPC wrapper has no such escape hatch and
always clones a pinned clean framework commit.

## Policy configs

- `configs/pfp_benchmark_run.cafa3.json`: published CAFA3 compatibility,
  including the legacy MF header alias.
- `configs/pfp_benchmark_run.temporal.json`: contemporary temporal benchmark;
  protein IDs are globally disjoint between development and test, while exact
  sequences are disjoint within each ontology. Cross-ontology sequence overlap
  is reported but allowed because PFP trains BP, CC and MF separately.
- `configs/pfp_benchmark_run.homology.json`: protein IDs and exact sequences are
  globally disjoint across all three splits. Domain-owned homology validation
  evidence is mandatory.

The homology config does not replace the homology builder's cluster validation.
Pass its `validation_report.json`, `output_manifest.json` and completion marker
with repeated `--benchmark-evidence` options. The generic runner requires a
passed validation report, verifies the completion marker binds the manifest,
then verifies that the manifest binds both that report and the exact nine CSV
bytes selected for training.
It also invokes the homology builder's own full publication validator against
the benchmark directory, covering cluster/split metadata, scientific
fingerprint derivation, manifests, term-universe recount, and attrition policy.
The PFP boundary additionally requires `benchmark_scope=dissertation-production`,
`production_eligible=true`, and `fixture_mode=false`; a valid fixture or
diagnostic pilot is evidence, not a trainable dissertation benchmark.

Embedding assembly state is supplied with repeated `--embedding-evidence`.
For a required dissertation gate, pass `coverage.json`, `contract.json`,
`targets.tsv`, and `pair_status.tsv` together. The runner verifies the state
contract hash, exact nine CSV hashes, protein IDs and sequence hashes, modality
policy, accepted IDs, per-pair accepted embedding hashes, and observed cache
counts. It also hashes every valid array before execution and repeats that
check after execution, so the final report identifies the exact cache bytes
used and rejects both pre-run substitution and in-run mutation.

## Local smoke usage

Preparation only performs no model training:

```bash
bash scripts/model_execution/run_pfp_benchmark.sh \
  --benchmark-id contemporary-2025-2026-supervisor \
  --benchmark-dir /absolute/path/to/nine/csvs \
  --obo-file /absolute/path/to/go-basic.obo \
  --pfp-root /absolute/path/to/PFP \
  --work-dir /absolute/path/to/new/work \
  --output-dir /absolute/path/to/new/output \
  --config configs/pfp_benchmark_run.temporal.json \
  --execution-mode prepare-only
```

The local test suite uses tiny fixtures and fake PFP boundaries. It never
downloads data or starts real training:

```bash
PYTHONDONTWRITEBYTECODE=1 \
  python3 -m unittest discover -s scripts/model_execution/tests -v
```

## CAFA3 v1.5 acceptance rung

The orchestration-only control uses Zijian's published CSVs, embeddings,
checkpoints and their paired prepared data. The `--reference-data-dir` must be
the top-level `data/` directory extracted directly from Zijian's published
`mmfp_data_splits.tar.gz` (Zenodo record 19498341), not data produced by an
earlier local run. The frozen archive is catalogued at:

```text
/SAN/bioinf/bmpfp/reference_artifacts/zijian_mmfp/mmfp_data_splits.tar.gz
```

The CAFA3 policy verifies the published archive's catalogued size and MD5, then
checks that all 30 reference files are byte-identical to its `data/` members.
The preparation report records those member SHA-256 values and a combined
reference fingerprint, making both the comparand and its published provenance
auditable.

```bash
qsub hpc_jobs/active/hpc_pfp_benchmark.sh \
  --benchmark-id cafa3-v1.5-published-artifacts \
  --benchmark-dir /path/to/published/nine-csvs \
  --embedding-cache-root /path/to/published/embedding_cache \
  --checkpoint-root /path/to/published/results/full_model \
  --reference-data-dir /path/to/directly-extracted-mmfp-data-splits/data \
  --reference-source-archive /SAN/bioinf/bmpfp/reference_artifacts/zijian_mmfp/mmfp_data_splits.tar.gz \
  --obo-file /path/to/cafa3/go.obo \
  --ia-file-dir /path/to/published/prepared/data \
  --results-root /SAN/bioinf/bmpfp/model_runs/cafa3_v1_5 \
  --config "$PWD/configs/pfp_benchmark_run.cafa3.json" \
  --execution-mode eval-only \
  --expected-metrics "$PWD/configs/pfp_cafa3_published_metrics.json" \
  --require-reference-match
```

This isolates the new wrapper from embedding regeneration. Only after v1.5
matches should the same runner be used with regenerated CAFA3 embeddings (v2).
The already observed fresh-retraining variation, approximately 0.003 Fmax, is
recorded separately from the tighter evaluation-only tolerance.
The CAFA3 policy requires and hash-records all three published
`BPO/CCO/MFO_ia.txt` files. It never silently recomputes IA for this acceptance
rung.

`42` is the default seed, not a forced value. Both the local runner and HPC
wrapper accept `--seed N`, forward it to training and evaluation, and record it
in the final report.

## Optional prediction artifacts

Add `--capture-predictions` to either runner when a separate label-cohort
sensitivity analysis is planned. The option is off by default. It observes the
prediction and truth arrays from PFP's existing strict CAFA evaluation pass,
records their content hashes together with protein/term order, checkpoint and
IA hashes, and publishes them under `evaluation/prediction_artifacts/`. It does
not change canonical metrics and does not perform a second inference pass.

The artifact is intentionally larger than the ordinary metric report. Capture
it only for runs that will be analysed, then pass its manifest to
`scripts/diagnostics/evaluate_pfp_label_sensitivity.py`. That analysis is
separately staged and can never overwrite the canonical run.

## Existing embedding-state evidence upgrade

States initialized before per-array evidence hashes were introduced must not be
re-initialized with a newer framework commit. After every retry job targeting
the state has finished and final coverage has been read, submit the dedicated
additive upgrade instead:

```bash
qsub hpc_jobs/active/hpc_embedding_state_evidence_upgrade.sh \
  --state-root /SAN/.../retry_state \
  --confirm-retries-finished
```

The command preserves the existing contract, accepted membership, failure
ledger and cumulative cache. It verifies the contracted baseline archive,
hashes every accepted baseline and delta array under the state lock, refreshes
`pair_status.tsv`, and fails if accepted counts differ before and after. Expect
it to take longer than the original initialization because it reads every
accepted array. Do not submit it while a retry generation job is still running;
the lock can serialize merges but cannot detect a job that has not reached its
merge step yet.

## Contemporary training

First consolidate the archive-backed retry state after all retry jobs finish:

```bash
qsub hpc_jobs/active/hpc_finalize_contemporary_embedding_state.sh \
  --state-root /SAN/.../embeddings/contemporary/.../retry_state \
  --benchmark-dir /SAN/.../benchmarks/contemporary/... \
  --obo-file /SAN/.../frozen_inputs/ontology/2025-02-06/go-basic.obo \
  --final-root /SAN/.../embeddings/contemporary/.../finalized_pfp_cache \
  --confirm-retries-finished \
  --retire-source-embeddings
```

Then use the completed benchmark, the frozen t0 benchmark ontology that defines
its term universe, the consolidated archive and its self-contained evidence.
Paths below are examples, not defaults:

```bash
qsub hpc_jobs/active/hpc_pfp_benchmark.sh \
  --benchmark-id contemporary-2025-01-to-2026-02-supervisor \
  --benchmark-dir /SAN/.../benchmarks/contemporary/... \
  --benchmark-evidence /SAN/.../benchmarks/contemporary/.../build_manifest.json \
  --embedding-cache-archive /SAN/.../finalized_pfp_cache/contemporary_embedding_cache.tar.gz \
  --embedding-evidence /SAN/.../finalized_pfp_cache/evidence/coverage.json \
  --embedding-evidence /SAN/.../finalized_pfp_cache/evidence/contract.json \
  --embedding-evidence /SAN/.../finalized_pfp_cache/evidence/targets.tsv \
  --embedding-evidence /SAN/.../finalized_pfp_cache/evidence/pair_status.tsv \
  --require-embedding-evidence \
  --obo-file /SAN/.../frozen_inputs/ontology/2025-02-06/go-basic.obo \
  --results-root /SAN/bioinf/bmpfp/model_runs/contemporary/full \
  --config "$PWD/configs/pfp_benchmark_run.temporal.json" \
  --execution-mode train-eval \
  --modality-mode full
```

Repeat with a separate result root and `--modality-mode sequence-only` for the
coverage-robust control. Never point two jobs at the same result directory.
The archive is safely extracted into job-owned scratch before the same strict
cache/evidence validator runs; directory-backed cache input remains supported.

## Homology integration

No further PFP code changes are needed for homology. Each completed threshold
supplies its own nine CSVs and domain validation artifacts to the same wrapper:

```bash
qsub hpc_jobs/active/hpc_pfp_benchmark.sh \
  --benchmark-id homology-identity-30 \
  --benchmark-dir /SAN/.../homology/identity_30 \
  --benchmark-evidence /SAN/.../identity_30/validation_report.json \
  --benchmark-evidence /SAN/.../identity_30/output_manifest.json \
  --benchmark-evidence /SAN/.../identity_30/RUN_COMPLETE.json \
  --embedding-cache-root /SAN/.../embeddings/homology/identity_30/cache \
  --embedding-evidence /SAN/.../embeddings/homology/identity_30/coverage.json \
  --embedding-evidence /SAN/.../embeddings/homology/identity_30/contract.json \
  --embedding-evidence /SAN/.../embeddings/homology/identity_30/targets.tsv \
  --embedding-evidence /SAN/.../embeddings/homology/identity_30/pair_status.tsv \
  --require-embedding-evidence \
  --obo-file /SAN/.../go-basic.obo \
  --results-root /SAN/bioinf/bmpfp/model_runs/homology/identity_30 \
  --config "$PWD/configs/pfp_benchmark_run.homology.json" \
  --execution-mode train-eval
```

Thresholds should be submitted as independent jobs, or through a separate
task-indexed launcher that maps `SGE_TASK_ID` to one threshold, after one pilot
has validated the paths and result contract. This model wrapper performs one
explicit benchmark run; it does not map array indices, submit, or recursively
schedule jobs itself.

## Metrics and interpretation

Ordinary CAFA Fmax is the primary cross-benchmark metric. Weighted Fmax and
Smin are still generated and required, but they are secondary
benchmark-specific metrics because their information-accretion values depend
on each benchmark's training labels and ontology snapshot. Coverage for every
modality must be reported beside every full-model result.
