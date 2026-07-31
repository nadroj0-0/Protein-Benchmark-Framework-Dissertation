# Diagnostics

## Audit a PFP working copy

`audit_pfp_working_copy.py` compares an existing local PFP working directory
with a fresh temporary clone of `psipred/PFP`. It does not modify the supplied
working directory, upload anything, or retain the public clone.

Activate the Python environment that was used for PFP, then run:

```bash
python scripts/diagnostics/audit_pfp_working_copy.py /path/to/working/PFP \
  > pfp_working_copy_audit.md
```

The redirection is optional and is the only persistent output. Without it, the
Markdown report is printed to the terminal. A non-zero exit means that the
working copy differs materially from the current public release.

The report separates:

- public tracked files modified or missing locally;
- local Git-tracked files absent from the public release;
- names of local untracked entries, without reading or printing their contents;
- the active Python environment used to run the audit, including a complete
  `pip freeze --all`.

Local-only files are evidence of local development, not automatic evidence
that those files should have been included in the published repository.

## Benchmark-agnostic label-space audit

`audit_pfp_label_space.py` audits any benchmark implementing the exact PFP
nine-CSV contract. It is not CAFA3-specific: use the same command for the
contemporary benchmark, each homology threshold, and future benchmark builds.
It validates the CSVs against the supplied ontology, including ancestor closure
and root connectivity, reports root-only targets, annotation depth, term support
and label concentration, and can independently verify one or more prepared PFP
data directories against the source CSVs. Both Unicode and the object-dtype name
arrays emitted by upstream PFP are accepted, but every loaded name must be a
plain string.

```bash
python scripts/diagnostics/audit_pfp_label_space.py \
  --benchmark-id contemporary-2025-2026 \
  --benchmark-dir /absolute/path/to/nine-csvs \
  --obo-file /absolute/path/to/go.obo \
  --config configs/pfp_benchmark_run.temporal.json \
  --prepared-data framework=/absolute/path/to/prepared_data \
  --output-dir /absolute/path/to/new/label-audit
```

For CAFA3, pass the direct published nine CSVs and the directly extracted
author-prepared data as separate evidence. The legacy singular `protein`
header is accepted only through the CAFA3 config and is recorded as an alias;
the source CSV is never edited. High root-only prevalence is a reported
finding, not a hard failure threshold.

Compare any number of completed audits without assuming matching protein or
term universes:

```bash
python scripts/diagnostics/compare_pfp_label_audits.py \
  --report /path/to/cafa3/label_space_audit.json \
  --report /path/to/contemporary/label_space_audit.json \
  --report /path/to/homology-30/label_space_audit.json \
  --output-dir /absolute/path/to/new/label-audit-comparison
```

Every audit is staged and published only after validation. Its manifests bind
the exact CSV, OBO, policy, optional IA files and optional prepared-data
evidence. Input files are hashed before parsing and checked again afterward.

## Root-only evaluation sensitivity

Prediction sensitivity is deliberately separate from canonical PFP results.
First opt in while evaluating a checkpoint by adding `--capture-predictions`
to `run_pfp_benchmark.sh` or `hpc_pfp_benchmark.sh`. This observes the arrays
already produced by PFP's normal CAFA evaluation; it does not rerun inference.
The completed run publishes compressed prediction/truth arrays, protein and GO
term order, checkpoint hashes, the exact IA files, both code revisions, and
copies of the preparation and embedding-validation reports under
`evaluation/prediction_artifacts/`.

Then run the standalone analysis:

```bash
python scripts/diagnostics/evaluate_pfp_label_sensitivity.py \
  --prediction-manifest /path/to/prediction_artifacts/prediction_artifact_manifest.json \
  --obo-file /absolute/path/to/the-same-go.obo \
  --output-dir /absolute/path/to/new/root-exclusion-sensitivity
```

The analysis first reproduces the canonical strict cafaeval result from the
captured artifact and fails if it drifts. It then reports:

- strict CAFA metrics after excluding targets with no positive non-root term;
- both a re-optimized Fmax and a result fixed at the canonical threshold;
- a strict-cafaeval root-only prediction baseline; and
- a clearly labelled flat non-root diagnostic with no GO propagation.

It never retrains a model or overwrites canonical output. Captured benchmark
rows, cafaeval-evaluable targets and all-zero rows are reported separately.
The original-threshold cohort result avoids retuning after changing the test
cohort. When modes are compared, each mode's original threshold is shown
explicitly; this is not described as a shared-threshold comparison.

After running the available modality modes, compare every additive/full mode
with the sequence-only baseline using repeated `--report` arguments:

```bash
python scripts/diagnostics/compare_pfp_label_sensitivity.py \
  --report /path/to/full/root_exclusion_sensitivity.json \
  --report /path/to/sequence-only/root_exclusion_sensitivity.json \
  --report /path/to/sequence-text/root_exclusion_sensitivity.json \
  --report /path/to/sequence-structure/root_exclusion_sensitivity.json \
  --report /path/to/sequence-ppi/root_exclusion_sensitivity.json \
  --output-dir /absolute/path/to/new/sensitivity-comparison
```

The comparator calculates each mode-minus-sequence delta only when the benchmark
fingerprint, source CSVs, seed, config, OBO, PFP revision, protein order, truth,
GO-term order, IA bytes, sequence-embedding content and finalized embedding
evidence match. Non-evaluable aspects remain visible with an explicit status.
It does not treat Fmax values from different benchmark label spaces as directly
comparable model rankings.

## IA and Xu-specificity evaluation

Daniel Buchan's proposed information-content analysis can be run directly from
the same immutable prediction artifacts. The analyzer reads the exact
information-accretion file used by the canonical evaluation. It can also
calculate Xu et al.'s topology-only semantic totipotency from the exact frozen
`is_a + part_of` OBO graph. IA and Xu are always separate panels: IA remains
the canonical CAFA weight; Xu is only a specificity stratifier.

```bash
python scripts/diagnostics/evaluate_pfp_information_content.py \
  --prediction-manifest /path/to/prediction_artifacts/prediction_artifact_manifest.json \
  --obo /path/to/frozen-go.obo \
  --specificity-measure all_separate \
  --specificity-measure xu_neglog_totipotency \
  --positive-bins 4 \
  --bootstrap-replicates 10000 \
  --bootstrap-seed 0 \
  --output-dir /absolute/path/to/new/specificity-analysis
```

The output names Jaccard set agreement explicitly rather than calling it
ordinary accuracy. Unweighted protein-centric metrics are primary; IA-weighted
diagnostics are reported alongside them. Roots are excluded from flat bins.
Raw Xu `T` is lower for more specific terms. `-log2(T)` is retained only as an
exploratory display transform and is not attributed to Xu et al. or called IA.
Bin-specific optima are descriptive; the common fixed result is visibly
labelled `descriptive_test_oracle_fixed`. No binwise result is canonical CAFA
`Smin`.

## Policy-bound temporal cohort states

For Zijian's published CAFA3 benchmark, use the organizer classifications
instead of reconstructing organizer-private historical state from a later GOA
snapshot:

```bash
python scripts/diagnostics/build_cafa3_knowledge_state_census.py \
  --published-csv-dir /path/to/zijian-nine-csvs \
  --official-cafa-archive /path/to/benchmark20171115.tar \
  --output-dir /path/to/new/cafa3-knowledge-census
```

This classifies test rows with CAFA3's official `type1` (no-knowledge) and
`type2` (limited-knowledge) lists, audits their alignment with the published
CSVs, retains the same state labels for organizer `too_few` targets while
marking them separately, and cross-tabulates each state against root-only/non-root
observed truth.

`benchmark20171115.tar` is the organizer benchmark bundled inside the CAFA3
report's public Figshare supplement (`10.6084/m9.figshare.8135393.v3`). Do not
substitute DeepGOPlus's processed `data-cafa.tar.gz`; it contains ground truth
but not the organizer `type1`/`type2` target lists required for this census.
The official lists apply only to challenge targets, so training and validation
receive root-state counts but no invented CAFA knowledge-state label.

The accepted contemporary `supervisor` benchmark is global qualifying
no-knowledge, so its existing test predictions have no known-protein comparator.
`build_temporal_annotation_ledger.py` therefore builds a new Layer-B state
artifact from independently prepared historical inputs. It requires direct
terms and propagated closures for both timepoints plus explicit protein
presence. This prevents an absent protein being mislabelled as no-knowledge and
enforces closure-before-difference.

Every annotation input is a long-form TSV with the exact header
`protein_id<TAB>aspect<TAB>go_term`. Protein scope and presence files have the
single header `protein_id`.

```bash
python scripts/diagnostics/build_temporal_annotation_ledger.py \
  --t0-direct-annotations /path/to/direct_t0.tsv \
  --t1-direct-annotations /path/to/direct_t1.tsv \
  --t0-closure-annotations /path/to/closure_t0.tsv \
  --t1-closure-annotations /path/to/closure_t1.tsv \
  --t0-protein-presence /path/to/t0_present.tsv \
  --t1-protein-presence /path/to/t1_present.tsv \
  --exposure-table /path/to/development_exposure.tsv \
  --protein-scope /path/to/protein_ids.tsv \
  --t0-snapshot 2025-03-08 \
  --t1-snapshot 2026-06-17 \
  --evidence-policy-id supervisor_snapshot_membership \
  --graph-policy-id cafa_narrow_is_a_part_of \
  --relationship is_a \
  --relationship part_of \
  --benchmark-id contemporary-2025-2026 \
  --output-dir /absolute/path/to/new/temporal-ledger
```

The tool labels `no_qualifying`, `cross_ontology_known`,
`same_aspect_partial`, `root_only` and `unknown`; it does not conflate
same-aspect partial knowledge with CAFA limited knowledge. It records direct
and closure terms, then forms retained-known, gained and lost closure sets.
It deliberately does not parse raw GAF, resolve GO IDs, filter evidence or
perform propagation. Those Layer-A operations must be completed and
hash-bound upstream rather than inferred from filenames or reverse-engineered
from propagated CSV truth.

If any scoped protein has qualifying `t0` knowledge, `--exposure-table` is
mandatory. It records train/validation membership by ID and exact sequence,
declared homology-cluster overlap, modality availability and feature temporal
policy. This prevents a seen-protein annotation-extension sensitivity from
being described as unseen generalization.

After predictions exist for the same protein scope, evaluate the cohorts and
retained/gained partitions with:

```bash
python scripts/diagnostics/evaluate_pfp_knowledge_cohorts.py \
  --prediction-manifest /path/to/prediction_artifact_manifest.json \
  --temporal-ledger-dir /path/to/temporal-ledger \
  --truth-graph-policy-id cafa_narrow_is_a_part_of \
  --bootstrap-replicates 10000 \
  --output-dir /absolute/path/to/new/knowledge-cohort-analysis
```

On the accepted `supervisor` test this deliberately returns
`not_evaluable_empty_cohort` for known comparators. On a separately generated
partial cohort it reports retained-known recovery, acquisition-conditioned and
deployment-like gained-term performance, and the `t0` annotation-copy
baseline. These are flat diagnostics. They become strict unseen results only
after a separately authorized disjoint retraining design.

## Validation prediction capture foundation

`evaluate_pfp_checkpoints.py` keeps `--evaluation-split test` as its default.
An explicit `--evaluation-split valid` evaluates the saved checkpoint on the
prepared validation split and records `evaluation_split=valid` in the result
and prediction-artifact manifest. Validation artifacts must use a separate new
artifact directory. Existing test-only sensitivity analyzers reject them; they
are intended as inputs to a separately specified calibration analysis.

## Post-selection validation calibration

`calibrate_pfp_predictions.py` fits only on a captured validation artifact and
then evaluates temporal test transport once without refitting:

```bash
python scripts/diagnostics/calibrate_pfp_predictions.py \
  --validation-prediction-manifest /path/to/valid/prediction_artifact_manifest.json \
  --test-prediction-manifest /path/to/test/prediction_artifact_manifest.json \
  --obo /path/to/frozen-go.obo \
  --positive-ia-bins 4 \
  --output-dir /absolute/path/to/new/post-selection-calibration
```

The fitted quantity is the estimated probability of membership in the
benchmark-observed qualifying propagated `t1` label set. It is not biological
truth and is never a p-value. The model has a positive-slope Platt backbone,
regularized IA-bin and supported term intercepts, and the deterministic
fallback order `term_shrinkage -> aspect_mode_ia_bin -> aspect_mode_platt`.
It preserves raw, post-propagation and calibrated values, reliability and
temporal-shift metrics, support, fallbacks and a hierarchy audit.

This first implementation is explicitly
`post_selection_validation_calibration`: the same validation population
previously influenced checkpoint selection and early stopping. A primary
independent calibration claim requires a disjoint `V_select`/`V_cal` design
and usually retraining. Per-prediction uncertainty intervals remain a separate
publication gate; the point estimates must not be wrapped into a user-facing
tool until a protein-cluster-aware interval method has been validated.

Canonical Fmax, weighted Fmax and Smin can be compared independently from the
completed PFP run reports:

```bash
python scripts/diagnostics/compare_pfp_modality_runs.py \
  --run-report /path/to/sequence-only/reports/run_report.json \
  --run-report /path/to/sequence-text/reports/run_report.json \
  --run-report /path/to/sequence-structure/reports/run_report.json \
  --run-report /path/to/sequence-ppi/reports/run_report.json \
  --run-report /path/to/full/reports/run_report.json \
  --prediction-manifest sequence-only=/path/to/sequence-only/prediction_artifact_manifest.json \
  --prediction-manifest sequence-text=/path/to/sequence-text/prediction_artifact_manifest.json \
  --prediction-manifest sequence-structure=/path/to/sequence-structure/prediction_artifact_manifest.json \
  --prediction-manifest sequence-ppi=/path/to/sequence-ppi/prediction_artifact_manifest.json \
  --prediction-manifest full=/path/to/full/prediction_artifact_manifest.json \
  --output-dir /absolute/path/to/new/modality-comparison
```

Canonical reports must be from `train-eval` runs. A prediction capture may come
from the same run or a later `eval-only` replay, but its exact checkpoint hash,
canonical metrics, prepared benchmark, configuration, IA and active embedding
content must bind back to the canonical retraining run. Framework commit drift
is rejected unless it has been audited and explicitly acknowledged with
`--allow-framework-commit-drift`.
