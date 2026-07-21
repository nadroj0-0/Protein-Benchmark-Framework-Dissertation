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
