# Validation Scripts

This directory contains validation workflows for checking benchmark-generation
code against historical reference artefacts.

## CAFA3 Historical Validation

Main runner:

```bash
bash scripts/validation/run_cafa3_historical_validation.sh
```

HPC wrapper:

```bash
qsub hpc_jobs/active/hpc_cafa3_historical_validation.sh
```

The HPC wrapper is intentionally thin. It follows the existing cluster-script
pattern:

1. create a job-specific scratch directory;
2. clone the framework into scratch;
3. call the real implementation under `scripts/validation/`;
4. copy reports and logs back to `$HOME/cafa3_historical_validation_reports/`;
5. remove the scratch directory.

Raw GOA, UniProt and GO downloads are kept in scratch only. They are not copied
back to home.

## Outputs

Each run copies a small report bundle to:

```text
~/cafa3_historical_validation_reports/<job-or-timestamp>/
```

Expected report files:

```text
cafa3_historical_validation_report.md
csv_comparison.tsv
pickle_comparison.tsv
protein_overlap.tsv
go_term_overlap.tsv
run_manifest.md
logs/
```

The scratch run directory contains the heavy raw/generated/reference files while
the job is running:

```text
raw/
generated/
reference/
reports/
logs/
```

By default scratch is deleted at the end of the run. Set `KEEP_SCRATCH=1` only
for debugging.

## Optional Inputs

The workflow always downloads the historical raw CAFA3 snapshots and the MMFP
reference split tarball.

DeepGOPlus reference pickle artefacts are optional because the exact public URL
may not always be discoverable automatically. To enable pickle comparison, set
one of:

```bash
export DEEPGOPLUS_PICKLES_DIR=/path/to/extracted/deepgoplus/cafa3/files
export DEEPGOPLUS_PICKLES_URL=https://example/path/to/deepgoplus-cafa3.tar.gz
```

If neither is set and no URL is discoverable from available local materials, the
workflow continues with CSV-only comparison and records the skip in the report.

## Python Environment

The validation runner uses the contemporary benchmark builder directly from:

```text
benchmark_builders/contemporary_cafa/src
```

The active Python environment must provide:

```text
numpy
pandas
```

The workflow does not modify PFP and does not install packages automatically.

## Scratch Policy

The cluster wrapper and main runner both clean scratch by default. This is
deliberate: raw database snapshots and extracted build artefacts are large and
should not be left behind after the reports have been copied home.
