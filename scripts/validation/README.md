# Validation Scripts

This directory contains validation workflows for checking benchmark-generation
code against historical reference artefacts.

## CAFA3 Historical Validation

There are now three historical validation workflows.

### DeepGOPlus Pickle-Generation Validation

This validates the recovered `cafa3_data.py` layer:

```text
Official CAFA3/DeepGOPlus files
    -> train_data.pkl
    -> test_data.pkl
    -> terms.pkl
    -> compare against released DeepGOPlus pickles
```

Main runner:

```bash
bash scripts/validation/run_cafa3_deepgoplus_pickle_generation_validation.sh
```

HPC wrapper:

```bash
qsub hpc_jobs/active/hpc_cafa3_deepgoplus_pickle_generation_validation.sh
```

By default the runner downloads the public CAFA archive candidates:

```text
https://deepgo.cbrc.kaust.edu.sa/data/data-cafa.tar.gz
https://deepgo.cbrc.kaust.edu.sa/data/deepgoplus-cafa.tar.gz
```

It can also use a local extracted copy:

```bash
export DEEPGOPLUS_CAFA_DIR="/path/to/data-cafa"
```

The runner copies reports/logs back to:

```text
~/cafa3_deepgoplus_pickle_generation_reports/<job-or-timestamp>/
```

Scratch is removed by default.

### DeepGOPlus/TEMPROT Validation

This is the preferred validation for the recovered PFP-facing pipeline:

```text
DeepGOPlus released pickles
    -> TEMPROT-style ontology CSV export
    -> PFP-compatible 9 CSVs
    -> compare against Zenodo 7409660
```

Main runner:

```bash
bash scripts/validation/run_cafa3_deepgoplus_validation.sh
```

HPC wrapper:

```bash
qsub hpc_jobs/active/hpc_cafa3_deepgoplus_validation.sh
```

This path downloads the DeepGOPlus CAFA archive by default. It tries the known
CAFA archive first and then the older DeepGOPlus-named archive as a fallback:

```text
https://deepgo.cbrc.kaust.edu.sa/data/data-cafa.tar.gz
https://deepgo.cbrc.kaust.edu.sa/data/deepgoplus-cafa.tar.gz
```

It can also use a local extracted copy:

```bash
export DEEPGOPLUS_PICKLES_DIR="/path/to/data-cafa"
```

The runner copies reports/logs back to:

```text
~/cafa3_deepgoplus_validation_reports/<job-or-timestamp>/
```

Scratch is removed by default.

### Raw Snapshot Audit

This is the heavier audit path that attempts to regenerate the historical
intermediates from raw 2017 UniProt/GOA/GO snapshots. It is useful for
methodology forensics, but it is not expected to be bit-for-bit identical unless
all official CAFA/DeepGOPlus intermediate curation steps are reproduced.

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

The active HPC wrapper uses the scratch-heavy GOA mode by default:

```bash
DECOMPRESS_GOA=1
USE_PIGZ=1
```

and requests `200G` of scratch. This keeps the decompressed `.gaf` files in
scratch for faster parsing, then removes them during cleanup.

## Outputs

Each run copies a small report bundle to:

```text
~/cafa3_historical_validation_reports/<job-or-timestamp>/
```

Expected report files:

```text
cafa3_historical_validation_report.md
cafa3_deepgoplus_pickle_generation_report.md
csv_comparison.tsv
pickle_comparison.tsv
pickle_generation_comparison.tsv
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

## GOA Parsing Performance

GOA files are kept compressed by default and streamed by the benchmark builder.
The builder filters GOA rows early against the loaded UniProt accession universe,
then evidence code, `NOT`, aspect, and taxon. It also logs progress during large
GAF parses.

Useful knobs:

```bash
# Default: use pigz for streaming gzip decompression if available.
export USE_PIGZ=1

# Disable pigz and use Python/gzip fallback.
export USE_PIGZ=0

# Default for direct/local runner use: keep .gaf.gz compressed on disk and stream it.
export DECOMPRESS_GOA=0

# HPC wrapper default: decompress GOA to .gaf in scratch before parsing. This
# uses more scratch but can help when gzip decompression is the limit.
export DECOMPRESS_GOA=1

# Default progress heartbeat.
export GOA_PROGRESS_INTERVAL=1000000
```

These settings affect scratch/runtime only. They do not change the benchmark
definition.

## Optional Inputs

The workflow always downloads the historical raw CAFA3 snapshots and the nine
canonical CAFA3 CSV reference files from Zenodo record
`https://zenodo.org/records/7409660`. These are the CSV benchmark interface that
PFP ultimately consumes.

Zijian's Zenodo record `https://zenodo.org/records/19498341` contains MMFP/PFP
artefacts generated from that benchmark, but it is not used as the canonical CSV
comparison source here.

DeepGOPlus reference pickle artefacts are optional. By default the workflow tries
the public CAFA archives:

```text
https://deepgo.cbrc.kaust.edu.sa/data/data-cafa.tar.gz
https://deepgo.cbrc.kaust.edu.sa/data/deepgoplus-cafa.tar.gz
```

To override the default, set one of:

```bash
export DEEPGOPLUS_PICKLES_DIR=/path/to/extracted/deepgoplus/cafa3/files
export DEEPGOPLUS_PICKLES_URL=https://deepgo.cbrc.kaust.edu.sa/data/data-cafa.tar.gz
```

If the pickle archive cannot be downloaded or extracted, the workflow continues
with CSV-only comparison and records the skip in the report.

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
