# PFP Embedding Inventory

This package inventories an existing PFP embedding cache against a nine-CSV
benchmark. It reports whether each protein/modality pair is present and whether
the available provenance is sufficient for reuse. It never generates or copies
embeddings.

## Policies

`paper-faithful` is the submitted default. `maximize-coverage` remains an
explicit diagnostic capability.

Sequence reuse requires exact sequence identity. Text, structure, and PPI
reuse additionally require authenticated benchmark, archive, cache, and PFP
evidence. Missing or uncertain provenance is reported as regeneration or
manual review rather than assumed safe.

The source benchmark and target benchmark have separate contracts. This is
important because a published CAFA3 cache can be a valid source even when its
population and split rules differ from the new target benchmark.

## Run

Direct package usage:

```bash
python -m pip install -e embedding_inventory
python scripts/verification/inventory_embeddings.py --help
```

The complete portable wrapper can stage authenticated published inputs:

```bash
bash scripts/verification/run_contemporary_embedding_inventory.sh \
  --benchmark-dir /absolute/path/to/nine_csvs \
  --work-dir /absolute/path/to/new_work_dir \
  --output-dir /absolute/path/to/new_output_dir \
  --policy paper-faithful \
  --artifact-catalog /absolute/path/to/artifact_paths.tsv
```

The wrapper pins the external PFP reference commit. Explicit source benchmark,
archive, and PFP checkout paths take precedence over catalogue or download
fallbacks.

## Output

The publication contains deterministic machine-readable summaries, per-pair
decisions and reasons, modality ID lists, a generation FASTA, checksums, and a
completion marker. Reports distinguish observed bytes from policy decisions.

Inventory output is advisory to generation. The separate
`benchmark_reuse_planner` establishes the coarse exact-ID/exact-sequence
partition, and `resolve_embedding_reuse_sources.py` authenticates actual source
archives for pair-resolved NK+LK and homology assembly.

All outputs must be written outside the source checkout.
