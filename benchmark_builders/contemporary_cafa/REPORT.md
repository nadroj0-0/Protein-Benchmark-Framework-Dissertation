# Implementation Report

## Scope

Version 0.2.0 hardens the raw-snapshot front end needed for the full 2025 to
2026 benchmark. PFP remains unchanged. The previously validated historical
CAFA-file to DeepGOPlus-pickle to TEMPROT/PFP-CSV paths remain intact.

## Corrected temporal methodology

- Test candidates originate from the `t0` target population, never the `t1`
  population.
- Every test row uses its `t0` sequence.
- Every protein with a qualifying `t0` annotation is excluded before sequence
  joins, term filtering and train/validation splitting.
- GAF annotation dates are retained during streaming and `t1` rows dated on or
  before the `t0` cutoff are removed as backfill.
- UniProt primary and secondary accessions provide the cross-release identity
  map. Missing, ambiguous and many-to-one mappings are excluded and reported.
- Sequence changes default to exclusion; `use-t0` and `error` are explicit CLI
  alternatives.
- GO IDs are resolved against their matching source ontology and then mapped
  into the frozen `t0` prediction graph.
- Final CAFA3 protein-binding-only MFO targets are removed under the default
  profile policy.

## Policy flexibility

Named profiles separate:

1. reconstructed CAFA3 policy;
2. contemporary CAFA3-style broad training policy;
3. the supervisor's same-organism and extended-evidence policy.

Training and target taxon scopes are independent. Training and target reviewed
status are also independent, so a Swiss-Prot training population can coexist
with a Swiss-Prot plus TrEMBL target population.

## Reproducibility and QC

Snapshot protein rows, annotations and GO columns are sorted before the seeded
DeepGOPlus split. Repeated fixture builds produce byte-identical CSV and pickle
artefacts.

Strict production QC checks temporal membership, annotation leakage, protein
and exact-sequence overlap, CSV schemas, duplicate IDs, binary labels and empty
ontology/split outputs. The builder writes input/output checksums, a complete
policy/environment manifest, per-target flow, exclusion reasons, evidence and
taxon summaries, and aggregate statistics.

## Operations

The reusable runner is:

```text
scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
```

The UCL Grid Engine wrapper is:

```text
hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh
```

The wrapper stages frozen inputs in node-local scratch, runs the reusable shell
entrypoint, copies the complete run directory to durable storage, and clears
scratch through an EXIT/SIGINT/SIGTERM trap. It does not run PFP training or
evaluation.

## Remaining scientific decision

Both requested policies are implemented. The supervisor still needs to confirm
which profile should be the primary dissertation result. No code change is
required after that decision; the selected profile is supplied at runtime.
