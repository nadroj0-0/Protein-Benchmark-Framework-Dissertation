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
- Source-resolution diagnostics are persisted before strict QC. A raw ID that
  is unresolved only in the nearest retained source product may use an exact
  live match in the frozen graph; `consider` terms are never selected.
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

The same package now supports a historical-only official-target catalogue.
Released CAFA target IDs and sequences remain authoritative, while conservative
UniProt mappings are reported as mapped, unmapped or ambiguous. The historical
runner selects UniProtKB 2016_08 as the nearest public pre-freeze Swiss-Prot
snapshot and retains the former February-2017 route as an explicit legacy
comparison. These controls do not change either contemporary profile.

## Operations

The reusable runner is:

```text
scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
```

The UCL Grid Engine wrapper is:

```text
hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh
```

The wrapper stages available frozen inputs in node-local scratch and lets the
reusable entrypoint download any missing inputs there. Large TrEMBL sources are
stream-filtered to the CAFA3 taxa rather than stored wholesale. It copies the
complete run directory to durable storage and clears scratch through an
EXIT/SIGINT/SIGTERM trap. It does not run PFP training or evaluation.

The public GO archive does not contain a standalone product directory for the
2025-03-07 ontology declared by GOA 225. The build freezes labels to the last
preceding release (2025-02-06) and uses 2025-03-16 only for source-ID
normalisation. This explicit approximation is recorded in the build manifest
and must be stated in the dissertation methodology.

The exact GOA 225 rows affected by that archive gap, plus valid t1-only terms
outside the frozen graph, are written to dedicated TSVs before any strict-QC
failure. Production outputs remain incomplete until the full cluster rerun
finishes and validates all CSVs.

## Remaining scientific decision

Both requested policies are implemented. The supervisor still needs to confirm
which profile should be the primary dissertation result. No code change is
required after that decision; the selected profile is supplied at runtime.
