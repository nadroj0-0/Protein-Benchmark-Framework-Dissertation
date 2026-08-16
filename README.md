# Protein Benchmark Framework

Scheduler-neutral source code used to build, validate, train, and analyse the
benchmarks reported in the dissertation. Generated databases, embeddings,
checkpoints, results, figures, and dissertation prose are intentionally not
tracked.

The submitted workflows are:

1. CAFA3 historical validation and regenerated-embedding reproduction.
2. Contemporary global no-knowledge (global-NK) benchmarking.
3. Contemporary ontology-specific no-knowledge plus limited-knowledge
   (NK+LK) benchmarking.
4. Framework-generated UniRef50 homology benchmarks at 30, 25, 20, 15, 10,
   and 5 percent identity.

Exact policies, external revisions, workflow order, and analysis parameters
are indexed in `configs/dissertation_submission_workflows.json`.

## Layout

```text
benchmark_builders/       benchmark construction packages
benchmark_forensics/      benchmark overlap and coverage checks
benchmark_reuse_planner/  exact reuse/regenerate planning
embedding_inventory/      provenance-aware cache inventory
configs/                  submitted workflow and validation policies
scripts/data_acquisition/ frozen input acquisition
scripts/embeddings/       embedding generation and cache assembly
scripts/model_execution/  validated PFP preparation, training, evaluation
scripts/diagnostics/      reported analyses
scripts/validation/       CAFA3 validation routes
```

## Environment

The external PFP model and embedding stack uses Python 3.9.23 and the accepted
package versions recorded in `requirements.txt`. Install it into an isolated
environment appropriate for the host's CPU/CUDA platform. PyTorch Geometric
binary wheels must match the installed PyTorch and CUDA build. The standalone
homology builder uses Python 3.10 or newer and pins its smaller dependency set
inside its own package.

PFP is an external dependency and is strictly pinned to:

```text
1e04fd6d6d3c40458fd41ec1a881ed6e24de768e
```

The CAFA Assessment Tool is pinned to:

```text
d72f0a5abb66d3224bd808e2015b55f1c9d18340
```

The submitted homology workflow validates MMseqs2 release `18-8cc5c` and full
commit `8cc5ce367b5638c4306c2d7cfc652dd099a4643f`. On compatible Linux AVX2 hosts,
the input populator extracts and catalogues the authenticated executable. An
authenticated catalogue entry takes precedence; without one, explicit
`MMSEQS_BIN` remains available for another compatible installation before the
launcher falls back to `PATH`. Every selected binary must report the same full
identity. The dissertation launcher does not permit overriding that expected
identity.

Framework Git metadata is optional report-only provenance. It is never an
execution gate or part of scientific identity, so a clean clone, edited clone,
or source archive can run the framework.

## Inputs

Create a persistent input store and catalogue:

```bash
bash scripts/data_acquisition/populate_san_frozen_inputs.sh \
  --root /absolute/path/to/pfp_inputs \
  --profile all
```

The historical filename is retained for compatibility; the script has no SAN
dependency. Input identities, release guards, and known checksums are declared
in `scripts/data_acquisition/san_frozen_inputs.tsv`. Workflows resolve explicit
paths first and then `ARTIFACT_CATALOG`. Use
`configs/artifact_paths.example.tsv` as the machine-local catalogue template.

## Benchmark Construction

The contemporary global-NK route defaults to the submitted `supervisor`
profile:

```bash
DB_ROOT=/absolute/path/to/pfp_inputs \
ARTIFACT_CATALOG=/absolute/path/to/artifact_paths.tsv \
bash scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh
```

For NK+LK, set `PROFILE=supervisor-nk-lk`. See
`benchmark_builders/contemporary_cafa/README.md` for the temporal contract.

The homology launcher defaults to the submitted UniRef50, sensitivity-4,
Swiss-Prot-plus-TrEMBL profile and builds all six thresholds:

```bash
ARTIFACT_CATALOG=/absolute/path/to/artifact_paths.tsv \
OUTPUT_ROOT=/absolute/path/to/benchmark_outputs \
TEMP_DIR=/absolute/path/to/work \
bash scripts/benchmark_generation/run_homology_cluster_benchmark.sh
```

It requires authenticated common-preprocessing and cluster-cache roots from
the input catalogue. See `benchmark_builders/homology_cluster/README.md`.

## Embeddings

CAFA3 generation is embedding-only and finishes by comparing regenerated and
published arrays:

```bash
bash scripts/reproduction/run_cafa3_full_from_scratch_reproduction.sh --help
```

Contemporary global-NK uses the paper-faithful inventory/reuse plan and the
explicit temporal-text route. The accepted cache applies the `2025-03-08`
cutoff to every target through the retained full-text replacement and archive
finalization path.

NK+LK and homology use pair-resolved source ledgers. Generate sequence, text,
structure, and PPI runs separately, then assemble them with authenticated
completion markers. NK+LK uses
`configs/contemporary_nk_lk_embedding_generation.json`; homology uses
`configs/homology_embedding_generation.json`.

See `scripts/embeddings/README.md` for the exact direct commands.

## Model Execution

All submitted training and evaluation uses the generic runner:

```bash
bash scripts/model_execution/run_pfp_benchmark.sh --help
```

It validates benchmark, ontology, embedding, and external PFP contracts;
trains BPO, CCO, and MFO separately with upstream `train.py --single`; and
strictly evaluates with `norm=cafa`, `prop=max`, and `no_orphans=false`.
Prediction capture supports `--evaluation-split valid` and
`--evaluation-split test`.

The fresh CAFA3 embedding driver does not train or evaluate models. Hand its
validated cache to this runner with `configs/pfp_benchmark_run.cafa3.json`.
Global-NK, NK+LK, and homology use their corresponding
`configs/pfp_benchmark_run.*.json` files.

## Validation

Before expensive execution:

```bash
python scripts/verification/verify_splits.py --data-dir /path/to/data --strict
python scripts/verification/verify_embeddings.py \
  --data-dir /path/to/data \
  --config configs/cafa3.json \
  --strict
```

CAFA3's released-pickle, official-artifact, and historical reconstruction
checks remain separate under `scripts/validation/`. Reported root exclusion,
specificity, knowledge-cohort, modality, relationship-policy, confidence, and
homology progression analyses are under `scripts/diagnostics/`.

Outputs should always be written outside the source checkout. Publication
steps use staged writes, hashes, manifests, and completion markers; existing
completed outputs are not silently overwritten.
