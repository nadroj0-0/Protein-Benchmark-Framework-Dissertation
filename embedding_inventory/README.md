# PFP embedding inventory and reuse planner

This package reads the nine `{bp,cc,mf}-{training,validation,test}.csv` files,
inventories an existing PFP cache in place, validates every selected NumPy
array, and emits a binary reuse/regeneration plan. It never generates, downloads,
extracts, copies, or links embeddings.

## What is being decided

The planner keeps four questions separate:

1. Is a file physically present and readable?
2. Is it a finite floating vector of the exact PFP dimension?
3. Does the target correspond to its source ID or complete sequence?
4. Is that correspondence scientifically compatible for this modality and
   benchmark context?

Consequently, `present-valid`, `missing`, `unreadable`, `wrong-dimension`,
`non-finite`, `sequence-mismatch`, and `provenance-unknown` are factual states.
The action is always exactly `reuse` or `regenerate`. Missing, unreadable,
ambiguous, incompatible, and unproven arrays all require regeneration. This is
compatible with PFP's
`MultiModalDataset._load_embedding_raw()` behavior: PFP substitutes a zero
vector and mask `0.0`.

The modality rules are intentionally asymmetric:

- **ProtT5** is sequence-only. Reuse requires the exact complete amino-acid
  sequence, checked by SHA-256. Different IDs can reuse one source only when
  their full hashes are identical. Approximate or fuzzy matching is absent.
- **Temporal text** depends on the description and time context. The exact
  published CAFA3 text bytes can be reused for the verified published artifact;
  that does not establish a generation recipe for a 2025 t0 benchmark. The
  pooling/provenance gap remains unresolved, so new-context text embeddings
  receive action `regenerate`.
- **Structure** requires compatible ID, exact sequence, and structure
  source/version evidence. Same-sequence cross-ID reuse is not inferred.
- **PPI** requires the same STRING release and extractor policy (the published
  PFP source is STRING v12). In the contemporary config, a valid authenticated
  published array can be reused for the same source and target UniProt
  accession. Cross-ID PPI reuse is unsupported until an authenticated STRING
  mapping artifact is integrated.

Aliases are accepted only from an explicit TSV and are never fuzzy. ProtT5 may
use an alias only after exact full-sequence equality. Text, structure, and PPI
aliases are retained as diagnostics but cannot authorize reuse without a future
authenticated external mapping/input artifact.

## Separate target and source contracts

Schema version 3 requires both `target_benchmark_contract` and
`source_benchmark_contract` and defines the binary reuse/regenerate action
contract; earlier schemas are rejected rather than silently reinterpreted. The
target contract is the leakage gate. The source contract
only validates the identities, sequences, and memberships represented by a
pre-existing cache, so it may deliberately allow a split policy that would be
unacceptable for a new target.

Available overlap policies include:

- `global-evaluation-disjoint`: the union of training+validation is disjoint
  from test across all ontologies;
- `global-disjoint`: training, validation, and test are pairwise disjoint
  globally;
- `per-ontology-disjoint`: legacy per-ontology checking;
- `allow`: no split-overlap assertion.

`configs/embedding_inventory.contemporary.json` rejects a protein ID appearing
in train/validation and test even across ontologies. The homology template also
rejects an identical sequence under a different ID across that boundary.
Nine CSVs contain no homology-cluster assignments, so the planner does **not**
claim to validate cluster-level separation. That requires an external cluster
assignment manifest and a future explicit checker.

## Deterministic artifact proof

An `artifact-scoped` label is not authority. Cache authentication first requires:

- the source benchmark equals the pinned canonical benchmark fingerprint;
- the full embedding-cache catalog fingerprint, modality counts, file count,
  and logical byte count match;
- all available published archive SHA-256 values match;
- the immutable PFP reference commit and pinned workflow/source file hashes
  match.

Authentication establishes what the cache is; it does not by itself authorize
reuse for a changed target. An `artifact-scoped` record receives action `reuse`
only when the target also has the pinned canonical fingerprint, the policy is
`paper-faithful`, and the valid array is selected by direct protein ID.

The benchmark fingerprint hashes canonical JSON lines sorted by protein ID.
Each line includes the exact protein ID, full sequence SHA-256, sequence length,
and sorted `(ontology, split)` memberships. It is stable under irrelevant CSV
row order, while any ID, sequence, split, or ontology membership change changes
the identity. GO label vectors are validated as binary while streaming but are
not embedding-identity inputs.

The cache catalog hashes every configured `.npy` byte stream using globally
sorted lines:

```text
cache/relative/path<TAB>logical-byte-size<TAB>file-sha256<LF>
```

Directory names, artifact names, aliases, and counts alone cannot pass this
gate. Cross-benchmark ProtT5 can still transfer by exact sequence; text,
structure, and PPI do not inherit canonical artifact scope.

Every CLI run records `run_provenance.json` and a concise Markdown version with
UTC time, exact command, repository commit/dirty state, Python/NumPy/package
versions, config and alias hashes, all 18 target/source CSV roles and hashes
(shared files are hashed once), semantic benchmark fingerprints, cache catalog,
archive/reference checks, policy, report level, and runtime paths/options.
Execution outside Git is recorded rather than failing.

## Reporting and storage

`--report-level compact` is the default. It writes gzip-compressed compact
protein/action tables without repeating sequences, plus:

- `embedding_summary.{json,md}` and `run_provenance.{json,md}`;
- `benchmark_proteins.tsv.gz`, `protein_embedding_summary.tsv.gz`;
- full-field, gzip row-level `embedding_inventory.tsv.gz`;
- primary `reuse.tsv` and `regenerate.tsv` action manifests;
- `reuse/{prott5,text,structure,ppi}.txt` reusable ID lists;
- `regenerate/{prott5,text,structure,ppi}.txt` generation ID lists and
  `regenerate/prott5.fasta`;
- `regenerate_reasons.tsv` with compact factual/reason counts;
- `cache_extras.tsv.gz`, `exact_sequence_reuse.tsv`, and `errors.tsv` as
  supporting diagnostics;

`--report-level full` adds `benchmark_proteins_full.tsv.gz`; it is the only general table that
repeats complete sequences. Gzip output is deterministic. Output must be a new
directory and cannot be nested inside benchmark inputs or the cache.
`output_manifest.json` hashes every payload and `RUN_COMPLETE.json` is written
last. Reports are built in a sibling staging directory and atomically renamed,
so a failed run cannot look complete. Use a temporary directory for tests and
a deliberate archive for persistent integrations. The CLI rejects output
beneath the immutable artifact/PFP root. Reports are ignored by Git.

## Commands

Run from the repository root. `--artifact-root` is optional when the cache is
`PFP_reference_clone/data/embedding_cache`; otherwise point it at the directory
containing the three archives and pinned reference clone.

Canonical CAFA3, exact published artifact:

```bash
python3 scripts/verification/inventory_embeddings.py \
  --benchmark-dir "/path/to/cafa3_raw" \
  --source-benchmark-dir "/path/to/cafa3_raw" \
  --embedding-cache "/path/to/PFP_reference_clone/data/embedding_cache" \
  --artifact-root "/path/to/PFP_reference_clone" \
  --config configs/embedding_inventory.cafa3_published.json \
  --policy paper-faithful \
  --report-level compact \
  --output-dir "/path/to/results/cafa3_paper_faithful"
```

Historical/reconstructed target against the canonical source cache:

```bash
python3 scripts/verification/inventory_embeddings.py \
  --benchmark-dir "/path/to/historical/generated" \
  --source-benchmark-dir "/path/to/cafa3_raw" \
  --embedding-cache "/path/to/PFP_reference_clone/data/embedding_cache" \
  --artifact-root "/path/to/PFP_reference_clone" \
  --config configs/embedding_inventory.cafa3_published.json \
  --policy maximize-coverage \
  --report-level compact \
  --output-dir "/path/to/results/historical_plan"
```

Future contemporary temporal benchmark:

```bash
python3 scripts/verification/inventory_embeddings.py \
  --benchmark-dir "/path/to/contemporary_2025_2026/generated" \
  --source-benchmark-dir "/path/to/cafa3_raw" \
  --embedding-cache "/path/to/PFP_reference_clone/data/embedding_cache" \
  --config configs/embedding_inventory.contemporary.json \
  --policy maximize-coverage \
  --report-level compact \
  --output-dir "/path/to/results/contemporary_plan"
```

For a real contemporary benchmark on UCL Grid Engine, the framework also
provides a staged integration workflow. It downloads and authenticates the
three published embedding archives in node-local scratch, downloads and
authenticates the canonical CAFA3 source CSVs, runs this inventory, copies only
compact reports/lists home, and clears scratch:

```bash
qsub -v BENCHMARK_DIR=/path/to/contemporary/run/outputs \
  hpc_jobs/active/hpc_contemporary_embedding_inventory.sh
```

`BENCHMARK_DIR` is never hard-coded. Individual roles can instead be supplied
with `BP_TRAINING_CSV`, `BP_VALIDATION_CSV`, `BP_TEST_CSV`, and the equivalent
`CC_*` and `MF_*` variables. A directory may be combined with one or more
per-file overrides. Optional local caches are `SOURCE_BENCHMARK_DIR` and
`PUBLISHED_EMBEDDING_ARCHIVE_DIR`; by default both source CSVs and embeddings
are downloaded during the job.

The reusable non-scheduler implementation is:

```bash
bash scripts/verification/run_contemporary_embedding_inventory.sh \
  --benchmark-dir /path/to/contemporary/run/outputs \
  --work-dir /path/to/disposable/work \
  --output-dir /path/to/new/result
```

The result root contains `job_summary.md`, acquisition provenance,
`physical_coverage/{valid,not_valid,missing}_{modality}.txt` convenience lists, and an
`inventory/` directory with the normal manifests and completion marker.
`regenerate/prott5.fasta` is immediately actionable because ProtT5 transfer
is approved only for an exact complete-sequence match. The contemporary config
regenerates text and structure unless the exact original description or
AlphaFold structure input is positively proven. It permits direct-ID PPI reuse
only for the fixed STRING v12 published extractor identity. Missing, invalid,
ambiguous, unavailable, and unknown-provenance cases all receive action
`regenerate`; their detailed reasons remain auditable.

Future homology benchmark uses the same executable and nine CSV interface:

```bash
python3 scripts/verification/inventory_embeddings.py \
  --benchmark-dir "/path/to/homology/generated" \
  --source-benchmark-dir "/path/to/cafa3_raw" \
  --embedding-cache "/path/to/PFP_reference_clone/data/embedding_cache" \
  --config configs/embedding_inventory.homology.example.json \
  --policy maximize-coverage \
  --report-level compact \
  --output-dir "/path/to/results/homology_plan"
```

The command is unchanged for future nine-CSV benchmarks; the target contract
and modality evidence config express the benchmark-specific scientific rules.
The homology command checks global ID and exact-sequence separation, not cluster
membership.

## Tests

```bash
cd embedding_inventory
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

The opt-in real test defaults to a temporary output directory. Supply an
explicit persistent archive only when reports should be retained:

```bash
PFP_RUN_REAL_INTEGRATION=1 \
PFP_INVENTORY_REAL_OUTPUT_ROOT="/path/to/supplementary/results/run" \
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest tests/test_real_integration.py -v
```

Internal validation proves schema consistency, exact set/hash relations,
array shape/finiteness, and the recorded artifact gate. It does not establish
external biological correctness, cluster separation without assignments, or a
missing temporal-text pooling procedure.
