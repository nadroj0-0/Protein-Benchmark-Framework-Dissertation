# PFP embedding inventory and reuse planner

This package reads the nine `{bp,cc,mf}-{training,validation,test}.csv` files,
inventories an existing PFP cache in place, validates every selected NumPy
array, and emits reuse/generation/masking plans. It never generates, downloads,
extracts, copies, or links embeddings.

## What is being decided

The planner keeps four questions separate:

1. Is a file physically present and readable?
2. Is it a finite floating vector of the exact PFP dimension?
3. Does the target correspond to its source ID or complete sequence?
4. Is that correspondence scientifically compatible for this modality and
   benchmark context?

Consequently, `present-valid`, `missing`, `unreadable`, `wrong-dimension`,
`non-finite`, `sequence-mismatch`, and `provenance-unknown` are factual states;
`reuse`, `generate`, `leave-masked`, `unavailable`, and `manual-review` are
requested actions. Missing or unreadable arrays are compatible with PFP's
`MultiModalDataset._load_embedding_raw()` behavior: PFP substitutes a zero
vector and mask `0.0`.

The modality rules are intentionally asymmetric:

- **ProtT5** is sequence-only. Reuse requires the exact complete amino-acid
  sequence, checked by SHA-256. Different IDs can reuse one source only when
  their full hashes are identical. Approximate or fuzzy matching is absent.
- **Temporal text** depends on the description and time context. The exact
  published CAFA3 text bytes can be reused for the verified published artifact;
  that does not establish a generation recipe for a 2025 t0 benchmark. The
  pooling/provenance gap remains unresolved, so new-context transfer and text
  generation stay manual.
- **Structure** requires compatible ID, exact sequence, and structure
  source/version evidence. Same-sequence cross-ID reuse is not inferred.
- **PPI** requires a defensible protein-to-node mapping and the same STRING
  release (the published PFP source is STRING v12). A same-named file is not
  mapping evidence.

Aliases are accepted only from an explicit TSV. They are never fuzzy, must
carry the modality evidence described by `--help`/the tests, and cannot create
artifact ownership.

## Separate target and source contracts

Schema version 2 requires both `target_benchmark_contract` and
`source_benchmark_contract`; schema version 1 is rejected rather than silently
reinterpreted. The target contract is the leakage gate. The source contract
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

An `artifact-scoped` label is not authority. Exact published-artifact reuse is
enabled only when every configured proof passes:

- target and source equal the pinned canonical benchmark fingerprint;
- target and source fingerprints equal each other;
- the full embedding-cache catalog fingerprint, modality counts, file count,
  and logical byte count match;
- all available published archive SHA-256 values match;
- the immutable PFP reference commit and pinned workflow/source file hashes
  match; and
- the selected policy is `paper-faithful`.

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
- `reuse_manifest.tsv.gz`, `generation_manifest.tsv`,
  `manual_review.tsv.gz`, `cache_extras.tsv.gz`;
- per-modality reuse, missing, masked, generation, unavailable, and manual
  lists;
- `exact_sequence_reuse.tsv`, `manual_review_reasons.tsv`, and `errors.tsv`;
- actionable `generate_prott5.fasta` when ProtT5 generation is requested.

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
