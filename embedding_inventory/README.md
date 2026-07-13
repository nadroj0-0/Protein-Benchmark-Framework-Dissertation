# PFP embedding inventory and reuse planner

This package inventories an existing PFP embedding cache against any directory
containing the nine benchmark CSVs:

```text
{bp,cc,mf}-{training,validation,test}.csv
```

It does not generate, copy, link, download, or extract embeddings. It validates
the benchmark and arrays, classifies scientific reuse eligibility per modality,
then writes manifests for a later materialisation or generation stage.

## Architecture

- `benchmark.py` strictly parses all nine CSVs, builds the union protein table,
  validates repeated IDs/sequences/labels, and enforces the configured split
  overlap contract.
- `config.py` requires one explicit directory per modality and validates
  provenance and action configuration. Directory fallback lists are rejected.
- `inventory.py` validates every selected `.npy` shape and finite value, applies
  exact-ID, exact-sequence, and explicit-alias routes, and keeps factual state
  separate from requested action.
- `reports.py` writes full inventories, coverage summaries, cache-extra reports,
  and modality-specific reuse/generation/unavailability lists.
- `cli.py` composes the stages. `scripts/verification/inventory_embeddings.py`
  is the repository-local entry point and does not require package installation.

The source benchmark is mandatory. It defines the protein identities and
sequences represented by the embedding cache; the target benchmark defines the
proteins being planned. This prevents a same-named file from being treated as
identity or provenance evidence by itself.

## Decision model

Array validity and scientific eligibility are independent:

| Stage | Result |
|---|---|
| File absent | `missing` |
| NumPy load fails | `unreadable` |
| Physical shape is not exactly the configured 1-D shape | `wrong-dimension` |
| Any value is NaN/Inf | `non-finite` |
| Required complete sequence SHA-256 differs | `sequence-mismatch` |
| Source identity/context is not established | `provenance-unknown` |
| Configured contexts conflict | `provenance-incompatible` |
| Valid array and all modality gates pass | `present-valid` |

Reuse routes are modality-specific:

- ProtT5: exact ID plus exact complete sequence, or a deterministic source ID
  with the identical complete sequence SHA-256. Approximate similarity is never
  used.
- Structure: exact/explicitly aliased compatible ID, exact complete sequence,
  matching configured source identity, and per-protein structure source/version
  evidence. Sequence equality alone never enables cross-ID reuse.
- PPI: explicit per-protein mapping evidence plus matching configured STRING
  release/source identity. A same-named file is never sufficient.
- Text: explicit per-protein description/temporal evidence plus matching
  configured source identity. The CAFA3 role model also distinguishes current
  train/validation from historical test and rejects a mixed-role protein.

Aliases are optional and tab-separated:

```text
protein_id	source_protein_id	modality	mapping_route	source_identity	mapping_evidence
NEW_ID	SOURCE_ID	structure	curated-uniprot-alias-2026-07	esm-if1|AlphaFold|v4	structure-source:AlphaFold;structure-version:v4
```

`modality` is `prott5`, `text`, `structure`, `ppi`, or `*`. `source_identity`
must exactly match the configured modality identity before reuse. The evidence
field should identify the description digest/cutoff, structure source/version,
or STRING node/release as appropriate. Multiple distinct mappings for one
target/modality are never resolved heuristically: they become
`provenance-unknown` / `manual-review`. No identifier normalization, isoform
collapse, or fuzzy matching occurs. Alias targets absent from the benchmark are
rejected as likely typos.

Evidence uses semicolon-separated `key:value` fields and is machine-checked:

- ProtT5: `sequence-sha256`
- text: `description-sha256` (64 hex characters) and `temporal-context`
- structure: `structure-source` and `structure-version`
- PPI: `string-id` and `string-release`

Context/source/release values must exactly equal pipe-delimited tokens in the
configured source identity;
sequence evidence must equal the target's computed hash. These are still
trusted curator assertions: the planner validates their structure and internal
agreement but does not contact UniProt, AlphaFold, or STRING to re-derive them.

## Policies

- `paper-faithful`: valid eligible published arrays are reused; missing and
  unreadable arrays remain masked, matching PFP's zero-vector/mask behavior.
  Other invalid or uncertain arrays require manual review. It requests no new
  generation solely to fill publication-cache gaps.
- `maximize-coverage`: valid eligible arrays are reused. Missing/invalid arrays
  follow the modality's configured generation/unavailability action only when
  source/target provenance is established. Unknown or incompatible provenance
  remains manual review.

Provenance uncertainty always takes precedence over automatic generation or
reuse and is sent to manual review.

## Canonical CAFA3 command

The full published model explicitly uses
`exp_text_embeddings_temporal` (current descriptions for train/validation and
UniSave descriptions with a `2016-02-17` cutoff for test). The configuration
names that directory directly; it never falls back to
`exp_text_embeddings`.

The published `.npy` files do not carry per-protein text descriptions,
AlphaFold versions, STRING node mappings, or mapping routes. The conservative
CAFA3 config therefore reports present text/structure/PPI arrays as physically
valid but `provenance-unknown` / `manual-review`. ProtT5 is automatically
reusable only after exact complete-sequence verification. Missing arrays remain
masked in `paper-faithful`, so its generation manifest is empty.

```bash
python scripts/verification/inventory_embeddings.py \
  --benchmark-dir "/path/to/cafa3_raw" \
  --source-benchmark-dir "/path/to/cafa3_raw" \
  --embedding-cache "/path/to/PFP/data/embedding_cache" \
  --config configs/embedding_inventory.cafa3_published.json \
  --policy paper-faithful \
  --output-dir results/embedding_inventory/cafa3_paper_faithful
```

Output paths inside the row-level TSVs are cache-relative. The absolute cache
root is recorded once in `embedding_summary.json`, making the rows unambiguous
without repeating a long machine-local prefix hundreds of thousands of times.

## Contemporary and homology benchmarks

The executable and CSV contract do not change. Point `--benchmark-dir` at the
new nine-CSV directory, keep `--source-benchmark-dir` pointed at the benchmark
that defines the cache, and select a provenance configuration appropriate to
the scientific comparison:

```bash
python scripts/verification/inventory_embeddings.py \
  --benchmark-dir "/path/to/contemporary_2025_2026/generated" \
  --source-benchmark-dir "/path/to/cafa3_raw" \
  --embedding-cache "/path/to/PFP/data/embedding_cache" \
  --config configs/embedding_inventory.future.example.json \
  --policy maximize-coverage \
  --output-dir results/embedding_inventory/contemporary_2025_2026
```

The conservative future template permits only exact-sequence ProtT5 reuse or
generation. Text, structure, and PPI—including physically missing arrays—remain
manual review until their target extraction/description cutoff, structure
source/version, and STRING release/mapping identities and per-protein evidence
are supplied. A homology-cluster benchmark uses the same command; cluster
membership does not alter reuse rules and never substitutes for exact sequence
identity.

## Outputs

Every run writes:

- `benchmark_proteins.tsv`
- `protein_embedding_summary.tsv`
- `embedding_inventory.tsv`
- `embedding_summary.json`
- `embedding_summary.md`
- `reuse_manifest.tsv`
- `generation_manifest.tsv`
- `manual_review.tsv`
- `cache_extras.tsv`
- all modality-specific lists requested by the framework, including
  `generate_prott5.fasta`

Coverage is reported globally, by split, by ontology, and by ontology/split,
including physical presence, numerical validity, reuse, missing, generation,
masked, unavailable, manual-review, complete-four, and at-least-one counts.

Generated reports belong under ignored `results/embedding_inventory/` or an
external supplementary-results archive. They must not be committed.

## Testing

Fast synthetic and legacy-regression suite:

```bash
cd embedding_inventory
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
```

Opt-in local golden and historical integrations:

```bash
cd embedding_inventory
PFP_RUN_REAL_INTEGRATION=1 PYTHONDONTWRITEBYTECODE=1 \
  python -m unittest tests/test_real_integration.py -v
```

Set `PFP_INVENTORY_REAL_OUTPUT_ROOT` to redirect the full ignored reports.

## Scientific limitations

- Published arrays do not contain per-file sequence digests, text strings,
  structure versions/checksums, STRING IDs, or identifier mapping routes.
  Consequently, the shipped conservative configs do not automatically reuse
  text, structure, or PPI; a separately curated evidence mapping is required.
- The checked PFP text-generation scripts save multidimensional hidden states,
  while the published cache contains physical `(768,)` arrays. A shared model
  family/name alone is therefore not extraction compatibility evidence.
- Missing published structure/PPI arrays are conservatively unavailable unless
  a future config and explicit mapping/source manifest establishes generation
  or extraction feasibility.
- Internal checks prove CSV consistency, array shape/finiteness, exact set/hash
  relationships, and declared provenance gates. They do not independently prove
  biological correctness, upstream database truth, or downstream model quality.
