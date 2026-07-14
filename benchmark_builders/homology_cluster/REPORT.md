# Implementation report

## Implemented

- A separate `homology_cluster_benchmark` src-layout package; no temporal/PFP/imported private code
  was edited.
- Frozen local-or-URL inputs with byte size, resolved path, URL, release metadata, and embedded
  GOA/OBO checks. Hashes may be omitted only for preview/fixture use; production requires all five
  expected SHA-256 values and rechecks the inputs before publication.
- Streaming plain/gzip FASTA/GAF/idmapping parsing, explicit UniProt DAT secondary-alias handling,
  and disk-backed SQLite indexes for UniRef90/MMseqs memberships.
- Exact six-threshold MMseqs2 command construction with locked 0.8 longest-sequence coverage,
  literal identity, greedy set cover, reassign, `-s 7.5`, `-e 1e-4`, and production exact-
  version/executable identity enforcement before expensive parsing.
- Complete MMseqs TSV validation and the UniProt→UniRef90→MMseqs→split chain.
- Annotation-driven cluster retention with all retained UniRef members preserved.
- Deterministic `cluster-count-random` and bounded multi-candidate/swap `sequence-balanced`
  whole-cluster algorithms, including both split stages, ratio warnings, and production tolerance.
- Annotated-only label construction with own-label propagation, all three roots retained, literal
  development/test → complete-development `min_count` → 90:10 ordering, independent artifact
  recount, and test exclusion.
- Nine PFP CSVs, five DeepGOPlus-shaped pickles, all requested manifests/summaries/reports, hashes,
  atomic publication, post-publication verification, and completion marker last.
- Repeated-root validation and scoped aggregation of six independently scheduled outputs, driven
  by hashed child metadata and without copying their large payloads.
- A real shell entrypoint and a short UCL Grid Engine wrapper using an existing `mmfp` environment.
- Complete cross-dimensional annotation decision counts, bounded deterministic rejected-row
  samples, deterministic compressed accepted/retained/canonical membership audits, and no duplicate
  full MMseq publication.
- Reviewed frozen-input manifest binding, a recomputable hashed scientific-fingerprint payload,
  scope/eligibility correlation, hashed publication metadata, and exact marker agreement.
- Low-compression streamed GOA accepted scratch records, disk-backed UniRef/MMseq indexes and combined
  exact-sequence joins, cached GO ancestor closures, bounded working-map cleanup, and periodic
  progress/stage logs for long runs.
- Tiny fixtures/tests covering parsers, evidence/NOT, alt/obsolete/unknown terms, mappings,
  no label transfer, fully rejected/missing-sequence attrition, retention, two adversarial split
  vectors, determinism, giant-cluster ratio failure, roots in all nine CSVs, 49/50 plus 50 test-only
  term boundaries, compressed/bounded audits, rehashed provenance/scope forgery, six-root
  aggregation mismatch, atomic failure, both wrapper signals, exact and pre-copy capacity,
  exact Conda-path activation, exact MMseq versions, immutable-PFP ingestion, output coexistence,
  and guarded all-member behavior.

## Protected files inspected but not reused unsafely

- `benchmark_builders/contemporary_cafa/`: GOA/OBO/DeepGOPlus semantics and reporting contracts.
- `scripts/benchmark_generation/run_contemporary_temporal_benchmark.sh`: shell acquisition and
  command-array conventions.
- `hpc_jobs/active/hpc_contemporary_temporal_benchmark.sh`: supported Grid Engine/trap/copy patterns.
- immutable PFP `scripts/prepare_cafa3_data.py` and `mmfp/dataset.py`: raw/processed data contracts.
- `embedding_inventory/` and `configs/embedding_inventory.homology.example.json`: strict downstream
  inventory/provenance workflow.
- existing validation/reproduction scripts: regression boundaries and historical pickle contracts.

The temporal exporter was not imported because it performs protein-level validation splitting,
uses a different temporal term-population policy, and silently removes
cross-split exact sequences after export. Those behaviors conflict with this benchmark contract.
No protected private helper is imported. Validated schemas, scientific semantics, failure
boundaries, and shell/HPC conventions were adapted into the isolated package instead.

## Not implemented or run

- No full UniRef90/GOA/idmapping/UniProt data was downloaded or processed.
- No real MMseqs2 clustering was run locally; the executable is absent and the integration test
  skips explicitly.
- No embeddings were generated or copied.
- PFP was not trained or evaluated.
- No scheduler job was submitted.
- No safe resume/overwrite mode was added; publication refuses an existing final directory.
- No `all-cluster-members` PFP label policy was invented.
- No authoritative deleted-accession source was added, so an exact accession absent from
  `idmapping_selected` is explicitly `obsolete-status-unknown-absent-from-idmapping`, not asserted
  obsolete.
- No commit or push is part of this implementation task.

## Local validation results

The final package checks were:

```bash
cd benchmark_builders/homology_cluster
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -q -p no:cacheprovider
PYTHONPYCACHEPREFIX=/tmp/homology-pyc python3 -m compileall -q src tests
python3 -m ruff check --no-cache .

cd ../..
bash -n scripts/benchmark_generation/run_homology_cluster_benchmark.sh \
  hpc_jobs/active/hpc_homology_cluster_benchmark.sh
```

Results: `unittest` ran 66 cases: 65 passed and one real-MMseqs2 smoke skipped explicitly because
MMseqs2 is not installed. Pytest reported 65 passed and the same one skip. Compilation, Ruff, and
both shell syntax checks passed. ShellCheck was unavailable and therefore was not run.

The real shell entrypoint then published the tiny synthetic 30% assignment fixture. Its saved
validation report passed 18 checks with five expected risk warnings; the completion marker records
44 manifested payload files (46 files including manifest and marker), `benchmark_scope: fixture-only`, and
`production_eligible: false`. A second `validate` invocation passed.

Immutable PFP's unmodified `scripts/prepare_cafa3_data.py` consumed all nine CSVs into BPO/CCO/MFO
processed artifacts. The framework's unmodified `verify_splits.py --strict` then passed all three
splits for all three ontologies.

Protected fast regressions also passed: contemporary temporal builder 24/24; embedding inventory
54 passed with one explicit opt-in real-cache integration skip. No protected test wrote tracked
repository content.

## Risks before the first HPC run

1. Pin and verify exact frozen UniRef90, idmapping, UniProt sequence, GOA 234, and OBO files and
   expected hashes. Current-release endpoints are mutable.
2. Install or identify a modern MMseqs2 executable on compute nodes and record its exact version.
3. Run a small representative 30% smoke build to measure input DB, index, alignment/temp, SQLite,
   retained-manifest, copy-back, memory, and walltime costs.
4. Revisit provisional 64G memory / 200G scratch / 72h directives using measured evidence. The
   directives are known-valid syntax, not capacity claims.
5. Review the locked `-e 1e-4` scientifically after measured sensitivity analysis at 5–15%; any
   change requires a new reviewed benchmark contract rather than a runtime override.
6. Preserve the root-retaining PFP/TEMPROT compatibility lock unless formally revised.
7. Decide whether a frozen authoritative deleted-accession source should be added for positive
   obsolete classification.
8. Run existing embedding inventory on each published threshold before generating anything.
9. Repeat immutable PFP preparation on a representative real-MMseq output before full training.
10. Measure peak qualifying-population/frame memory; final supervised pandas frames must remain
    memory-resident to write the required historical pickle contract.

## Expected bottlenecks

- Frozen UniRef90 FASTA staging and MMseqs database creation.
- High-sensitivity prefilter/alignment at 5–15% identity.
- MMseqs temporary databases and cluster reassignment.
- Streaming the multi-gigabyte 22-column idmapping file.
- Hashing very large inputs and writing the single canonical compressed MMseq membership plus the
  distinct retained-member audit.
- Retained-member/output copy-back and filesystem inode/throughput pressure.
- Memory for qualifying accession/sequence/mapping state and final supervised dataframes. Large
  raw sources are streamed or disk-indexed, but full-production peak memory has not been measured.

## Scientific validation boundary

The implementation validates declared inputs/hashes/metadata, parser/filter behavior, mapping-chain
status, MMseq membership completeness, deterministic whole-cluster assignment, global ID and exact
sequence separation, protein-owned labels and the development-defined term universe, PFP/pickle
schemas, publication hashes, and marker
agreement. It does not prove biological optimality, exhaustive remote-homology edge discovery,
external source correctness beyond recorded evidence, full-scale runtime sufficiency, or successful
PFP training/evaluation.

Failure diagnostics under `_failed_runs/` have no completion marker and are not successful benchmark
publications. Fixture outputs are likewise explicitly non-production.
