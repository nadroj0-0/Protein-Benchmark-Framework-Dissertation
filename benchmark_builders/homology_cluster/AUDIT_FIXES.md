# Homology-cluster builder audit fixes

Status updated 2026-07-14 after local implementation and fixture validation. “Resolved” below means
the software contract and automated miniature tests pass; it does not mean a full UniRef90 or real
Grid Engine production run has validated the biological population or resource estimates.

## 1. UniProtKB population and provenance — resolved at interface/fixture level

The undifferentiated `uniprot_sequences` input was replaced with explicit
`sprot-only`, `trembl-only`, and `sprot-and-trembl` scopes plus separate Swiss-Prot and TrEMBL
declarations. Production has no default. UniRef90 remains the clustering scaffold; source scope
controls accession eligibility and supervised rows.

Frozen manifest schema v2 requires exactly five Swiss-Prot roles, five TrEMBL roles, or six
combined roles. It binds role, source population, release, URL, local filename, expected/observed
hash, size, acquisition action, embedded metadata, and origin review. Production sequence sources
must be DAT. FASTA is fixture/diagnostic-only and its missing secondary-accession information is
reported. DAT `ID` content is checked too: Swiss-Prot requires `Reviewed` and TrEMBL requires
`Unreviewed`, so a renamed or misdeclared product fails its source role.

Combined-mode ingestion is fixed-order and source-tagged. A disk-backed collision audit fails
duplicate primaries and conflicting sequences in strict production mode. Identical-sequence
secondary aliases attached to multiple primaries are counted and excluded if qualifying, avoiding
both arbitrary canonicalization and whole-run failure. Source-specific reports separate
primary/secondary selected-UniProt resolution, present UniRef90 mapping, and MMseqs2 assignment,
with explicit zero rows for selected sources. Tests cover all three scopes, required/irrelevant
sources, 5/5/6 manifests, role/scope/release/hash mismatches, determinism, primary/secondary
mappings, collisions, and out-of-scope retention denial.

Remaining boundary: no full official Swiss-Prot/TrEMBL bytes were processed. Daniel still needs to
choose the production scope.

## 2. Severe attrition gate — resolved at policy/fixture level

Production now requires a reviewed structured policy bound to source scope, releases, full commit,
and frozen-manifest hash. Eleven registered metrics cover selected-source and UniRef90 mapping,
qualifying annotation retention, retained-cluster members, raw/evaluable proteins, propagated and
ontology-specific evaluability, and both split deviations. Reports preserve definitions,
numerators, denominators, ratios, limits, and outcomes.

Exceeding a limit prevents production publication unless an exact structured override binds the
observed failures and run evidence. There is no boolean bypass, and the full-array launcher rejects
one array-wide override. Diagnostic pilot output always reports both non-production eligibility and
`production_authorized: false`. Policy validation occurs before expensive scientific work, and the
manifest/policy contract fails before large-input hashing; the policy hash is rechecked before
publication. Boundary, just-failing, wrong-binding, placeholder,
malformed-policy, malformed-override, and valid-override tests pass.

Remaining boundary: no biological limit values have been approved. They must come from reviewed
pilot evidence rather than software defaults.

## 3. Always clear HPC scratch — resolved for tested normal exit paths

The worker atomically creates a job/task-owned scratch directory and persistent claim. Success,
builder failure, both validation failures, INT, TERM, pre-copy failure, and persistent copy failure
attempt marker-free diagnostic publication and then safe cleanup. Copy failure no longer preserves
scratch. Cleanup failures propagate non-zero status and are not masked.

Empty, root, relative, symlinked, pre-existing, outside-base, or ownership-mismatched paths are
refused. A colliding invocation cannot overwrite an ownership marker, delete another task's
scratch, or move a prior successful publication to `.failed`. Tests cover successful/failing
cleanup, both signals, copy/validation failures, unsafe markers/paths, symlinks, pre-existing work,
persistent collisions, and idempotent absence handling.

Remaining boundary: SIGKILL and node loss cannot execute shell traps; no such guarantee is made.

## 4. Exact framework revision — resolved

Production and diagnostic pilots require exactly 40 lowercase hexadecimal characters. The worker
checks out that commit detached, compares exact `git rev-parse HEAD`, rejects a symbolic branch or
dirty tree, and records the requested/observed commit in task context, scientific fingerprint,
publication metadata, and completion marker. Tests cover missing, branch, tag, abbreviated,
uppercase/malformed, checkout failure, HEAD mismatch, dirty checkout, and an accepted full detached
SHA before later work.

## 5. Grid Engine parallel environment and array workflow — implemented; capacity pending pilot

The worker requests `-pe smp 2`. `NSLOTS` is authoritative, production requires two, and MMseqs2
threads must match. Task IDs map exactly `1→30`, `2→25`, `3→20`, `4→15`, `5→10`, `6→5`.
The pilot launcher uses `-t 1`; the full launcher uses `-t 1-6`. Both export a whitelist and shared
local checksum-pinned inputs with `NO_DOWNLOADS=1`. Six concurrently runnable tasks may request up
to 12 slots, but Grid Engine controls scheduling.

The full launcher requires a validated 30% diagnostic publication, reviewed attrition policy,
human approval, task context, and separately sourced resource measurements. Launcher-time evidence
hashes are rechecked by the queued worker, which reruns authorization from its detached revision.
Approval binds the reviewed-policy hash and one pilot run ID across marker, task context, and
measurement evidence. Authorization reproduces every saved pilot metric and requires the reviewed
production policy to accept the pilot observations; templates cannot pass as reviewed text.
Fake-`qsub` tests prove dry run makes no scheduler call, exact task/PE arguments are used, approval
failure blocks submission, and a fake scheduler failure propagates.

Remaining boundary: the provisional `64G / 200G / 72h` resources are valid syntax, not evidence of
sufficiency. Runtime, memory, scratch peak, output size, and copy-back require a real reviewed 30%
pilot.

## 6. Publish before production — resolved by the repository commit containing this file

Production rejects dirty or non-commit-addressable framework state. The implementation must be
pushed before a cluster clone can use it. Push success is operational evidence outside this file
and must only be reported after visible Git output; this document does not assume it.

## Current validation boundary

Local tests validate parser behavior, source/manifests, mapping and collision contracts, attrition
and approval bindings, whole-cluster/sequence leakage, output files, aggregation, scheduler command
construction, revision checks, publication isolation, and scratch cleanup using miniature fixtures
and fake executables. Immutable PFP ingestion is exercised when the local checkout is available.

During implementation, no real `qsub`, `qstat`, `qdel`, SSH, full source download, full MMseqs2
clustering, production HPC build, embedding generation, or PFP training/evaluation was performed.
The independent audit below later used read-only SSH environment probes only; it did not submit or
alter cluster work. Real-data mapping coverage, biological correctness, runtime, peak memory, peak
scratch, and copy-back capacity remain unproven.

## Independent post-implementation audit — 2026-07-14

The implementation was reviewed again after commit `a0ebe2c`. The six original engineering findings
above are resolved at code and miniature-fixture level. No new defect was found in the benchmark
transformations, cluster-level split enforcement, PFP export contract, provenance binding, atomic
publication, or owned-scratch cleanup logic.

### Verification repeated independently

- `PYTHONPATH=src python3 -m unittest discover -s tests -v`: 102 tests passed; one real-MMseqs2
  smoke test was skipped because MMseqs2 is not installed locally.
- Ruff with cache disabled, Python compilation with its bytecode cache under `/tmp`, `bash -n` over
  the five homology shell entrypoints, and `git diff --check` all passed.
- The immutable-PFP fixture integration test consumed all nine generated CSVs.
- A read-only UCL cluster probe found the existing `mmfp` environment uses Python 3.11.15, so the
  package's Python `>=3.10` declaration is compatible with the actual cluster environment.
- `qconf -spl` and `qconf -sp smp` confirmed that the `smp` parallel environment exists and uses
  `$pe_slots` allocation. This validates the PE name and slot-allocation shape, but not the provisional
  memory, scratch, or walltime requests.

### Outstanding blockers and decisions

1. **The real 30% pilot cannot currently start on the cluster.** MMseqs2 was not found on `PATH`, in
   the `mmfp` environment, through the available modules, or at the checked common shared paths. The
   exact executable/version pin cannot be completed until MMseqs2 is installed or an approved shared
   binary is identified.
2. **The shared frozen input collection is incomplete.** The remaining
   `~/protein_databases` tree contains the 2026_02 Swiss-Prot DAT and small metadata only. It does not
   currently contain the frozen UniRef90 FASTA, `idmapping_selected`, GOA 234 GAF, GO OBO, or TrEMBL
   DAT required by scopes that include TrEMBL. The production launcher deliberately forbids six
   independent downloads, so these inputs need persistent shared storage and reviewed hashes first.
3. **Resource values remain provisional.** `64G / 200G / 72h` has not been measured. The pilot must
   establish whether these requests are schedulable and sufficient before they are reused for all six
   array tasks.
4. **The full-array launcher currently makes the 30% pilot mandatory.** This matches Daniel's stated
   good-practice instruction to test one array item first and is conservative, but it is stricter than
   the later preference that a pilot should be strongly recommended rather than mechanically required.
   Resolve that policy before changing the launcher; it does not block the recommended pilot-first path.
5. **Scientific choices remain explicit rather than silently defaulted.** Daniel still needs to select
   the production UniProt source scope and one of the two implemented cluster split policies. Only
   `annotated-only` supervision is implemented. `all-cluster-members` continues to fail deliberately
   until a valid label policy for unannotated proteins is approved.

### Readiness verdict

- **Fixture/software validation:** ready.
- **Production launcher preview:** ready once the pinned input collection and MMseqs2 executable exist.
- **One-item 30% diagnostic pilot:** code-ready, but operationally blocked by missing MMseqs2 and
  incomplete shared inputs.
- **Six-task production array:** not yet authorized. It requires the successful pilot, reviewed input
  manifest and source scope, exact MMseqs2 pin, evidence-derived attrition limits, measured resources,
  and human approval described above.

Therefore the implementation should be pushed for reproducibility and cluster staging, but it must
not yet be described as production-validated or submitted as the full sensitivity array.

## Scratch-first runtime submission follow-up — 2026-07-15

The earlier readiness verdict described only the pre-staged, human-approved launcher. Two final
runtime entrypoints now remove those operational blockers without changing the scientific builder:

- `hpc_homology_cluster_runtime_pilot.sh` submits task 1 / 30% and performs automatic review only.
- `hpc_homology_cluster_runtime_array.sh` submits tasks 1-6 directly; the pilot is recommended but
  is not required or consumed as authorization evidence.

Their shared driver installs exact MMseqs2 `18-8cc5c` in job-owned scratch when needed, stages
provided inputs or downloads missing official UniProt 2026_02, GOA 234, GO 2026-06-15, UniRef90,
and idmapping bytes into scratch, validates release markers and embedded metadata, creates the
frozen manifest, delegates to the existing builder, validates every publication, and copies only
outputs/reports/logs home. The six-task wrapper caps concurrency at two because array tasks cannot
share node-local scratch and therefore independently download omitted inputs.

The combined-scope full-array scratch request was set to 1200 GB rather than the old provisional 200 GB. This is
based on the builder's conservative input, parser-index, MMseqs-work, and publication estimates;
it remains unmeasured until a real pilot. Runtime tests prove success copy-back, direct array
execution with no pilot variables, and the required behavior when home copy-back fails: non-zero
exit, incomplete-home cleanup, and unconditional owned-scratch deletion.

This follow-up resolves software installation, missing persistent input, and mechanically mandatory
pilot blockers for the scratch-first workflow. It does not claim that a full production task has
run, that 1200 GB is an empirically optimized request, or that automatic non-blocking attrition
bounds replace scientific review.

The subsequent diagnostic-pilot wrapper reduces the requested scratch to 300 GB and records actual
job-owned usage every 120 seconds plus explicit stage checkpoints. Pilot mode records rather than
enforces the speculative preflight estimate and uses neutral `1x` multipliers. The resulting peak
measurement is intended to replace the provisional full-array request after review.
