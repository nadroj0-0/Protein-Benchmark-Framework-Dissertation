# CAFA3 organiser replay experiment

This branch contains an isolated experiment for one boundary only:

```text
historical UniProt + GOA + GO snapshots
                  -> CAFA3 benchmark lists and leaf-only ground truth
```

It deliberately stops before DeepGOPlus, TEMPROT, the nine PFP CSVs,
embeddings, model training, and evaluation. Those downstream stages are already
validated separately and are not part of this experiment.

## Why this is a branch

The experiment lives on `codex/cafa3-organizer-replay`, based on the current
`fixed-evidence-codes` operational framework. The ordinary framework workflows
remain the default. The experiment adds a new entrypoint and does not replace
`cafa-benchmark-builder`.

## Targeted alignment

The replay reuses the framework's existing streaming UniProt parser, GAF
filtering and normalisation, official-target mapper, and GO graph. It aligns
the organiser-facing boundary with the best recovered evidence:

| Policy | Replay default | Evidence |
|---|---|---|
| Evidence codes | EXP, IDA, IPI, IMP, IGI, IEP, TAS, IC | final CAFA3 README |
| Type 1 | t1 aspect membership with no t0 annotation in any GO aspect | CAFA2 code |
| Type 2 | t1 aspect membership with no t0 annotation in that aspect, excluding type 1 | CAFA2 code |
| Backfill | remove t1 rows assigned on or before 2017-02-13 | final CAFA3 README |
| Protein binding | remove protein-binding-only MFO records from both snapshots before classification | CAFA2 code and final CAFA3 README |
| GO graph | `is_a + part_of` | CAFA2 `pfp_ontbuild` default |
| Leaf truth | propagate, then retain only the most specific surviving terms | CAFA2 `pfp_leafannot` |

The final CAFA3 code revision and glue script were not recovered. In
particular, exact NOT handling, exact historical GOA bytes, the final
relationship selection, and any manual corrections remain explicit
reproducibility boundaries. `--relationship-policy is-a` and
`--relationship-policy all` are therefore retained as named sensitivity modes;
they do not overwrite the CAFA2-aligned default.

## Run

Install the contemporary builder in an isolated environment, or expose its
source directory through `PYTHONPATH`, then run:

```bash
PYTHONPATH=benchmark_builders/contemporary_cafa/src \
python3 -m cafa_benchmark_builder.organizer_replay \
  --uniprot-t0 /path/to/t0-uniprot.dat.gz \
  --uniprot-t1 /path/to/t1-uniprot.dat.gz \
  --goa-t0 /path/to/goa_uniprot_all.gaf.163.gz \
  --goa-t1 /path/to/goa_uniprot_all.gaf.172.gz \
  --go-obo /path/to/frozen-benchmark-go.obo \
  --go-obo-t0 /path/to/t0-go.obo \
  --go-obo-t1 /path/to/t1-go.obo \
  --official-target-root /path/to/extracted/CAFA3_targets \
  --official-benchmark-dir /path/to/extracted/benchmark20171115 \
  --output-dir /path/to/new/replay-output
```

The target root must contain `Target files/` and `Mapping files/`, exactly as
released in `CAFA3_targets.tgz`. The output path must not already exist.

## Outputs

- `benchmark20171115/lists/{bpo,cco,mfo}_all_type{1,2,x}.txt`
- `benchmark20171115/lists/xxo_all_typex.txt`
- `benchmark20171115/groundtruth/leafonly_{BPO,CCO,MFO}.txt`
- `benchmark20171115/groundtruth/leafonly_all.txt` as a downstream compatibility aggregate
- numbered target-mapping, pre/post-filter annotation, and membership-flow TSV checkpoints
- set-level comparisons with the official final CAFA3 release
- a policy/evidence manifest, input hashes, output hashes, and completion marker

Aspect and leaf comparisons are set-based, so harmless output ordering
differences cannot inflate or hide agreement. The combined `xxo` list is
compared as a multiset because cross-aspect proteins repeat legitimately. The
report includes generated-only and official-only records, precision, recall,
and Jaccard for every all-species release artefact.

## Reference sources

- CAFA3 supplementary dataset: <https://doi.org/10.6084/m9.figshare.8135393.v3>
- CAFA2 benchmark/evaluation code: <https://github.com/yuxjiang/CAFA2>
- Public CAFA benchmark comparator: <https://github.com/nguyenngochuy91/CAFA_benchmark>

The official release is the output oracle. Matching it does not prove access to
the organisers' exact private historical execution; disagreement is retained as
evidence rather than tuned away without a source-backed policy reason.

## Grid execution

Submit the isolated four-condition replay matrix from a clean checkout of this
branch:

```bash
bash hpc_jobs/submit_cafa3_organizer_replay.sh
```

The job downloads the historical inputs once and runs:

1. assigned-date endpoint with `is_a + part_of` (primary CAFA2-aligned policy);
2. assigned-date endpoint with `is_a` only;
3. assigned-date endpoint with all parsed GO relationships; and
4. snapshot-membership endpoint with `is_a + part_of`.

Every condition is published separately under one matrix run and carries its
own `REPLAY_COMPLETE.json`. The matrix root is published only after all four
conditions complete.
