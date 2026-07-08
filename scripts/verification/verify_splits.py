#!/usr/bin/env python3
"""
verify_splits.py — contract gate for regenerated data splits.

Validates the OUTPUT of the data-prep stage (the protein-ID "spine" every
downstream stage keys off) using only INTERNAL CONSISTENCY — no reference to
any specific benchmark's values. These checks hold for ANY CAFA-style dataset
plugged into this pipeline, so the gate survives swapping the benchmark.

Per aspect/split it verifies the data contract:
  - names.npy exists, non-empty, no duplicate IDs
  - labels.npz row count == number of proteins
  - labels.npz col count == len(go_terms.json)
  - every protein in names has a non-empty sequence in sequences.json
  - sequences.json keys == names (no missing, no extra)
Per aspect, across splits:
  - no protein ID appears in more than one split (train/valid/test disjoint)
  - go_terms.json present and non-empty

Usage:
    python scripts/verification/verify_splits.py --data-dir data
    python scripts/verification/verify_splits.py --data-dir data --strict   # exit 1 on any FAIL
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import scipy.sparse as ssp

ASPECTS = ["BPO", "CCO", "MFO"]
SPLITS = ["train", "valid", "test"]

def load_names(p):
    return [str(x) for x in np.load(p, allow_pickle=True)]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--aspects", nargs="+", default=ASPECTS,
                    help="override for non-CAFA benchmarks")
    ap.add_argument("--splits", nargs="+", default=SPLITS)
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    data = Path(args.data_dir)
    fails, warns = [], []
    def fail(m): fails.append(m); print(f"  [FAIL] {m}")
    def ok(m):   print(f"  [ OK ] {m}")

    print(f"==> Verifying split contract under {data}\n")

    for a in args.aspects:
        print(f"--- {a} ---")
        # go_terms
        gt_path = data / f"{a}_go_terms.json"
        if not gt_path.exists():
            fail(f"{a}_go_terms.json missing"); n_terms = None
        else:
            go_terms = json.load(open(gt_path))
            if not go_terms:
                fail(f"{a}_go_terms.json empty"); n_terms = None
            else:
                n_terms = len(go_terms)
                if len(set(go_terms)) != n_terms:
                    fail(f"{a}_go_terms.json has duplicate GO terms")
                else:
                    ok(f"go_terms: {n_terms} unique")

        seen_ids = {}  # pid -> split, for cross-split disjointness
        for s in args.splits:
            names_p = data / f"{a}_{s}_names.npy"
            labels_p = data / f"{a}_{s}_labels.npz"
            seqs_p = data / f"{a}_{s}_sequences.json"

            if not names_p.exists():
                fail(f"{a}_{s}_names.npy missing"); continue
            names = load_names(names_p)
            n = len(names)
            if n == 0:
                fail(f"{a}_{s} has 0 proteins"); continue

            # duplicate IDs within split
            if len(set(names)) != n:
                fail(f"{a}_{s} has {n - len(set(names))} duplicate protein IDs")

            # labels shape contract
            if not labels_p.exists():
                fail(f"{a}_{s}_labels.npz missing")
            else:
                M = ssp.load_npz(labels_p)
                if M.shape[0] != n:
                    fail(f"{a}_{s} labels rows={M.shape[0]} != proteins={n}")
                if n_terms is not None and M.shape[1] != n_terms:
                    fail(f"{a}_{s} labels cols={M.shape[1]} != go_terms={n_terms}")

            # sequences contract
            if not seqs_p.exists():
                fail(f"{a}_{s}_sequences.json missing")
            else:
                seqs = json.load(open(seqs_p))
                missing = [p for p in names if p not in seqs]
                empty = [p for p in names if seqs.get(p, "") == ""]
                extra = [k for k in seqs if k not in set(names)]
                if missing: fail(f"{a}_{s}: {len(missing)} proteins have no sequence entry")
                if empty:   fail(f"{a}_{s}: {len(empty)} proteins have empty sequence")
                if extra:   warns.append(f"{a}_{s}: {len(extra)} extra seqs not in names")

            # cross-split disjointness
            for pid in names:
                if pid in seen_ids and seen_ids[pid] != s:
                    fail(f"{a}: protein {pid} in BOTH {seen_ids[pid]} and {s}")
                    break
                seen_ids[pid] = s

            if not any(f"{a}_{s}" in m for m in fails):
                ok(f"{s}: {n} proteins, contract holds")
        print()

    print("="*60)
    for w in warns: print(f"  [WARN] {w}")
    if fails:
        print(f"\nRESULT: FAIL — {len(fails)} contract violation(s). "
              f"Do NOT proceed to embedding generation.")
        sys.exit(1 if args.strict else 0)
    print("\nRESULT: PASS — split contract holds. Safe to generate embeddings.")
    sys.exit(0)

if __name__ == "__main__":
    main()
