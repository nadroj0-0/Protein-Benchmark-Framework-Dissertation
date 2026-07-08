#!/usr/bin/env python3
"""
generate_embeddings_fasta.py — build data/proteins.fasta from the split sequence
JSONs, keyed by the EXACT protein IDs used in the splits.

ProtT5 (extract_prott5_embeddings.py) writes one {record.id}.npy per FASTA
record; the dataset loader later reads {protein_id}.npy using IDs from
names.npy. Therefore FASTA record IDs MUST equal the split protein IDs verbatim
— any prefix/accession mismatch makes every embedding unresolvable downstream.

Deduplicates proteins that appear across multiple aspects/splits, and verifies
that a repeated protein has an identical sequence everywhere (a data-integrity
check, not an assumption).

Usage:
    python scripts/embeddings/generate_embeddings_fasta.py --data-dir data
    python scripts/embeddings/generate_embeddings_fasta.py --data-dir data --config configs/cafa3.json
"""
import argparse
import json
import sys
from pathlib import Path

DEFAULT_ASPECTS = ["BPO", "CCO", "MFO"]
DEFAULT_SPLITS = ["train", "valid", "test"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--config", default=None,
                    help="optional JSON with 'aspects'/'splits' (else CAFA3 defaults)")
    ap.add_argument("--out", default=None, help="output FASTA (default: <data-dir>/proteins.fasta)")
    ap.add_argument("--line-width", type=int, default=60, help="FASTA wrap width (0 = no wrap)")
    args = ap.parse_args()

    data = Path(args.data_dir)
    aspects, splits = DEFAULT_ASPECTS, DEFAULT_SPLITS
    if args.config:
        with open(args.config) as f:
            cfg = json.load(f)
        aspects = cfg.get("aspects", aspects)
        splits = cfg.get("splits", splits)

    out = Path(args.out) if args.out else data / "proteins.fasta"

    seqs: dict[str, str] = {}          # pid -> sequence
    conflicts: list[str] = []          # pids with differing sequences across files
    files_read = 0

    for a in aspects:
        for s in splits:
            f = data / f"{a}_{s}_sequences.json"
            if not f.exists():
                continue
            files_read += 1
            with open(f) as fh:
                d = json.load(fh)
            for pid, seq in d.items():
                if pid in seqs:
                    if seqs[pid] != seq:
                        conflicts.append(pid)
                else:
                    seqs[pid] = seq

    if files_read == 0:
        print(f"FAIL: no *_sequences.json found in {data} — run prepare_data first.")
        sys.exit(1)

    if conflicts:
        uniq = sorted(set(conflicts))
        print(f"FAIL: {len(uniq)} protein ID(s) have DIFFERENT sequences across files, e.g.:")
        for pid in uniq[:5]:
            print(f"    {pid}")
        print("  Deduplication would silently pick one. Investigate before proceeding.")
        sys.exit(1)

    empties = [pid for pid, s in seqs.items() if not s]
    if empties:
        print(f"FAIL: {len(empties)} protein(s) have empty sequences, e.g. {empties[:5]}")
        sys.exit(1)

    with open(out, "w") as fh:
        for pid, seq in seqs.items():
            fh.write(f">{pid}\n")
            if args.line_width and args.line_width > 0:
                for i in range(0, len(seq), args.line_width):
                    fh.write(seq[i:i + args.line_width] + "\n")
            else:
                fh.write(seq + "\n")

    print(f"==> Wrote {len(seqs)} unique proteins to {out} "
          f"(from {files_read} sequence files, no ID/sequence conflicts).")


if __name__ == "__main__":
    main()