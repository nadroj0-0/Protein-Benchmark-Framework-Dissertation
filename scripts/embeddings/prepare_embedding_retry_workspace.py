#!/usr/bin/env python3
"""Restrict a prepared PFP data directory to retry and control proteins."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Dict, Set

import numpy as np
from scipy import sparse


ASPECTS = ("BPO", "CCO", "MFO")
SPLITS = ("train", "valid", "test")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_ids(path: Path, modality: str) -> Set[str]:
    if not path.is_file():
        raise ValueError(f"Missing pair table: {path}")
    result: Set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["modality"] != modality:
                raise ValueError(
                    f"Pair table {path} contains {row['modality']}, expected {modality}"
                )
            result.add(row["protein_id"])
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--requested-pairs", type=Path, required=True)
    parser.add_argument("--control-pairs", type=Path, required=True)
    parser.add_argument(
        "--modality", choices=("sequence", "text", "structure", "ppi"), required=True
    )
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    requested = load_ids(args.requested_pairs, args.modality)
    controls = load_ids(args.control_pairs, args.modality)
    if requested & controls:
        raise SystemExit("Requested and subset-equivalence control IDs overlap")
    selected = requested | controls
    if not selected:
        raise SystemExit("Retry workspace selection is empty")

    all_sequences: Dict[str, str] = {}
    memberships = {}
    seen = set()
    for aspect in ASPECTS:
        for split in SPLITS:
            prefix = f"{aspect}_{split}"
            names_path = args.data_dir / f"{prefix}_names.npy"
            sequences_path = args.data_dir / f"{prefix}_sequences.json"
            labels_path = args.data_dir / f"{prefix}_labels.npz"
            if not names_path.is_file() or not sequences_path.is_file():
                raise SystemExit(f"Missing prepared PFP split view: {prefix}")
            names = np.load(names_path, allow_pickle=True)
            sequences = json.loads(sequences_path.read_text(encoding="utf-8"))
            control_indices = [
                index for index, value in enumerate(names) if str(value) in controls
            ]
            requested_indices = [
                index for index, value in enumerate(names) if str(value) in requested
            ]
            indices = np.asarray(control_indices + requested_indices, dtype=int)
            filtered_names = names[indices]
            filtered_sequences = {}
            for value in filtered_names:
                protein_id = str(value)
                sequence = sequences[protein_id]
                previous = all_sequences.get(protein_id)
                if previous is not None and previous != sequence:
                    raise SystemExit(f"Conflicting sequence for {protein_id}")
                all_sequences[protein_id] = sequence
                filtered_sequences[protein_id] = sequence
                seen.add(protein_id)
            np.save(names_path, filtered_names)
            sequences_path.write_text(
                json.dumps(filtered_sequences, sort_keys=True), encoding="utf-8"
            )
            if labels_path.is_file():
                labels = sparse.load_npz(labels_path)
                if labels.shape[0] != len(names):
                    raise SystemExit(f"Label/name row mismatch before retry filtering: {prefix}")
                sparse.save_npz(labels_path, labels[indices])
            memberships[prefix] = len(filtered_names)

    missing = sorted(selected - seen)
    if missing:
        raise SystemExit(
            f"Retry pair IDs are absent from prepared PFP data: {missing[:10]} "
            f"(total={len(missing)})"
        )

    fasta_path = args.data_dir / "proteins.fasta"
    with fasta_path.open("w", encoding="ascii") as handle:
        ordered_ids = sorted(controls) + sorted(requested)
        for protein_id in ordered_ids:
            if protein_id not in all_sequences:
                continue
            handle.write(f">{protein_id}\n{all_sequences[protein_id]}\n")

    report = {
        "schema_version": 1,
        "modality": args.modality,
        "requested_count": len(requested),
        "control_count": len(controls),
        "unique_protein_count": len(all_sequences),
        "memberships": memberships,
        "proteins_fasta_sha256": sha256_file(fasta_path),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
