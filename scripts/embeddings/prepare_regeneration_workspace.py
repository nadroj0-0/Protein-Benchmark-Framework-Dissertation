#!/usr/bin/env python3
"""Build the minimal PFP data workspace for a regeneration plan."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


ASPECTS = {"bp": "BPO", "cc": "CCO", "mf": "MFO"}
SPLITS = {"training": "train", "validation": "valid", "test": "test"}
EXPECTED_CSVS = tuple(
    f"{aspect}-{split}.csv"
    for aspect in ASPECTS
    for split in SPLITS
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_target_benchmark(plan_dir: Path, benchmark_dir: Path) -> dict[str, str]:
    manifest_path = plan_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise ValueError(f"Missing planner manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    target = manifest.get("benchmarks", {}).get("target", {})
    records = target.get("input_csvs", [])
    expected = {record["relative_path"]: record["sha256"] for record in records}
    if set(expected) != set(EXPECTED_CSVS):
        raise ValueError("Planner manifest does not describe the required nine target CSVs")

    observed = {}
    for name in EXPECTED_CSVS:
        path = benchmark_dir / name
        if not path.is_file():
            raise ValueError(f"Target benchmark is missing {name}: {benchmark_dir}")
        observed[name] = sha256_file(path)
        if observed[name] != expected[name]:
            raise ValueError(
                f"Target CSV does not match the reuse plan for {name}: "
                f"{observed[name]} != {expected[name]}"
            )
    return observed


def load_regeneration_rows(plan_dir: Path) -> list[dict[str, str]]:
    path = plan_dir / "regenerate_proteins.tsv"
    if not path.is_file():
        raise ValueError(f"Missing regeneration table: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Regeneration plan is empty")

    ids = set()
    for row in rows:
        protein_id = row["protein_id"]
        if (
            not protein_id
            or protein_id in {".", ".."}
            or Path(protein_id).name != protein_id
            or "/" in protein_id
            or "\\" in protein_id
            or any(character.isspace() for character in protein_id)
        ):
            raise ValueError(f"Unsafe protein ID for PFP workspace: {protein_id!r}")
        if protein_id in ids:
            raise ValueError(f"Duplicate regeneration protein: {protein_id}")
        ids.add(protein_id)
        if row["action"] != "regenerate":
            raise ValueError(f"Unexpected action for {protein_id}: {row['action']}")
        sequence = row["sequence"]
        observed = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if observed != row["sequence_sha256"]:
            raise ValueError(f"Sequence digest mismatch for {protein_id}")
    return rows


def split_for_row(row: dict[str, str]) -> str:
    memberships = json.loads(row["target_memberships"])
    splits = {name.removesuffix(".csv").rsplit("-", 1)[1] for name in memberships}
    if len(splits) != 1:
        raise ValueError(
            f"Protein {row['protein_id']} does not have one global target split: {memberships}"
        )
    split = next(iter(splits))
    if split not in SPLITS:
        raise ValueError(f"Unknown split for {row['protein_id']}: {split}")
    return split


def select_rows(
    rows: list[dict[str, str]], limit_per_split: int | None
) -> list[dict[str, str]]:
    if limit_per_split is None:
        return rows
    if limit_per_split <= 0:
        raise ValueError("--limit-per-split must be positive")
    selected = []
    counts = defaultdict(int)
    for row in rows:
        split = split_for_row(row)
        if counts[split] < limit_per_split:
            selected.append(row)
            counts[split] += 1
    missing = [split for split in SPLITS if counts[split] == 0]
    if missing:
        raise ValueError(f"Preflight selection has no proteins for: {', '.join(missing)}")
    return selected


def write_workspace(rows: list[dict[str, str]], data_dir: Path) -> dict:
    data_dir.mkdir(parents=True, exist_ok=True)
    # Preserve modality caches and downloaded runtime resources between a
    # preflight and the full run; replace only the split views consumed by the
    # PFP embedding scripts.
    for aspect in ASPECTS.values():
        for split in SPLITS.values():
            for suffix in ("names.npy", "sequences.json"):
                path = data_dir / f"{aspect}_{split}_{suffix}"
                if path.exists():
                    path.unlink()
    for name in ("proteins.fasta", "regeneration_workspace_manifest.json"):
        path = data_dir / name
        if path.exists():
            path.unlink()

    memberships: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    all_sequences = {}
    global_splits = defaultdict(int)

    for row in rows:
        protein_id = row["protein_id"]
        sequence = row["sequence"]
        all_sequences[protein_id] = sequence
        global_splits[split_for_row(row)] += 1
        for name in json.loads(row["target_memberships"]):
            stem = name.removesuffix(".csv")
            aspect, split = stem.split("-", 1)
            if aspect not in ASPECTS or split not in SPLITS:
                raise ValueError(f"Unknown target membership for {protein_id}: {name}")
            memberships[(ASPECTS[aspect], SPLITS[split])].append((protein_id, sequence))

    for aspect in ASPECTS.values():
        for split in SPLITS.values():
            entries = memberships[(aspect, split)]
            ids = [protein_id for protein_id, _ in entries]
            sequences = dict(entries)
            np.save(data_dir / f"{aspect}_{split}_names.npy", np.asarray(ids, dtype=object))
            (data_dir / f"{aspect}_{split}_sequences.json").write_text(
                json.dumps(sequences, sort_keys=True), encoding="utf-8"
            )

    fasta_path = data_dir / "proteins.fasta"
    with fasta_path.open("w", encoding="ascii") as handle:
        for protein_id in sorted(all_sequences):
            handle.write(f">{protein_id}\n{all_sequences[protein_id]}\n")

    manifest = {
        "schema_version": 1,
        "protein_count": len(all_sequences),
        "global_split_counts": dict(sorted(global_splits.items())),
        "membership_counts": {
            f"{aspect}_{split}": len(memberships[(aspect, split)])
            for aspect in ASPECTS.values()
            for split in SPLITS.values()
        },
        "proteins_fasta": {
            "path": "proteins.fasta",
            "sha256": sha256_file(fasta_path),
        },
    }
    (data_dir / "regeneration_workspace_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--target-benchmark-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--limit-per-split", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    checksums = validate_target_benchmark(args.plan_dir, args.target_benchmark_dir)
    rows = select_rows(load_regeneration_rows(args.plan_dir), args.limit_per_split)
    manifest = write_workspace(rows, args.data_dir)
    manifest["target_csv_sha256"] = checksums
    manifest["reuse_plan_dir"] = str(args.plan_dir.resolve())
    manifest["target_benchmark_dir"] = str(args.target_benchmark_dir.resolve())
    (args.data_dir / "regeneration_workspace_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
