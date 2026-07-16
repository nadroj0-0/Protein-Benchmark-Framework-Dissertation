#!/usr/bin/env python3
"""Temporarily reduce prepared CAFA3 split views for an embedding preflight."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import scipy.sparse as sparse


ASPECTS = ("BPO", "CCO", "MFO")
SPLITS = ("train", "valid", "test")
SUFFIXES = ("names.npy", "labels.npz", "sequences.json")
MANIFEST_NAME = "preflight_backup_manifest.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_paths(data_dir: Path) -> list[Path]:
    return [
        data_dir / f"{aspect}_{split}_{suffix}"
        for aspect in ASPECTS
        for split in SPLITS
        for suffix in SUFFIXES
    ]


def backup_paths(data_dir: Path) -> list[Path]:
    return [*split_paths(data_dir), data_dir / "proteins.fasta"]


def validate_inputs(data_dir: Path) -> None:
    missing = [path.name for path in backup_paths(data_dir) if not path.is_file()]
    if missing:
        raise ValueError("Missing prepared CAFA3 files: " + ", ".join(missing))


def load_split(data_dir: Path, aspect: str, split: str):
    names_path = data_dir / f"{aspect}_{split}_names.npy"
    labels_path = data_dir / f"{aspect}_{split}_labels.npz"
    sequences_path = data_dir / f"{aspect}_{split}_sequences.json"
    names = np.load(names_path, allow_pickle=True)
    labels = sparse.load_npz(labels_path)
    sequences = json.loads(sequences_path.read_text(encoding="utf-8"))
    if labels.shape[0] != len(names):
        raise ValueError(
            f"{aspect}_{split} label rows {labels.shape[0]} != names {len(names)}"
        )
    missing = [str(name) for name in names if str(name) not in sequences]
    if missing:
        raise ValueError(
            f"{aspect}_{split} has names without sequences: {', '.join(missing[:5])}"
        )
    return names, labels, sequences


def write_fasta(data_dir: Path) -> tuple[int, str]:
    sequences: dict[str, str] = {}
    for aspect in ASPECTS:
        for split in SPLITS:
            path = data_dir / f"{aspect}_{split}_sequences.json"
            current = json.loads(path.read_text(encoding="utf-8"))
            for protein_id, sequence in current.items():
                previous = sequences.setdefault(protein_id, sequence)
                if previous != sequence:
                    raise ValueError(
                        f"Conflicting sequences for preflight protein {protein_id}"
                    )

    fasta_path = data_dir / "proteins.fasta"
    with fasta_path.open("w", encoding="ascii") as handle:
        for protein_id in sorted(sequences):
            handle.write(f">{protein_id}\n{sequences[protein_id]}\n")
    return len(sequences), sha256_file(fasta_path)


def create_preflight(data_dir: Path, backup_dir: Path, limit: int) -> dict:
    if limit <= 0:
        raise ValueError("--limit-per-split must be positive")
    validate_inputs(data_dir)
    if backup_dir.exists():
        raise ValueError(f"Backup directory already exists: {backup_dir}")
    backup_dir.mkdir(parents=True)

    files = {}
    for source in backup_paths(data_dir):
        destination = backup_dir / source.name
        shutil.copy2(source, destination)
        files[source.name] = {
            "sha256": sha256_file(destination),
            "size_bytes": destination.stat().st_size,
        }

    counts = {}
    for aspect in ASPECTS:
        for split in SPLITS:
            names, labels, sequences = load_split(data_dir, aspect, split)
            if len(names) == 0:
                raise ValueError(f"Cannot preflight empty split: {aspect}_{split}")
            selected_names = names[:limit]
            selected_ids = [str(value) for value in selected_names]
            selected_sequences = {protein_id: sequences[protein_id] for protein_id in selected_ids}
            np.save(data_dir / f"{aspect}_{split}_names.npy", selected_names)
            sparse.save_npz(
                data_dir / f"{aspect}_{split}_labels.npz", labels[: len(selected_names)]
            )
            (data_dir / f"{aspect}_{split}_sequences.json").write_text(
                json.dumps(selected_sequences, sort_keys=True), encoding="utf-8"
            )
            counts[f"{aspect}_{split}"] = {
                "full": len(names),
                "preflight": len(selected_names),
            }

    protein_count, fasta_sha256 = write_fasta(data_dir)
    manifest = {
        "schema_version": 1,
        "data_dir": str(data_dir.resolve()),
        "files": files,
        "limit_per_split": limit,
        "split_counts": counts,
        "preflight_unique_proteins": protein_count,
        "preflight_fasta_sha256": fasta_sha256,
        "restored": False,
    }
    (backup_dir / MANIFEST_NAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def restore_full(data_dir: Path, backup_dir: Path) -> dict:
    manifest_path = backup_dir / MANIFEST_NAME
    if not manifest_path.is_file():
        raise ValueError(f"Missing preflight backup manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_files = {path.name for path in backup_paths(data_dir)}
    if set(manifest.get("files", {})) != expected_files:
        raise ValueError("Preflight backup manifest has an unexpected file set")

    for name, metadata in manifest["files"].items():
        source = backup_dir / name
        if not source.is_file() or sha256_file(source) != metadata["sha256"]:
            raise ValueError(f"Preflight backup failed authentication: {source}")
        shutil.copy2(source, data_dir / name)

    validate_inputs(data_dir)
    for path in backup_paths(data_dir):
        wanted = manifest["files"][path.name]["sha256"]
        if sha256_file(path) != wanted:
            raise ValueError(f"Restored split file failed authentication: {path}")

    protein_count = sum(
        1
        for line in (data_dir / "proteins.fasta").read_text(encoding="ascii").splitlines()
        if line.startswith(">")
    )
    fasta_sha256 = sha256_file(data_dir / "proteins.fasta")
    manifest["restored"] = True
    manifest["restored_unique_proteins"] = protein_count
    manifest["restored_fasta_sha256"] = fasta_sha256
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("create", "restore"):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--data-dir", type=Path, required=True)
        subparser.add_argument("--backup-dir", type=Path, required=True)
        if command == "create":
            subparser.add_argument("--limit-per-split", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "create":
        result = create_preflight(args.data_dir, args.backup_dir, args.limit_per_split)
    else:
        result = restore_full(args.data_dir, args.backup_dir)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
