#!/usr/bin/env python3
"""Bind a CAFA3 reconstruction to reusable and newly generated ProtT5 arrays."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


CSV_NAMES = tuple(
    f"{aspect}-{split}.csv"
    for aspect in ("bp", "cc", "mf")
    for split in ("training", "validation", "test")
)
SAFE_ID = re.compile(r"^[^\s/\\]+$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
MODALITIES = {
    "sequence": ("prott5", 1024),
    "text": ("exp_text_embeddings_temporal", 768),
    "structure": ("IF1", 512),
    "ppi": ("ppi", 512),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass


def atomic_json(path: Path, payload: object) -> None:
    atomic_write(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def validate_id(protein_id: str) -> None:
    if not SAFE_ID.fullmatch(protein_id) or protein_id in {".", ".."}:
        raise ValueError(f"Unsafe protein ID: {protein_id!r}")


def read_checksum_ledger(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = raw.split(maxsplit=1)
        if len(fields) != 2 or not SHA256.fullmatch(fields[0]):
            raise ValueError(f"Malformed checksum ledger row {line_number}: {path}")
        name = Path(fields[1]).name
        if name in CSV_NAMES:
            if name in result:
                raise ValueError(f"Checksum ledger repeats {name}")
            result[name] = fields[0]
    if set(result) != set(CSV_NAMES):
        raise ValueError("Checksum ledger does not bind exactly the nine CAFA3 CSVs")
    return result


def read_benchmark(
    benchmark_dir: Path, checksum_ledger: Path
) -> tuple[dict[str, str], dict[str, str]]:
    expected = read_checksum_ledger(checksum_ledger)
    observed: dict[str, str] = {}
    sequences: dict[str, str] = {}
    for name in CSV_NAMES:
        path = benchmark_dir / name
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Benchmark CSV is missing or unsafe: {path}")
        observed[name] = sha256_file(path)
        if observed[name] != expected[name]:
            raise ValueError(f"Benchmark CSV checksum mismatch: {name}")
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or len(reader.fieldnames) < 3:
                raise ValueError(f"Benchmark CSV has an invalid header: {name}")
            if reader.fieldnames[0] not in {"protein", "proteins"}:
                raise ValueError(f"Benchmark CSV has an invalid protein column: {name}")
            if reader.fieldnames[1] != "sequences":
                raise ValueError(
                    f"Benchmark CSV has an invalid sequence column: {name}"
                )
            protein_column = reader.fieldnames[0]
            for row in reader:
                protein_id = row[protein_column]
                sequence = row["sequences"]
                validate_id(protein_id)
                if not sequence:
                    raise ValueError(f"Empty sequence for {protein_id}")
                previous = sequences.setdefault(protein_id, sequence)
                if previous != sequence:
                    raise ValueError(f"Conflicting benchmark sequence for {protein_id}")
    if not sequences:
        raise ValueError("Benchmark contains no proteins")
    return sequences, observed


def read_source_targets(path: Path, expected_sha256: str) -> dict[str, str]:
    if sha256_file(path) != expected_sha256:
        raise ValueError("Canonical source target manifest checksum mismatch")
    result: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["protein_id", "sequence_sha256"]:
            raise ValueError("Canonical source target manifest has an invalid schema")
        for row in reader:
            protein_id = row["protein_id"]
            digest = row["sequence_sha256"]
            validate_id(protein_id)
            if not SHA256.fullmatch(digest):
                raise ValueError(f"Invalid source sequence hash for {protein_id}")
            if protein_id in result:
                raise ValueError(
                    f"Canonical source target manifest repeats {protein_id}"
                )
            result[protein_id] = digest
    return result


def inspect_array(path: Path, dimension: int) -> str:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"Embedding array is missing or unsafe: {path}")
    try:
        value = np.load(path, allow_pickle=False)
    except Exception as exc:
        raise ValueError(f"Cannot load embedding array {path}: {exc}") from exc
    if value.shape != (dimension,):
        raise ValueError(f"Unexpected embedding shape for {path}: {value.shape}")
    if value.dtype.kind != "f" or not np.isfinite(value).all():
        raise ValueError(f"Embedding array is not finite floating-point data: {path}")
    return sha256_file(path)


def write_actions(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ("protein_id", "sequence_sha256", "action", "source_protein_id")
    lines = ["\t".join(fields)]
    for row in rows:
        lines.append("\t".join(row[field] for field in fields))
    atomic_write(path, "\n".join(lines) + "\n")


def write_fasta(path: Path, proteins: list[tuple[str, str]]) -> None:
    lines: list[str] = []
    for protein_id, sequence in proteins:
        lines.extend((f">{protein_id}", sequence))
    atomic_write(path, "\n".join(lines) + ("\n" if lines else ""))


def prepare(args: argparse.Namespace) -> int:
    sequences, csv_hashes = read_benchmark(args.benchmark_dir, args.benchmark_checksums)
    source_targets = read_source_targets(
        args.source_targets, args.expected_source_targets_sha256
    )
    source_by_sequence: dict[str, list[str]] = defaultdict(list)
    for protein_id, digest in source_targets.items():
        source_by_sequence[digest].append(protein_id)
    for protein_ids in source_by_sequence.values():
        protein_ids.sort()

    sequence_dir = args.cache_root / MODALITIES["sequence"][0]
    if not sequence_dir.is_dir() or sequence_dir.is_symlink():
        raise ValueError(
            f"Canonical sequence cache is missing or unsafe: {sequence_dir}"
        )

    rows: list[dict[str, str]] = []
    copies: list[tuple[Path, Path]] = []
    missing: list[tuple[str, str]] = []
    counts: Counter[str] = Counter()
    for protein_id in sorted(sequences):
        sequence = sequences[protein_id]
        digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
        source_digest = source_targets.get(protein_id)
        if source_digest is not None and source_digest != digest:
            raise ValueError(
                f"Canonical source uses a different sequence for benchmark ID {protein_id}"
            )
        if source_digest == digest:
            source_id = protein_id
            action = "direct-reuse"
            inspect_array(sequence_dir / f"{source_id}.npy", 1024)
        elif source_by_sequence.get(digest):
            source_id = source_by_sequence[digest][0]
            action = "exact-sequence-reuse"
            source_path = sequence_dir / f"{source_id}.npy"
            destination = sequence_dir / f"{protein_id}.npy"
            inspect_array(source_path, 1024)
            if destination.exists():
                raise ValueError(
                    f"Unexpected pre-existing cross-ID destination: {destination}"
                )
            copies.append((source_path, destination))
        else:
            source_id = ""
            action = "generate"
            missing.append((protein_id, sequence))
        rows.append(
            {
                "protein_id": protein_id,
                "sequence_sha256": digest,
                "action": action,
                "source_protein_id": source_id,
            }
        )
        counts[action] += 1

    for source, destination in copies:
        shutil.copy2(source, destination)
    write_actions(args.actions_tsv, rows)
    write_fasta(args.missing_fasta, missing)
    target_manifest = ["protein_id\tsequence_sha256"]
    target_manifest.extend(
        f"{row['protein_id']}\t{row['sequence_sha256']}" for row in rows
    )
    atomic_write(args.target_manifest, "\n".join(target_manifest) + "\n")
    payload = {
        "schema_version": 1,
        "status": "prepared",
        "target_count": len(rows),
        "actions": dict(sorted(counts.items())),
        "benchmark_csv_sha256": dict(sorted(csv_hashes.items())),
        "benchmark_checksum_ledger_sha256": sha256_file(args.benchmark_checksums),
        "source_target_count": len(source_targets),
        "source_targets_sha256": sha256_file(args.source_targets),
        "actions_sha256": sha256_file(args.actions_tsv),
        "target_manifest_sha256": sha256_file(args.target_manifest),
        "missing_fasta_sha256": sha256_file(args.missing_fasta),
        "policy": {
            "direct_reuse": "Protein ID and sequence SHA-256 both match the authenticated source target manifest.",
            "exact_sequence_reuse": "ProtT5 only: sequence SHA-256 matches an authenticated source protein; auxiliary modalities are never copied across IDs.",
            "generate": "No authenticated source protein has the exact sequence; generate a fresh ProtT5 vector.",
        },
    }
    atomic_json(args.report, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


def read_actions(path: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["protein_id", "sequence_sha256", "action", "source_protein_id"]
        if reader.fieldnames != expected:
            raise ValueError("Action table has an invalid schema")
        for row in reader:
            protein_id = row["protein_id"]
            if protein_id in result:
                raise ValueError(f"Action table repeats {protein_id}")
            if row["action"] not in {
                "direct-reuse",
                "exact-sequence-reuse",
                "generate",
            }:
                raise ValueError(f"Action table has an invalid action for {protein_id}")
            result[protein_id] = row
    return result


def finalize(args: argparse.Namespace) -> int:
    sequences, csv_hashes = read_benchmark(args.benchmark_dir, args.benchmark_checksums)
    actions = read_actions(args.actions_tsv)
    if set(actions) != set(sequences):
        raise ValueError("Action table does not match the benchmark population")

    action_counts: Counter[str] = Counter()
    modality_counts: dict[str, int] = {}
    sequence_content = hashlib.sha256()
    for protein_id in sorted(sequences):
        row = actions[protein_id]
        digest = hashlib.sha256(sequences[protein_id].encode("utf-8")).hexdigest()
        if row["sequence_sha256"] != digest:
            raise ValueError(f"Action table sequence hash differs for {protein_id}")
        target = args.cache_root / "prott5" / f"{protein_id}.npy"
        target_sha = inspect_array(target, 1024)
        if row["action"] == "exact-sequence-reuse":
            source_id = row["source_protein_id"]
            validate_id(source_id)
            source = args.cache_root / "prott5" / f"{source_id}.npy"
            if inspect_array(source, 1024) != target_sha:
                raise ValueError(
                    f"Cross-ID ProtT5 copy differs from its source: {protein_id}"
                )
        sequence_content.update(protein_id.encode("utf-8"))
        sequence_content.update(b"\t")
        sequence_content.update(target_sha.encode("ascii"))
        sequence_content.update(b"\n")
        action_counts[row["action"]] += 1

    for modality, (directory, dimension) in MODALITIES.items():
        count = 0
        for protein_id in sorted(sequences):
            path = args.cache_root / directory / f"{protein_id}.npy"
            if path.exists():
                inspect_array(path, dimension)
                count += 1
        modality_counts[modality] = count
    if modality_counts["sequence"] != len(sequences):
        raise ValueError("Sequence cache is incomplete after generation")
    for modality in ("text", "structure", "ppi"):
        if modality_counts[modality] == 0:
            raise ValueError(f"Cache has zero valid {modality} arrays")

    payload = {
        "schema_version": 1,
        "status": "passed",
        "target_count": len(sequences),
        "actions": dict(sorted(action_counts.items())),
        "modality_coverage": {
            modality: {
                "valid": count,
                "missing": len(sequences) - count,
                "fraction": count / len(sequences),
            }
            for modality, count in sorted(modality_counts.items())
        },
        "sequence_array_content_sha256": sequence_content.hexdigest(),
        "benchmark_csv_sha256": dict(sorted(csv_hashes.items())),
        "actions_sha256": sha256_file(args.actions_tsv),
    }
    atomic_json(args.report, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)

    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--benchmark-dir", type=Path, required=True)
    prepare_parser.add_argument("--benchmark-checksums", type=Path, required=True)
    prepare_parser.add_argument("--source-targets", type=Path, required=True)
    prepare_parser.add_argument("--expected-source-targets-sha256", required=True)
    prepare_parser.add_argument("--cache-root", type=Path, required=True)
    prepare_parser.add_argument("--missing-fasta", type=Path, required=True)
    prepare_parser.add_argument("--target-manifest", type=Path, required=True)
    prepare_parser.add_argument("--actions-tsv", type=Path, required=True)
    prepare_parser.add_argument("--report", type=Path, required=True)
    prepare_parser.set_defaults(function=prepare)

    finalize_parser = commands.add_parser("finalize")
    finalize_parser.add_argument("--benchmark-dir", type=Path, required=True)
    finalize_parser.add_argument("--benchmark-checksums", type=Path, required=True)
    finalize_parser.add_argument("--cache-root", type=Path, required=True)
    finalize_parser.add_argument("--actions-tsv", type=Path, required=True)
    finalize_parser.add_argument("--report", type=Path, required=True)
    finalize_parser.set_defaults(function=finalize)
    return root


def main() -> int:
    args = parser().parse_args()
    if not SHA256.fullmatch(getattr(args, "expected_source_targets_sha256", "0" * 64)):
        raise SystemExit("--expected-source-targets-sha256 must be lowercase SHA-256")
    return args.function(args)


if __name__ == "__main__":
    raise SystemExit(main())
