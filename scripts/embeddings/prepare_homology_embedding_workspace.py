#!/usr/bin/env python3
"""Build one exact PFP workspace from a source-resolved homology ledger."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


MODALITIES = ("sequence", "text", "structure", "ppi")
ASPECTS = {"bp": "BPO", "cc": "CCO", "mf": "MFO"}
SPLITS = {"training": "train", "validation": "valid", "test": "test"}
EXPECTED_CSVS = tuple(
    f"{aspect}-{split}.csv" for aspect in ASPECTS for split in SPLITS
)


class WorkspaceError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def verify_ledger(ledger_dir: Path) -> dict:
    manifest_path = ledger_dir / "output_manifest.json"
    complete_path = ledger_dir / "RUN_COMPLETE.json"
    if not manifest_path.is_file() or not complete_path.is_file():
        raise WorkspaceError("Ledger lacks output_manifest.json or RUN_COMPLETE.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for record in manifest.get("files", []):
        path = ledger_dir / record["path"]
        if not path.is_file():
            raise WorkspaceError(f"Ledger payload is missing: {record['path']}")
        if path.stat().st_size != int(record["bytes"]):
            raise WorkspaceError(f"Ledger payload size changed: {record['path']}")
        if sha256_file(path) != record["sha256"]:
            raise WorkspaceError(f"Ledger payload hash changed: {record['path']}")
    complete = json.loads(complete_path.read_text(encoding="utf-8"))
    if complete.get("complete") is not True:
        raise WorkspaceError("Ledger completion marker is not complete")
    if complete.get("output_manifest_sha256") != sha256_file(manifest_path):
        raise WorkspaceError("Ledger completion marker has the wrong manifest hash")
    return manifest


def load_benchmark(benchmark_dir: Path) -> tuple[dict[str, str], dict[str, list[str]], dict]:
    sequences: dict[str, str] = {}
    memberships: dict[str, list[str]] = defaultdict(list)
    csv_records = []
    for name in EXPECTED_CSVS:
        path = benchmark_dir / name
        if not path.is_file():
            raise WorkspaceError(f"Target benchmark is missing {name}")
        seen: set[str] = set()
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or reader.fieldnames[:2] != [
                "proteins",
                "sequences",
            ]:
                raise WorkspaceError(f"{name} must begin with proteins,sequences")
            for line_number, row in enumerate(reader, start=2):
                protein_id = row["proteins"]
                sequence = row["sequences"]
                if not protein_id or Path(protein_id).name != protein_id:
                    raise WorkspaceError(f"Unsafe protein ID at {name}:{line_number}")
                try:
                    sequence.encode("ascii")
                except UnicodeEncodeError as error:
                    raise WorkspaceError(
                        f"Non-ASCII sequence for {protein_id} at {name}:{line_number}"
                    ) from error
                if protein_id in seen:
                    raise WorkspaceError(f"Duplicate protein in {name}: {protein_id}")
                seen.add(protein_id)
                previous = sequences.get(protein_id)
                if previous is not None and previous != sequence:
                    raise WorkspaceError(f"Conflicting benchmark sequence for {protein_id}")
                sequences[protein_id] = sequence
                memberships[protein_id].append(name)
        csv_records.append(
            {"name": name, "rows": len(seen), "sha256": sha256_file(path)}
        )
    for protein_id in memberships:
        memberships[protein_id].sort()
    return sequences, dict(memberships), {"csvs": csv_records}


def load_plan_rows(ledger_dir: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    for action in ("reuse", "regenerate"):
        path = ledger_dir / f"{action}_proteins.tsv"
        if not path.is_file():
            raise WorkspaceError(f"Ledger action table is missing: {path.name}")
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "protein_id",
                "sequence",
                "sequence_sha256",
                "action",
                "regenerate_modalities",
                "target_memberships",
            }
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise WorkspaceError(f"Incomplete action-table schema: {path}")
            for row in reader:
                protein_id = row["protein_id"]
                if protein_id in rows:
                    raise WorkspaceError(f"Protein occurs in both action tables: {protein_id}")
                if row["action"] != action:
                    raise WorkspaceError(f"Action mismatch for {protein_id}")
                if sequence_sha256(row["sequence"]) != row["sequence_sha256"]:
                    raise WorkspaceError(f"Sequence hash mismatch for {protein_id}")
                modalities = json.loads(row["regenerate_modalities"])
                memberships = json.loads(row["target_memberships"])
                if modalities != sorted(set(modalities)) or not set(modalities) <= set(MODALITIES):
                    raise WorkspaceError(f"Invalid regenerate modalities for {protein_id}")
                if memberships != sorted(set(memberships)):
                    raise WorkspaceError(f"Invalid target memberships for {protein_id}")
                rows[protein_id] = row
    return rows


def validate_against_benchmark(
    plan_rows: dict[str, dict[str, str]],
    sequences: dict[str, str],
    memberships: dict[str, list[str]],
) -> None:
    if set(plan_rows) != set(sequences):
        missing = sorted(set(sequences) - set(plan_rows))[:5]
        extra = sorted(set(plan_rows) - set(sequences))[:5]
        raise WorkspaceError(
            f"Ledger/benchmark target mismatch: missing={missing}, extra={extra}"
        )
    for protein_id, sequence in sequences.items():
        row = plan_rows[protein_id]
        if row["sequence"] != sequence:
            raise WorkspaceError(f"Ledger sequence differs for {protein_id}")
        if json.loads(row["target_memberships"]) != memberships[protein_id]:
            raise WorkspaceError(f"Ledger memberships differ for {protein_id}")


def effective_temporal_role(protein_id: str, names: list[str]) -> str:
    splits = {name.removesuffix(".csv").rsplit("-", 1)[1] for name in names}
    unknown = splits - set(SPLITS)
    if unknown:
        raise WorkspaceError(f"Unknown splits for {protein_id}: {sorted(unknown)}")
    if "test" in splits:
        return "test"
    if "validation" in splits:
        return "validation"
    if "training" in splits:
        return "training"
    raise WorkspaceError(f"Protein has no benchmark membership: {protein_id}")


def select_rows(
    plan_rows: dict[str, dict[str, str]], modality: str, limit_per_split: int | None
) -> list[dict[str, str]]:
    selected = [
        row
        for row in plan_rows.values()
        if modality in json.loads(row["regenerate_modalities"])
    ]
    selected.sort(key=lambda row: row["protein_id"])
    if not selected:
        raise WorkspaceError(f"Ledger selects no {modality} pairs for regeneration")
    if limit_per_split is None:
        return selected
    if limit_per_split <= 0:
        raise WorkspaceError("--limit-per-split must be positive")
    counts: dict[str, int] = defaultdict(int)
    bounded = []
    for row in selected:
        split = effective_temporal_role(
            row["protein_id"], json.loads(row["target_memberships"])
        )
        if counts[split] < limit_per_split:
            bounded.append(row)
            counts[split] += 1
    missing = sorted(set(SPLITS) - set(counts))
    if missing:
        raise WorkspaceError(f"Preflight selection lacks global splits: {missing}")
    return bounded


def clear_split_views(data_dir: Path) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    for aspect in ASPECTS.values():
        for split in SPLITS.values():
            for suffix in ("names.npy", "sequences.json"):
                (data_dir / f"{aspect}_{split}_{suffix}").unlink(missing_ok=True)
    (data_dir / "proteins.fasta").unlink(missing_ok=True)


def write_workspace(rows: list[dict[str, str]], data_dir: Path) -> dict:
    clear_split_views(data_dir)
    membership_rows: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    sequences: dict[str, str] = {}
    temporal_roles: dict[str, int] = defaultdict(int)
    mixed_aspect_split_proteins = 0
    for row in rows:
        protein_id = row["protein_id"]
        sequence = row["sequence"]
        sequences[protein_id] = sequence
        names = json.loads(row["target_memberships"])
        observed_splits = {
            name.removesuffix(".csv").rsplit("-", 1)[1] for name in names
        }
        mixed_aspect_split_proteins += int(len(observed_splits) > 1)
        temporal_roles[effective_temporal_role(protein_id, names)] += 1
        for name in names:
            aspect, split = name.removesuffix(".csv").split("-", 1)
            membership_rows[(ASPECTS[aspect], SPLITS[split])].append(
                (protein_id, sequence)
            )

    for aspect in ASPECTS.values():
        for split in SPLITS.values():
            entries = membership_rows[(aspect, split)]
            np.save(
                data_dir / f"{aspect}_{split}_names.npy",
                np.asarray([protein_id for protein_id, _ in entries], dtype=object),
            )
            (data_dir / f"{aspect}_{split}_sequences.json").write_text(
                json.dumps(dict(entries), sort_keys=True), encoding="utf-8"
            )

    fasta_path = data_dir / "proteins.fasta"
    with fasta_path.open("w", encoding="ascii") as handle:
        for protein_id in sorted(sequences):
            handle.write(f">{protein_id}\n{sequences[protein_id]}\n")
    return {
        "protein_count": len(sequences),
        "effective_temporal_role_counts": dict(sorted(temporal_roles.items())),
        "mixed_aspect_split_proteins": mixed_aspect_split_proteins,
        "shared_text_temporal_role_policy": "historical-if-test-in-any-aspect",
        "proteins_fasta_sha256": sha256_file(fasta_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger-dir", type=Path, required=True)
    parser.add_argument("--target-benchmark-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--modality", choices=MODALITIES, required=True)
    parser.add_argument("--limit-per-split", type=int)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    ledger_manifest = verify_ledger(args.ledger_dir)
    sequences, memberships, benchmark = load_benchmark(args.target_benchmark_dir)
    plan_rows = load_plan_rows(args.ledger_dir)
    validate_against_benchmark(plan_rows, sequences, memberships)
    rows = select_rows(plan_rows, args.modality, args.limit_per_split)
    report = {
        "schema_version": 1,
        "modality": args.modality,
        "selection": "preflight" if args.limit_per_split is not None else "full",
        "ledger_dir": str(args.ledger_dir.resolve()),
        "ledger_output_manifest_sha256": sha256_file(
            args.ledger_dir / "output_manifest.json"
        ),
        "ledger_payload_file_count": ledger_manifest["payload_file_count"],
        "target_benchmark_dir": str(args.target_benchmark_dir.resolve()),
        **benchmark,
        **write_workspace(rows, args.data_dir),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, json.JSONDecodeError, WorkspaceError) as error:
        raise SystemExit(f"ERROR: {error}") from error
