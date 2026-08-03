#!/usr/bin/env python3
"""Bind a hydrated PFP embedding archive to an exact nine-CSV benchmark."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import tarfile
import tempfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np


EXPECTED_CSVS = tuple(
    f"{aspect}-{split}.csv"
    for aspect in ("bp", "cc", "mf")
    for split in ("training", "validation", "test")
)
PFP_COMMIT = "1e04fd6d6d3c40458fd41ec1a881ed6e24de768e"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return sha256_text(payload)


def safe_protein_id(protein_id: str) -> bool:
    return bool(
        protein_id
        and protein_id not in {".", ".."}
        and PurePosixPath(protein_id).name == protein_id
        and not any(character.isspace() for character in protein_id)
    )


def atomic_write(path: Path, content: str) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_json(path: Path, value: Any) -> None:
    atomic_write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def benchmark_files(benchmark_dir: Path) -> list[dict[str, Any]]:
    result = []
    for name in EXPECTED_CSVS:
        path = benchmark_dir / name
        if not path.is_file():
            raise ValueError(f"Missing benchmark CSV: {path}")
        result.append(
            {"name": name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        )
    return result


def load_targets(benchmark_dir: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for name in EXPECTED_CSVS:
        path = benchmark_dir / name
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValueError(f"Benchmark CSV has no header: {path}")
            id_column = "proteins" if "proteins" in reader.fieldnames else "protein"
            if id_column not in reader.fieldnames or "sequences" not in reader.fieldnames:
                raise ValueError(f"Benchmark CSV lacks protein/sequence columns: {path}")
            for line_number, row in enumerate(reader, start=2):
                protein_id = row[id_column]
                sequence = row["sequences"]
                if not safe_protein_id(protein_id):
                    raise ValueError(f"Unsafe protein ID at {path}:{line_number}: {protein_id!r}")
                if not sequence:
                    raise ValueError(f"Empty sequence at {path}:{line_number}")
                try:
                    sequence.encode("ascii")
                except UnicodeEncodeError as error:
                    raise ValueError(
                        f"Non-ASCII sequence at {path}:{line_number}"
                    ) from error
                previous = sequences.get(protein_id)
                if previous is not None and previous != sequence:
                    raise ValueError(f"Conflicting sequences for protein {protein_id}")
                sequences[protein_id] = sequence
    if not sequences:
        raise ValueError("Benchmark contains no proteins")
    return {
        protein_id: sha256_text(sequence)
        for protein_id, sequence in sorted(sequences.items())
    }


def load_policy(config_path: Path) -> dict[str, Any]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    modalities = config.get("modalities")
    expected = {"sequence", "text", "structure", "ppi"}
    if not isinstance(modalities, dict) or set(modalities) != expected:
        raise ValueError("Run config must define exactly four modalities")
    policy_modalities = {}
    seen_directories: set[str] = set()
    for modality in sorted(expected):
        specification = modalities[modality]
        directory = Path(str(specification.get("directory", ""))).name
        dimension = int(specification.get("dimension", 0))
        if not directory or directory in seen_directories or dimension <= 0:
            raise ValueError(f"Invalid run-config modality: {modality}")
        seen_directories.add(directory)
        policy_modalities[modality] = {
            "cache_directory": directory,
            "dimension": dimension,
        }
    return {
        "schema_version": 1,
        "name": "observed-hydrated-pfp-cache",
        "modalities": policy_modalities,
    }


def inspect_array(data: bytes, dimension: int, member_name: str) -> None:
    try:
        value = np.load(io.BytesIO(data), allow_pickle=False)
    except Exception as error:
        raise ValueError(f"Unreadable array {member_name}: {error}") from error
    if value.shape != (dimension,):
        raise ValueError(
            f"Wrong array shape for {member_name}: {tuple(value.shape)} != {(dimension,)}"
        )
    if value.dtype.kind != "f":
        raise ValueError(f"Unsupported array dtype for {member_name}: {value.dtype}")
    if not np.isfinite(value).all():
        raise ValueError(f"Non-finite array values in {member_name}")
    with np.errstate(over="ignore", invalid="ignore"):
        converted = value.astype(np.float32)
    if not np.isfinite(converted).all():
        raise ValueError(f"Array becomes non-finite as float32: {member_name}")


def scan_archive(
    archive_path: Path,
    targets: dict[str, str],
    policy: dict[str, Any],
) -> tuple[dict[tuple[str, str], str], dict[str, int]]:
    directory_to_modality = {
        specification["cache_directory"]: modality
        for modality, specification in policy["modalities"].items()
    }
    accepted: dict[tuple[str, str], str] = {}
    ignored_arrays = 0
    inspected_members = 0
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile():
                continue
            parts = PurePosixPath(member.name).parts
            if len(parts) < 2 or not parts[-1].endswith(".npy"):
                continue
            modality = directory_to_modality.get(parts[-2])
            if modality is None:
                continue
            protein_id = parts[-1][:-4]
            if protein_id not in targets:
                ignored_arrays += 1
                continue
            key = (protein_id, modality)
            if key in accepted:
                raise ValueError(f"Archive repeats target/modality pair {key}")
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValueError(f"Cannot read archive member: {member.name}")
            data = extracted.read()
            if len(data) != member.size:
                raise ValueError(f"Short archive member read: {member.name}")
            inspect_array(data, int(policy["modalities"][modality]["dimension"]), member.name)
            accepted[key] = sha256_bytes(data)
            inspected_members += 1
    sequence_missing = sorted(
        protein_id for protein_id in targets if (protein_id, "sequence") not in accepted
    )
    if sequence_missing:
        sample = ", ".join(sequence_missing[:10])
        raise ValueError(
            f"Sequence embeddings are missing for {len(sequence_missing)} targets; sample: {sample}"
        )
    return accepted, {
        "accepted_arrays": len(accepted),
        "ignored_non_target_arrays": ignored_arrays,
        "inspected_members": inspected_members,
    }


def targets_content(targets: dict[str, str]) -> str:
    lines = ["protein_id\tsequence_sha256"]
    lines.extend(f"{protein_id}\t{digest}" for protein_id, digest in targets.items())
    return "\n".join(lines) + "\n"


def pair_status_content(
    targets: dict[str, str],
    accepted: dict[tuple[str, str], str],
) -> str:
    lines = [
        "protein_id\tmodality\tstate\tsequence_sha256\tembedding_sha256\t"
        "attempts\tlatest_reason\tlatest_detail"
    ]
    for protein_id, sequence_sha in targets.items():
        for modality in ("sequence", "text", "structure", "ppi"):
            embedding_sha = accepted.get((protein_id, modality), "")
            lines.append(
                "\t".join(
                    (
                        protein_id,
                        modality,
                        "accepted" if embedding_sha else "needs_retry",
                        sequence_sha,
                        embedding_sha,
                        "0",
                        "" if embedding_sha else "missing_from_hydrated_archive",
                        "",
                    )
                )
            )
    return "\n".join(lines) + "\n"


def publish(args: argparse.Namespace) -> dict[str, Any]:
    benchmark_dir = args.benchmark_dir.resolve()
    archive_path = args.archive.resolve()
    config_path = args.config.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")
    if not archive_path.is_file() or not config_path.is_file():
        raise ValueError("Archive and config must be regular files")
    if len(args.framework_commit) != 40 or len(args.pfp_commit) != 40:
        raise ValueError("Framework and PFP commits must be full 40-character hashes")

    policy = load_policy(config_path)
    csv_records = benchmark_files(benchmark_dir)
    targets = load_targets(benchmark_dir)
    archive_sha = sha256_file(archive_path)
    if args.archive_sha256 and archive_sha != args.archive_sha256:
        raise ValueError(
            f"Archive SHA-256 mismatch: expected {args.archive_sha256}, observed {archive_sha}"
        )
    accepted, scan_summary = scan_archive(archive_path, targets, policy)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.staging-", dir=output_dir.parent))
    try:
        target_text = targets_content(targets)
        atomic_write(stage / "targets.tsv", target_text)
        atomic_write(stage / "pair_status.tsv", pair_status_content(targets, accepted))
        contract = {
            "schema_version": 1,
            "benchmark_id": args.benchmark_id,
            "benchmark_csvs": csv_records,
            "targets": {
                "count": len(targets),
                "manifest_sha256": sha256_text(target_text),
            },
            "pfp_commit": args.pfp_commit,
            "framework_commit": args.framework_commit,
            "policy": policy,
            "policy_sha256": canonical_sha256(policy),
            "environment": None,
            "source_files": [
                {
                    "label": "hydrated_embedding_archive",
                    "name": archive_path.name,
                    "path": str(archive_path),
                    "sha256": archive_sha,
                    "size_bytes": archive_path.stat().st_size,
                },
                {
                    "label": "model_execution_config",
                    "name": config_path.name,
                    "path": str(config_path),
                    "sha256": sha256_file(config_path),
                    "size_bytes": config_path.stat().st_size,
                },
            ],
            "runtime": {"bound_at_utc": utc_now()},
        }
        contract["contract_sha256"] = canonical_sha256(contract)
        atomic_json(stage / "contract.json", contract)

        coverage = {}
        for modality in ("sequence", "text", "structure", "ppi"):
            accepted_count = sum(1 for protein_id in targets if (protein_id, modality) in accepted)
            coverage[modality] = {
                "accepted": accepted_count,
                "needs_retry": len(targets) - accepted_count,
                "fraction": accepted_count / len(targets),
            }
        atomic_json(
            stage / "coverage.json",
            {
                "schema_version": 1,
                "contract_sha256": contract["contract_sha256"],
                "target_count": len(targets),
                "coverage": coverage,
            },
        )
        summary = {
            "schema_version": 1,
            "complete": True,
            "benchmark_id": args.benchmark_id,
            "benchmark_dir": str(benchmark_dir),
            "archive": str(archive_path),
            "archive_sha256": archive_sha,
            "target_count": len(targets),
            "coverage": coverage,
            "scan": scan_summary,
            "contract_sha256": contract["contract_sha256"],
        }
        atomic_json(stage / "summary.json", summary)
        files = []
        for name in (
            "contract.json",
            "coverage.json",
            "targets.tsv",
            "pair_status.tsv",
            "summary.json",
        ):
            path = stage / name
            files.append(
                {"path": name, "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
            )
        atomic_json(stage / "output_manifest.json", {"schema_version": 1, "files": files})
        manifest = stage / "output_manifest.json"
        atomic_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "output_manifest": {
                    "path": manifest.name,
                    "size_bytes": manifest.stat().st_size,
                    "sha256": sha256_file(manifest),
                },
            },
        )
        os.replace(stage, output_dir)
        return summary
    except Exception:
        if stage.exists():
            for child in stage.iterdir():
                child.unlink()
            stage.rmdir()
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--archive-sha256")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--framework-commit", required=True)
    parser.add_argument("--pfp-commit", default=PFP_COMMIT)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(publish(args), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
