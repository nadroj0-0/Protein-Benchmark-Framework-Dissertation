#!/usr/bin/env python3
"""Maintain one provenance-bound, cumulative embedding cache.

The state directory is intentionally independent of PFP. PFP writes into a
disposable scratch cache; this tool validates and atomically publishes only
accepted arrays into persistent storage. Missing or invalid protein/modality
pairs remain in one compact retry ledger regardless of their failure reason.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import math
import os
import shutil
import sys
import tarfile
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np


EXPECTED_CSVS = tuple(
    f"{aspect}-{split}.csv"
    for aspect in ("bp", "cc", "mf")
    for split in ("training", "validation", "test")
)
ASSEMBLY_MODALITIES = {
    "prott5": "sequence",
    "text": "text",
    "structure": "structure",
    "ppi": "ppi",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def safe_protein_id(protein_id: str) -> bool:
    return bool(
        protein_id
        and protein_id not in {".", ".."}
        and Path(protein_id).name == protein_id
        and "/" not in protein_id
        and "\\" not in protein_id
        and not any(character.isspace() for character in protein_id)
    )


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_json(path: Path, value: object) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=str(destination.parent)
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        shutil.copyfile(source, temporary)
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_copy_stream(source, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=str(destination.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            shutil.copyfileobj(source, handle, length=1024 * 1024)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


@contextmanager
def state_lock(state_root: Path) -> Iterator[None]:
    state_root.mkdir(parents=True, exist_ok=True)
    lock_path = state_root / ".state.lock"
    with lock_path.open("a+", encoding="ascii") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def load_policy(path: Path) -> dict:
    policy = json.loads(path.read_text(encoding="utf-8"))
    if policy.get("schema_version") != 1:
        raise ValueError("Embedding policy must use schema_version 1")
    modalities = policy.get("modalities")
    if not isinstance(modalities, dict) or not modalities:
        raise ValueError("Embedding policy has no modalities")
    required = {"sequence", "text", "structure", "ppi"}
    if set(modalities) != required:
        raise ValueError(f"Embedding policy modalities must be exactly {sorted(required)}")
    cache_directories: Set[str] = set()
    for modality, specification in modalities.items():
        directory = specification.get("cache_directory")
        dimension = specification.get("dimension")
        if not isinstance(directory, str) or not directory:
            raise ValueError(f"Missing cache_directory for {modality}")
        if directory in cache_directories:
            raise ValueError(f"Repeated cache directory: {directory}")
        cache_directories.add(directory)
        if not isinstance(dimension, int) or dimension <= 0:
            raise ValueError(f"Invalid dimension for {modality}: {dimension}")
        if "min_accepted_count" not in specification and "min_accepted_fraction" not in specification:
            raise ValueError(f"No acceptance threshold for {modality}")
    return policy


def load_targets(data_dir: Path) -> Dict[str, str]:
    sequence_files = sorted(data_dir.glob("*_sequences.json"))
    if not sequence_files:
        raise ValueError(f"No prepared sequence JSON files under {data_dir}")
    sequences: Dict[str, str] = {}
    for path in sequence_files:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError(f"Sequence JSON is not an object: {path}")
        for raw_id, raw_sequence in payload.items():
            protein_id = str(raw_id)
            sequence = str(raw_sequence)
            if not safe_protein_id(protein_id):
                raise ValueError(f"Unsafe protein ID: {protein_id!r}")
            try:
                sequence.encode("ascii")
            except UnicodeEncodeError as error:
                raise ValueError(f"Non-ASCII sequence for {protein_id}") from error
            previous = sequences.get(protein_id)
            if previous is not None and previous != sequence:
                raise ValueError(f"Conflicting sequences for {protein_id}")
            sequences[protein_id] = sequence
    return {
        protein_id: sha256_text(sequence)
        for protein_id, sequence in sorted(sequences.items())
    }


def load_target_tables(paths: Sequence[Path]) -> Dict[str, str]:
    if not paths:
        raise ValueError("At least one target table is required")
    targets: Dict[str, str] = {}
    for path in paths:
        if not path.is_file():
            raise ValueError(f"Missing target table: {path}")
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"protein_id", "sequence_sha256"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(
                    f"Target table {path} must contain {sorted(required)}"
                )
            for row in reader:
                protein_id = row["protein_id"]
                sequence_sha256 = row["sequence_sha256"]
                if not safe_protein_id(protein_id):
                    raise ValueError(f"Unsafe protein ID: {protein_id!r}")
                if len(sequence_sha256) != 64 or any(
                    character not in "0123456789abcdef" for character in sequence_sha256
                ):
                    raise ValueError(
                        f"Invalid sequence SHA-256 for {protein_id}: {sequence_sha256}"
                    )
                if row.get("sequence"):
                    observed = sha256_text(row["sequence"])
                    if observed != sequence_sha256:
                        raise ValueError(f"Sequence digest mismatch for {protein_id}")
                if protein_id in targets:
                    raise ValueError(f"Protein occurs in multiple target tables: {protein_id}")
                targets[protein_id] = sequence_sha256
    return dict(sorted(targets.items()))


def benchmark_files(benchmark_dir: Path) -> List[dict]:
    records = []
    for name in EXPECTED_CSVS:
        path = benchmark_dir / name
        if not path.is_file():
            raise ValueError(f"Missing benchmark CSV: {path}")
        records.append(
            {
                "name": name,
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        )
    return records


def parse_labelled_paths(values: Sequence[str]) -> List[dict]:
    records = []
    seen: Set[str] = set()
    for value in values:
        if "=" not in value:
            raise ValueError(f"Expected LABEL=PATH, got {value!r}")
        label, raw_path = value.split("=", 1)
        path = Path(raw_path)
        if not label or label in seen:
            raise ValueError(f"Invalid or repeated source label: {label!r}")
        if not path.is_file():
            raise ValueError(f"Missing provenance source {label}: {path}")
        seen.add(label)
        records.append(
            {
                "label": label,
                "name": path.name,
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        )
    return sorted(records, key=lambda record: record["label"])


def parse_runtime_values(values: Sequence[str]) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Expected KEY=VALUE, got {value!r}")
        key, item = value.split("=", 1)
        if not key or key in result:
            raise ValueError(f"Invalid or repeated runtime key: {key!r}")
        result[key] = item
    return dict(sorted(result.items()))


def provenance_file(path: Path, include_path: bool = False) -> dict:
    if not path.is_file():
        raise ValueError(f"Missing provenance file: {path}")
    record = {
        "name": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if include_path:
        record["path"] = str(path.resolve())
    return record


def targets_tsv(targets: Mapping[str, str]) -> str:
    lines = ["protein_id\tsequence_sha256"]
    lines.extend(f"{protein_id}\t{digest}" for protein_id, digest in targets.items())
    return "\n".join(lines) + "\n"


def load_target_manifest(state_root: Path) -> Dict[str, str]:
    path = state_root / "targets.tsv"
    if not path.is_file():
        raise ValueError(f"Missing state target manifest: {path}")
    result: Dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            protein_id = row["protein_id"]
            if protein_id in result:
                raise ValueError(f"Duplicate state target: {protein_id}")
            result[protein_id] = row["sequence_sha256"]
    return result


def build_contract(args: argparse.Namespace, policy: dict, targets: Mapping[str, str]) -> dict:
    environment = None
    if args.environment_report:
        environment_path = Path(args.environment_report)
        if not environment_path.is_file():
            raise ValueError(f"Missing environment report: {environment_path}")
        environment = {
            "name": environment_path.name,
            "sha256": sha256_file(environment_path),
            "size_bytes": environment_path.stat().st_size,
        }
    target_content = targets_tsv(targets)
    baseline = None
    if args.baseline_archive or args.baseline_assembly_report:
        if not args.baseline_archive or not args.baseline_assembly_report:
            raise ValueError(
                "--baseline-archive and --baseline-assembly-report must be used together"
            )
        baseline = {
            "archive": provenance_file(Path(args.baseline_archive), include_path=True),
            "assembly_report": provenance_file(
                Path(args.baseline_assembly_report), include_path=True
            ),
        }
    contract = {
        "schema_version": 1,
        "benchmark_id": args.benchmark_id,
        "benchmark_csvs": benchmark_files(Path(args.benchmark_dir)),
        "targets": {
            "count": len(targets),
            "manifest_sha256": sha256_text(target_content),
        },
        "pfp_commit": args.pfp_commit,
        "framework_commit": args.framework_commit,
        "policy": policy,
        "policy_sha256": sha256_file(Path(args.policy)),
        "environment": environment,
        "source_files": parse_labelled_paths(args.source_file),
        "runtime": parse_runtime_values(args.runtime_value),
    }
    if baseline is not None:
        contract["baseline"] = baseline
    contract["contract_sha256"] = sha256_text(canonical_json(contract))
    return contract


def contract_comparison_payload(
    contract: Mapping[str, object], *, strict_framework_commit: bool
) -> dict:
    """Return contract content with non-comparable metadata removed."""
    payload = dict(contract)
    payload.pop("contract_sha256", None)
    if not strict_framework_commit:
        payload.pop("framework_commit", None)
    return payload


def contract_difference_paths(
    expected: object, observed: object, path: str = "$"
) -> List[str]:
    """Describe structural differences without dumping large contract values."""
    if isinstance(expected, Mapping) and isinstance(observed, Mapping):
        differences: List[str] = []
        for key in sorted(set(expected) | set(observed)):
            child_path = f"{path}.{key}"
            if key not in expected:
                differences.append(f"{child_path}:unexpected")
            elif key not in observed:
                differences.append(f"{child_path}:missing")
            else:
                differences.extend(
                    contract_difference_paths(expected[key], observed[key], child_path)
                )
        return differences
    if isinstance(expected, list) and isinstance(observed, list):
        differences = []
        if len(expected) != len(observed):
            differences.append(f"{path}.length:{len(expected)}!={len(observed)}")
        for index, (expected_item, observed_item) in enumerate(zip(expected, observed)):
            differences.extend(
                contract_difference_paths(
                    expected_item, observed_item, f"{path}[{index}]"
                )
            )
        return differences
    if expected != observed:
        return [f"{path}:value"]
    return []


def load_contract(state_root: Path) -> dict:
    path = state_root / "contract.json"
    if not path.is_file():
        raise ValueError(f"State is not initialized: {state_root}")
    return json.loads(path.read_text(encoding="utf-8"))


def baseline_archive_path(contract: Mapping[str, object]) -> Optional[Path]:
    baseline = contract.get("baseline")
    if not baseline:
        return None
    return Path(baseline["archive"]["path"])  # type: ignore[index]


def baseline_assembly_report_path(
    contract: Mapping[str, object],
) -> Optional[Path]:
    baseline = contract.get("baseline")
    if not baseline:
        return None
    return Path(baseline["assembly_report"]["path"])  # type: ignore[index]


def build_baseline_index(
    report_path: Path,
    targets: Mapping[str, str],
    policies: Mapping[str, dict],
) -> Tuple[str, Set[str], dict]:
    seen: Dict[str, int] = {}
    modality_bits = {
        modality: 1 << index for index, modality in enumerate(sorted(policies))
    }
    expected_mask = sum(modality_bits.values())
    available_counts = {modality: 0 for modality in policies}
    archive_members: Set[str] = set()
    lines = ["protein_id\tmodality\tarchive_member"]
    row_count = 0

    opener = gzip.open if report_path.suffix == ".gz" else open
    with opener(report_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "modality", "status", "dimension"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Baseline assembly report must contain {sorted(required)}: {report_path}"
            )
        for row in reader:
            row_count += 1
            protein_id = row["protein_id"]
            report_modality = row["modality"]
            if protein_id not in targets:
                raise ValueError(
                    f"Baseline assembly report contains unknown protein: {protein_id}"
                )
            if report_modality not in ASSEMBLY_MODALITIES:
                raise ValueError(
                    f"Unknown baseline assembly modality: {report_modality}"
                )
            modality = ASSEMBLY_MODALITIES[report_modality]
            specification = policies[modality]
            bit = modality_bits[modality]
            if seen.get(protein_id, 0) & bit:
                raise ValueError(
                    f"Repeated baseline pair: {(protein_id, modality)}"
                )
            seen[protein_id] = seen.get(protein_id, 0) | bit
            if int(row["dimension"]) != int(specification["dimension"]):
                raise ValueError(
                    f"Baseline dimension mismatch for {(protein_id, modality)}"
                )
            status = row["status"]
            if status not in {"available", "missing"}:
                raise ValueError(
                    f"Unknown baseline status for {(protein_id, modality)}: {status}"
                )
            if status == "available":
                member = (
                    "data/embedding_cache/"
                    f"{specification['cache_directory']}/{protein_id}.npy"
                )
                if member in archive_members:
                    raise ValueError(f"Repeated baseline archive member: {member}")
                archive_members.add(member)
                lines.append(f"{protein_id}\t{modality}\t{member}")
                available_counts[modality] += 1

    expected_rows = len(targets) * len(policies)
    if row_count != expected_rows:
        raise ValueError(
            f"Baseline assembly report has {row_count} rows, expected {expected_rows}"
        )
    incomplete = [
        protein_id
        for protein_id in targets
        if seen.get(protein_id, 0) != expected_mask
    ]
    if incomplete:
        raise ValueError(
            "Baseline assembly report does not contain all modalities for "
            f"{len(incomplete)} proteins; sample={incomplete[:5]}"
        )
    summary = {
        "report_rows": row_count,
        "available_pairs": len(archive_members),
        "available_by_modality": available_counts,
    }
    return "\n".join(lines) + "\n", archive_members, summary


def verify_baseline_archive_members(
    archive_path: Path, expected: Set[str]
) -> Tuple[dict, Dict[str, str]]:
    missing = set(expected)
    observed_npy = 0
    member_hashes: Dict[str, str] = {}
    with tarfile.open(archive_path, mode="r|gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".npy"):
                continue
            observed_npy += 1
            if member.name not in expected:
                raise ValueError(
                    f"Baseline archive contains an unreported array: {member.name}"
                )
            source = archive.extractfile(member)
            if source is None:
                raise ValueError(f"Cannot read baseline archive member: {member.name}")
            digest = hashlib.sha256()
            with source:
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    digest.update(chunk)
            member_hashes[member.name] = digest.hexdigest()
            missing.discard(member.name)
    if missing:
        sample = sorted(missing)[:5]
        raise ValueError(
            f"Baseline archive is missing {len(missing)} reported arrays; sample={sample}"
        )
    if observed_npy != len(expected):
        raise ValueError(
            "Baseline archive contains repeated array members: "
            f"observed={observed_npy};expected={len(expected)}"
        )
    return (
        {"observed_arrays": observed_npy, "expected_arrays": len(expected)},
        member_hashes,
    )


def bind_baseline_index_hashes(index_text: str, member_hashes: Mapping[str, str]) -> str:
    rows = list(csv.DictReader(index_text.splitlines(), delimiter="\t"))
    lines = ["protein_id\tmodality\tarchive_member\tembedding_sha256"]
    for row in rows:
        member = row["archive_member"]
        digest = member_hashes.get(member)
        if digest is None:
            raise ValueError(f"Baseline archive hash is missing for: {member}")
        lines.append(
            f"{row['protein_id']}\t{row['modality']}\t{member}\t{digest}"
        )
    return "\n".join(lines) + "\n"


def load_baseline_accepted(state_root: Path) -> Dict[str, Set[str]]:
    result: Dict[str, Set[str]] = {}
    path = state_root / "baseline_accepted.tsv"
    if not path.is_file():
        return result
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result.setdefault(row["modality"], set()).add(row["protein_id"])
    return result


def load_baseline_members(state_root: Path) -> Dict[Tuple[str, str], str]:
    result: Dict[Tuple[str, str], str] = {}
    path = state_root / "baseline_accepted.tsv"
    if not path.is_file():
        return result
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[(row["protein_id"], row["modality"])] = row["archive_member"]
    return result


def load_baseline_hashes(state_root: Path) -> Dict[Tuple[str, str], str]:
    result: Dict[Tuple[str, str], str] = {}
    path = state_root / "baseline_accepted.tsv"
    if not path.is_file():
        return result
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "embedding_sha256" not in reader.fieldnames:
            raise ValueError(
                "Baseline index lacks embedding hashes; run upgrade-evidence-hashes "
                "after all retry jobs have finished"
            )
        for row in reader:
            digest = row["embedding_sha256"]
            if len(digest) != 64 or any(
                character not in "0123456789abcdef" for character in digest
            ):
                raise ValueError(
                    "Invalid baseline embedding SHA-256 for "
                    f"{(row['protein_id'], row['modality'])}"
                )
            result[(row["protein_id"], row["modality"])] = digest
    return result


def modality_policy(contract: Mapping[str, object]) -> Mapping[str, dict]:
    return contract["policy"]["modalities"]  # type: ignore[index]


def validate_array(path: Path, expected_dimension: int) -> Tuple[np.ndarray, str]:
    try:
        array = np.load(path, allow_pickle=False)
    except Exception as error:
        raise ValueError(f"cannot_load:{type(error).__name__}:{error}") from error
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError(f"non_numeric_dtype:{array.dtype}")
    array = np.squeeze(array)
    if array.ndim != 1 or array.shape[0] != expected_dimension:
        raise ValueError(f"wrong_shape:{tuple(array.shape)}:expected=({expected_dimension},)")
    if not np.isfinite(array).all():
        raise ValueError("non_finite_values")
    return array, sha256_file(path)


def cache_path(state_root: Path, specification: Mapping[str, object], protein_id: str) -> Path:
    return state_root / "cache" / str(specification["cache_directory"]) / f"{protein_id}.npy"


def load_failure_ledger(state_root: Path) -> Dict[Tuple[str, str], dict]:
    path = state_root / "failure_ledger.tsv"
    if not path.is_file():
        return {}
    result: Dict[Tuple[str, str], dict] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["protein_id"], row["modality"])
            result[key] = row
    return result


def load_pairs(path: Optional[Path], targets: Mapping[str, str], modalities: Iterable[str]) -> Set[Tuple[str, str]]:
    if path is None:
        return {(protein_id, modality) for protein_id in targets for modality in modalities}
    result: Set[Tuple[str, str]] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["protein_id"], row["modality"])
            if key[0] not in targets:
                raise ValueError(f"Requested pair contains unknown protein: {key[0]}")
            if key[1] not in modalities:
                raise ValueError(f"Requested pair contains unknown modality: {key[1]}")
            if key in result:
                raise ValueError(f"Repeated requested pair: {key}")
            result.add(key)
    return result


def modality_exit_statuses(path: Optional[Path]) -> Dict[str, int]:
    if path is None or not path.is_file():
        return {}
    result: Dict[str, int] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["modality"]] = int(row["exit_status"])
    return result


def alphafold_reasons(path: Optional[Path]) -> Dict[str, Tuple[str, str]]:
    if path is None or not path.is_file():
        return {}
    reasons: Dict[str, Tuple[str, str]] = {}
    section = ""
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if line in {
            "FOUND IN ALPHAFOLD",
            "NOT FOUND IN ALPHAFOLD",
            "NO UNIPROT MAPPING",
            "API ERRORS",
        }:
            section = line
            continue
        if not line or line.startswith("=") or line.startswith("CAFA_ID"):
            continue
        fields = raw_line.split("\t")
        protein_id = fields[0].strip()
        if not safe_protein_id(protein_id):
            continue
        if section == "NOT FOUND IN ALPHAFOLD":
            reasons[protein_id] = ("alphafold_not_found", "HTTP 404 or empty AlphaFold result")
        elif section == "NO UNIPROT MAPPING":
            reasons[protein_id] = ("alphafold_no_mapping", "No UniProt mapping")
        elif section == "API ERRORS":
            detail = fields[-1].strip() if len(fields) > 1 else "AlphaFold API error"
            reasons[protein_id] = ("alphafold_api_error", detail)
    return reasons


def alphafold_prefetch_reasons(path: Optional[Path]) -> Dict[str, Tuple[str, str]]:
    if path is None or not path.is_file():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {
        str(protein_id): ("alphafold_pdb_download_failed", str(detail))
        for protein_id, detail in payload.get("download_failures", {}).items()
    }


def accepted_ids(
    state_root: Path,
    targets: Mapping[str, str],
    specification: Mapping[str, object],
    modality: Optional[str] = None,
) -> Set[str]:
    directory = state_root / "cache" / str(specification["cache_directory"])
    result: Set[str] = set()
    if directory.is_dir():
        for path in directory.glob("*.npy"):
            if path.stem not in targets:
                raise ValueError(
                    f"Persistent cache contains an out-of-contract protein: {path}"
                )
            result.add(path.stem)
    if modality is not None:
        baseline = load_baseline_accepted(state_root).get(modality, set())
        unknown = baseline - set(targets)
        if unknown:
            raise ValueError(
                "Baseline index contains out-of-contract proteins: "
                f"{sorted(unknown)[:5]}"
            )
        result.update(baseline)
    return result


def acceptance_threshold(specification: Mapping[str, object], target_count: int) -> int:
    if "min_accepted_count" in specification:
        value = int(str(specification["min_accepted_count"]))
    else:
        value = math.ceil(
            float(str(specification["min_accepted_fraction"])) * target_count
        )
    if value < 0 or value > target_count:
        raise ValueError(f"Acceptance threshold {value} is outside 0..{target_count}")
    return value


def refresh_reports(state_root: Path) -> dict:
    evidence_marker = state_root / "EVIDENCE_HASHES_COMPLETE.json"
    if evidence_marker.exists():
        evidence_marker.unlink()
    contract = load_contract(state_root)
    targets = load_target_manifest(state_root)
    policies = modality_policy(contract)
    failures = load_failure_ledger(state_root)
    coverage: Dict[str, dict] = {}
    accepted_by_modality: Dict[str, Set[str]] = {}
    baseline_hashes = load_baseline_hashes(state_root)
    for modality, specification in policies.items():
        accepted_set = accepted_ids(state_root, targets, specification, modality)
        accepted_by_modality[modality] = accepted_set
        threshold = acceptance_threshold(specification, len(targets))
        coverage[modality] = {
            "accepted": len(accepted_set),
            "needs_retry": len(targets) - len(accepted_set),
            "target_count": len(targets),
            "fraction": len(accepted_set) / len(targets) if targets else 0.0,
            "required_accepted": threshold,
            "gate_passed": len(accepted_set) >= threshold,
        }

    status_lines = [
        "protein_id\tmodality\tstate\tsequence_sha256\tembedding_sha256\t"
        "attempts\tlatest_reason\tlatest_detail"
    ]
    retry_lines = [
        "protein_id\tmodality\tsequence_sha256\tattempts\tlatest_reason\tlatest_detail"
    ]
    for protein_id, sequence_sha256 in targets.items():
        for modality in sorted(policies):
            key = (protein_id, modality)
            pair_accepted = protein_id in accepted_by_modality[modality]
            failure = failures.get(key, {})
            attempts = failure.get("attempts", "0")
            reason = failure.get("latest_reason", "") if not pair_accepted else ""
            detail = failure.get("latest_detail", "") if not pair_accepted else ""
            if not pair_accepted and not reason:
                reason = "not_attempted"
            embedding_sha256 = ""
            if pair_accepted:
                embedding_sha256 = baseline_hashes.get(key, "")
                if not embedding_sha256:
                    destination = cache_path(state_root, policies[modality], protein_id)
                    _, embedding_sha256 = validate_array(
                        destination, int(policies[modality]["dimension"])
                    )
            values = [
                protein_id,
                modality,
                "accepted" if pair_accepted else "needs_retry",
                sequence_sha256,
                embedding_sha256,
                str(attempts),
                str(reason).replace("\t", " ").replace("\n", " "),
                str(detail).replace("\t", " ").replace("\n", " "),
            ]
            status_lines.append("\t".join(values))
            if not pair_accepted:
                retry_lines.append("\t".join(values[:2] + values[3:4] + values[5:]))

    overall_passed = all(item["gate_passed"] for item in coverage.values())
    summary = {
        "schema_version": 1,
        "state_root": str(state_root.resolve()),
        "contract_sha256": contract["contract_sha256"],
        "refreshed_at": utc_now(),
        "target_count": len(targets),
        "coverage": coverage,
        "embedding_gate_passed": overall_passed,
    }
    atomic_write_text(state_root / "pair_status.tsv", "\n".join(status_lines) + "\n")
    atomic_write_text(state_root / "needs_retry.tsv", "\n".join(retry_lines) + "\n")
    atomic_write_json(state_root / "coverage.json", summary)
    passed_marker = state_root / "EMBEDDING_GATE_PASSED.json"
    incomplete_marker = state_root / "GENERATION_INCOMPLETE.json"
    if overall_passed:
        atomic_write_json(passed_marker, summary)
        if incomplete_marker.exists():
            incomplete_marker.unlink()
    else:
        atomic_write_json(incomplete_marker, summary)
        if passed_marker.exists():
            passed_marker.unlink()
    return summary


def command_initialize(args: argparse.Namespace) -> dict:
    state_root = Path(args.state_root)
    policy = load_policy(Path(args.policy))
    if args.data_dir:
        targets = load_targets(Path(args.data_dir))
    else:
        targets = load_target_tables([Path(path) for path in args.target_table])
    contract = build_contract(args, policy, targets)
    baseline_index = None
    baseline_validation = None
    if args.baseline_archive:
        baseline_index, archive_members, report_summary = build_baseline_index(
            Path(args.baseline_assembly_report), targets, policy["modalities"]
        )
        archive_summary, member_hashes = verify_baseline_archive_members(
            Path(args.baseline_archive), archive_members
        )
        baseline_index = bind_baseline_index_hashes(baseline_index, member_hashes)
        baseline_validation = {
            "schema_version": 1,
            "assembly_report": report_summary,
            "archive": archive_summary,
        }
    with state_lock(state_root):
        contract_path = state_root / "contract.json"
        targets_path = state_root / "targets.tsv"
        if contract_path.exists() or targets_path.exists():
            existing = load_contract(state_root)
            existing_payload = contract_comparison_payload(
                existing, strict_framework_commit=args.strict_framework_commit
            )
            requested_payload = contract_comparison_payload(
                contract, strict_framework_commit=args.strict_framework_commit
            )
            if existing_payload != requested_payload:
                differences = contract_difference_paths(existing_payload, requested_payload)
                raise ValueError(
                    "Persistent embedding state contract mismatch; use a new state root. "
                    f"existing={existing.get('contract_sha256')} "
                    f"requested={contract.get('contract_sha256')} "
                    f"differing_fields={','.join(differences[:20])}"
                )
            if existing != contract:
                print(
                    "WARNING: current framework commit differs from the state "
                    "initializer; continuing because framework commit matching is "
                    "disabled and every scientific contract field matched exactly",
                    file=sys.stderr,
                )
                contract = existing
            existing_targets = load_target_manifest(state_root)
            if existing_targets != targets:
                raise ValueError("Persistent target manifest changed without a contract change")
        else:
            atomic_write_json(contract_path, contract)
            atomic_write_text(targets_path, targets_tsv(targets))
            for specification in policy["modalities"].values():
                (state_root / "cache" / specification["cache_directory"]).mkdir(
                    parents=True, exist_ok=True
                )
        if baseline_index is not None:
            baseline_path = state_root / "baseline_accepted.tsv"
            if baseline_path.is_file():
                existing_index = baseline_path.read_text(encoding="utf-8")
                if existing_index != baseline_index:
                    existing_rows = list(
                        csv.DictReader(existing_index.splitlines(), delimiter="\t")
                    )
                    replacement_rows = list(
                        csv.DictReader(baseline_index.splitlines(), delimiter="\t")
                    )
                    legacy_existing = [
                        (row["protein_id"], row["modality"], row["archive_member"])
                        for row in existing_rows
                    ]
                    legacy_replacement = [
                        (row["protein_id"], row["modality"], row["archive_member"])
                        for row in replacement_rows
                    ]
                    if legacy_existing != legacy_replacement:
                        raise ValueError("Persistent baseline index changed")
                    atomic_write_text(baseline_path, baseline_index)
            else:
                atomic_write_text(baseline_path, baseline_index)
            atomic_write_json(
                state_root / "baseline_validation.json", baseline_validation
            )
        return refresh_reports(state_root)


def failure_rows_tsv(rows: Mapping[Tuple[str, str], Mapping[str, object]]) -> str:
    columns = [
        "protein_id",
        "modality",
        "attempts",
        "first_failed_at",
        "last_failed_at",
        "latest_reason",
        "latest_detail",
        "latest_attempt_id",
    ]
    lines = ["\t".join(columns)]
    for key in sorted(rows):
        row = rows[key]
        lines.append(
            "\t".join(
                str(row.get(column, "")).replace("\t", " ").replace("\n", " ")
                for column in columns
            )
        )
    return "\n".join(lines) + "\n"


def copy_small_reports(state_root: Path, report_dir: Optional[Path], summary: dict, merge: Optional[dict] = None) -> None:
    if report_dir is None:
        return
    report_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "contract.json",
        "targets.tsv",
        "baseline_validation.json",
        "coverage.json",
        "needs_retry.tsv",
        "failure_ledger.tsv",
        "EMBEDDING_GATE_PASSED.json",
        "GENERATION_INCOMPLETE.json",
    ):
        source = state_root / name
        if source.is_file():
            shutil.copy2(source, report_dir / name)
    atomic_write_json(report_dir / "state_summary.json", summary)
    if merge is not None:
        atomic_write_json(report_dir / "merge_summary.json", merge)


def command_merge(args: argparse.Namespace) -> dict:
    state_root = Path(args.state_root)
    generated_root = Path(args.generated_cache_root)
    with state_lock(state_root):
        contract = load_contract(state_root)
        targets = load_target_manifest(state_root)
        policies = modality_policy(contract)
        requested = load_pairs(
            Path(args.requested_pairs) if args.requested_pairs else None,
            targets,
            policies,
        )
        allowed_extra = load_pairs(
            Path(args.allowed_extra_pairs) if args.allowed_extra_pairs else Path(args.requested_pairs) if args.requested_pairs else None,
            targets,
            policies,
        )
        if not args.allowed_extra_pairs:
            allowed_extra = set(requested)
        statuses = modality_exit_statuses(
            Path(args.modality_status) if args.modality_status else None
        )
        structure_reasons = alphafold_reasons(
            Path(args.alphafold_report) if args.alphafold_report else None
        )
        structure_reasons.update(
            alphafold_prefetch_reasons(
                Path(args.alphafold_prefetch_report)
                if args.alphafold_prefetch_report
                else None
            )
        )
        failures = load_failure_ledger(state_root)
        accepted_before = {
            modality: accepted_ids(state_root, targets, specification, modality)
            for modality, specification in policies.items()
        }
        timestamp = utc_now()
        accepted_count = 0
        already_accepted = 0
        failed_count = 0
        invalid_count = 0
        generated_seen: Set[Tuple[str, str]] = set()

        for modality, specification in policies.items():
            directory = generated_root / specification["cache_directory"]
            if not directory.is_dir():
                continue
            for path in directory.glob("*.npy"):
                key = (path.stem, modality)
                if path.stem not in targets:
                    raise ValueError(f"Generated cache contains unknown target: {path}")
                if key not in requested and key not in allowed_extra:
                    raise ValueError(f"Generated cache contains unrequested pair: {key}")
                generated_seen.add(key)

        for protein_id, modality in sorted(requested):
            specification = policies[modality]
            destination = cache_path(state_root, specification, protein_id)
            key = (protein_id, modality)
            if protein_id in accepted_before[modality]:
                if destination.is_file():
                    validate_array(destination, int(specification["dimension"]))
                already_accepted += 1
                failures.pop(key, None)
                continue

            source = generated_root / specification["cache_directory"] / f"{protein_id}.npy"
            reason = ""
            detail = ""
            if source.is_file():
                try:
                    _, source_sha = validate_array(source, int(specification["dimension"]))
                    atomic_copy(source, destination)
                    _, destination_sha = validate_array(
                        destination, int(specification["dimension"])
                    )
                    if destination_sha != source_sha:
                        raise ValueError("published_sha256_mismatch")
                    accepted_count += 1
                    failures.pop(key, None)
                    continue
                except ValueError as error:
                    invalid_count += 1
                    reason = "invalid_generated_array"
                    detail = str(error)
            else:
                if modality == "structure" and protein_id in structure_reasons:
                    reason, detail = structure_reasons[protein_id]
                elif statuses.get(modality, 0) != 0:
                    reason = "generator_exit_status"
                    detail = str(statuses[modality])
                else:
                    reason = "missing_after_generation"
                    detail = "No validated array was produced"

            failed_count += 1
            previous = failures.get(key)
            failures[key] = {
                "protein_id": protein_id,
                "modality": modality,
                "attempts": int(previous.get("attempts", 0)) + 1 if previous else 1,
                "first_failed_at": previous.get("first_failed_at", timestamp) if previous else timestamp,
                "last_failed_at": timestamp,
                "latest_reason": reason,
                "latest_detail": detail,
                "latest_attempt_id": args.attempt_id,
            }

        atomic_write_text(state_root / "failure_ledger.tsv", failure_rows_tsv(failures))
        summary = refresh_reports(state_root)
        merge_summary = {
            "schema_version": 1,
            "attempt_id": args.attempt_id,
            "merged_at": timestamp,
            "requested_pairs": len(requested),
            "generated_files_seen": len(generated_seen),
            "newly_accepted": accepted_count,
            "already_accepted": already_accepted,
            "failed": failed_count,
            "invalid": invalid_count,
            "embedding_gate_passed": summary["embedding_gate_passed"],
            "coverage": summary["coverage"],
        }
        atomic_write_json(state_root / "last_merge.json", merge_summary)
        copy_small_reports(
            state_root,
            Path(args.report_dir) if args.report_dir else None,
            summary,
            merge_summary,
        )
        return merge_summary


def command_pending(args: argparse.Namespace) -> dict:
    state_root = Path(args.state_root)
    with state_lock(state_root):
        summary = refresh_reports(state_root)
        source = state_root / "needs_retry.tsv"
        rows = []
        with source.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if args.modality and row["modality"] != args.modality:
                    continue
                rows.append(row)
        columns = [
            "protein_id",
            "modality",
            "sequence_sha256",
            "attempts",
            "latest_reason",
            "latest_detail",
        ]
        lines = ["\t".join(columns)]
        lines.extend("\t".join(row[column] for column in columns) for row in rows)
        atomic_write_text(Path(args.output), "\n".join(lines) + "\n")
        return {
            "output": str(Path(args.output).resolve()),
            "modality": args.modality,
            "pair_count": len(rows),
            "embedding_gate_passed": summary["embedding_gate_passed"],
        }


def load_planner_global_splits(plan_dir: Path) -> Dict[str, str]:
    split_by_id: Dict[str, str] = {}
    for action in ("reuse", "regenerate"):
        path = plan_dir / f"{action}_proteins.tsv"
        if not path.is_file():
            raise ValueError(f"Missing planner action table: {path}")
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"protein_id", "target_memberships"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(f"Invalid planner action-table schema: {path}")
            for row in reader:
                protein_id = row["protein_id"]
                memberships = json.loads(row["target_memberships"])
                splits = {
                    name.removesuffix(".csv").rsplit("-", 1)[1]
                    for name in memberships
                }
                if len(splits) != 1:
                    raise ValueError(
                        f"Protein {protein_id} does not have one global target split: "
                        f"{memberships}"
                    )
                split = next(iter(splits))
                if split not in {"training", "validation", "test"}:
                    raise ValueError(f"Unknown target split for {protein_id}: {split}")
                if protein_id in split_by_id:
                    raise ValueError(f"Protein occurs in both planner actions: {protein_id}")
                split_by_id[protein_id] = split
    return split_by_id


def balanced_controls(
    accepted: Sequence[str], split_by_id: Mapping[str, str], count: int
) -> Tuple[List[str], Dict[str, int]]:
    split_order = ("training", "validation", "test")
    buckets = {
        split: [
            protein_id
            for protein_id in accepted
            if split_by_id.get(protein_id) == split
        ]
        for split in split_order
    }
    missing_plan_ids = sorted(set(accepted) - set(split_by_id))
    if missing_plan_ids:
        raise ValueError(
            f"Accepted controls are absent from planner tables: {missing_plan_ids[:10]} "
            f"(total={len(missing_plan_ids)})"
        )
    empty_splits = [split for split in split_order if not buckets[split]]
    if empty_splits:
        raise ValueError(
            "Cannot balance diagnostic controls; no accepted proteins for: "
            + ", ".join(empty_splits)
        )
    if count < len(split_order):
        raise ValueError(
            f"Balanced controls require at least {len(split_order)} proteins"
        )

    selected: List[str] = []
    offsets = {split: 0 for split in split_order}
    while len(selected) < count:
        progressed = False
        for split in split_order:
            offset = offsets[split]
            if offset >= len(buckets[split]):
                continue
            selected.append(buckets[split][offset])
            offsets[split] += 1
            progressed = True
            if len(selected) == count:
                break
        if not progressed:
            break
    counts = {
        split: sum(split_by_id[item] == split for item in selected)
        for split in split_order
    }
    return selected, counts


def command_controls(args: argparse.Namespace) -> dict:
    if args.count <= 0:
        raise ValueError("--count must be positive")
    state_root = Path(args.state_root)
    with state_lock(state_root):
        contract = load_contract(state_root)
        targets = load_target_manifest(state_root)
        policies = modality_policy(contract)
        if args.modality not in policies:
            raise ValueError(f"Unknown modality: {args.modality}")
        accepted = sorted(
            accepted_ids(
                state_root,
                targets,
                policies[args.modality],
                args.modality,
            )
        )
        split_counts = None
        if args.balance_global_splits:
            if args.plan_dir is None:
                raise ValueError("--balance-global-splits requires --plan-dir")
            selected, split_counts = balanced_controls(
                accepted,
                load_planner_global_splits(Path(args.plan_dir)),
                args.count,
            )
        else:
            selected = accepted[: args.count]
        lines = ["protein_id\tmodality\tsequence_sha256"]
        lines.extend(
            f"{protein_id}\t{args.modality}\t{targets[protein_id]}"
            for protein_id in selected
        )
        atomic_write_text(Path(args.output), "\n".join(lines) + "\n")
        result = {
            "output": str(Path(args.output).resolve()),
            "modality": args.modality,
            "requested_count": args.count,
            "selected_count": len(selected),
        }
        if split_counts is not None:
            result["global_split_counts"] = split_counts
        return result


def materialize_pairs(
    state_root: Path,
    output_root: Path,
    pairs: Set[Tuple[str, str]],
) -> dict:
    contract = load_contract(state_root)
    targets = load_target_manifest(state_root)
    policies = modality_policy(contract)
    baseline_members = load_baseline_members(state_root)
    copied = 0
    skipped = 0
    baseline_copied = 0
    delta_copied = 0
    baseline_requested: Dict[str, Tuple[str, str]] = {}

    for protein_id, modality in pairs:
        if protein_id not in targets or modality not in policies:
            raise ValueError(f"Unknown materialization pair: {(protein_id, modality)}")
        member = baseline_members.get((protein_id, modality))
        source = cache_path(state_root, policies[modality], protein_id)
        if member and source.is_file():
            raise ValueError(
                f"Pair exists in both baseline and retry delta: {(protein_id, modality)}"
            )
        if member:
            baseline_requested[member] = (protein_id, modality)
        elif not source.is_file():
            raise ValueError(f"Pair is not accepted: {(protein_id, modality)}")

    if baseline_requested:
        archive_path = baseline_archive_path(contract)
        if archive_path is None or not archive_path.is_file():
            raise ValueError(f"Baseline archive is unavailable: {archive_path}")
        missing_members = set(baseline_requested)
        with tarfile.open(archive_path, mode="r|gz") as archive:
            for archive_member in archive:
                pair = baseline_requested.get(archive_member.name)
                if pair is None:
                    continue
                if not archive_member.isfile():
                    raise ValueError(
                        f"Baseline member is not a regular file: {archive_member.name}"
                    )
                protein_id, modality = pair
                specification = policies[modality]
                destination = (
                    output_root
                    / str(specification["cache_directory"])
                    / f"{protein_id}.npy"
                )
                if destination.is_file():
                    validate_array(destination, int(specification["dimension"]))
                    skipped += 1
                else:
                    source_handle = archive.extractfile(archive_member)
                    if source_handle is None:
                        raise ValueError(
                            f"Cannot read baseline member: {archive_member.name}"
                        )
                    with source_handle:
                        atomic_copy_stream(source_handle, destination)
                    validate_array(destination, int(specification["dimension"]))
                    copied += 1
                    baseline_copied += 1
                missing_members.discard(archive_member.name)
        if missing_members:
            raise ValueError(
                "Baseline archive is missing requested members: "
                f"{sorted(missing_members)[:5]}"
            )

    for protein_id, modality in sorted(pairs):
        if (protein_id, modality) in baseline_members:
            continue
        specification = policies[modality]
        source = cache_path(state_root, specification, protein_id)
        _, source_sha = validate_array(source, int(specification["dimension"]))
        destination = (
            output_root
            / str(specification["cache_directory"])
            / f"{protein_id}.npy"
        )
        if destination.is_file():
            _, destination_sha = validate_array(
                destination, int(specification["dimension"])
            )
            if destination_sha != source_sha:
                raise ValueError(f"Materialization conflict: {destination}")
            skipped += 1
            continue
        atomic_copy(source, destination)
        _, destination_sha = validate_array(
            destination, int(specification["dimension"])
        )
        if destination_sha != source_sha:
            raise ValueError(f"Materialization checksum mismatch: {destination}")
        copied += 1
        delta_copied += 1

    return {
        "output_cache_root": str(output_root.resolve()),
        "requested_pairs": len(pairs),
        "copied": copied,
        "baseline_copied": baseline_copied,
        "delta_copied": delta_copied,
        "already_present": skipped,
    }


def command_materialize(args: argparse.Namespace) -> dict:
    state_root = Path(args.state_root)
    with state_lock(state_root):
        contract = load_contract(state_root)
        targets = load_target_manifest(state_root)
        pairs = load_pairs(
            Path(args.pairs), targets, modality_policy(contract)
        )
        result = materialize_pairs(state_root, Path(args.output_cache_root), pairs)
        if args.report:
            atomic_write_json(Path(args.report), result)
        return result


def command_hydrate(args: argparse.Namespace) -> dict:
    state_root = Path(args.state_root)
    output_root = Path(args.output_cache_root)
    with state_lock(state_root):
        contract = load_contract(state_root)
        targets = load_target_manifest(state_root)
        policies = modality_policy(contract)
        pairs = {
            (protein_id, modality)
            for modality, specification in policies.items()
            for protein_id in accepted_ids(
                state_root, targets, specification, modality
            )
        }
        result = materialize_pairs(state_root, output_root, pairs)
        if args.preserve_evidence:
            coverage_path = state_root / "coverage.json"
            if not coverage_path.is_file():
                raise ValueError(
                    "--preserve-evidence requires an existing coverage.json"
                )
            summary = json.loads(coverage_path.read_text(encoding="utf-8"))
            if summary.get("contract_sha256") != contract.get("contract_sha256"):
                raise ValueError("Existing coverage.json is not bound to the state contract")
        else:
            summary = refresh_reports(state_root)
        result["embedding_gate_passed"] = summary["embedding_gate_passed"]
        result["coverage"] = summary["coverage"]
        if args.report:
            atomic_write_json(Path(args.report), result)
        return result


def command_summary(args: argparse.Namespace) -> dict:
    state_root = Path(args.state_root)
    with state_lock(state_root):
        summary = refresh_reports(state_root)
        copy_small_reports(
            state_root,
            Path(args.report_dir) if args.report_dir else None,
            summary,
        )
        return summary


def accepted_membership(
    state_root: Path,
    targets: Mapping[str, str],
    policies: Mapping[str, dict],
) -> Dict[str, Set[str]]:
    return {
        modality: accepted_ids(state_root, targets, specification, modality)
        for modality, specification in policies.items()
    }


def command_upgrade_evidence_hashes(args: argparse.Namespace) -> dict:
    """Add array hashes to an existing state without changing its contract."""
    state_root = Path(args.state_root)
    with state_lock(state_root):
        contract = load_contract(state_root)
        targets = load_target_manifest(state_root)
        policies = modality_policy(contract)
        membership_before = accepted_membership(state_root, targets, policies)
        counts_before = {
            modality: len(protein_ids)
            for modality, protein_ids in membership_before.items()
        }

        baseline_path = state_root / "baseline_accepted.tsv"
        baseline_archive = baseline_archive_path(contract)
        assembly_report = baseline_assembly_report_path(contract)
        baseline_pairs = 0
        baseline_archive_summary = None
        baseline_index_changed = False
        if baseline_path.is_file():
            if baseline_archive is None or assembly_report is None:
                raise ValueError(
                    "State has a baseline index but its contract has incomplete baseline provenance"
                )
            if not baseline_archive.is_file():
                raise ValueError(f"Missing contracted baseline archive: {baseline_archive}")
            if not assembly_report.is_file():
                raise ValueError(
                    f"Missing contracted baseline assembly report: {assembly_report}"
                )

            baseline_contract = contract["baseline"]
            expected_archive = baseline_contract["archive"]
            observed_size = baseline_archive.stat().st_size
            expected_size = int(expected_archive["size_bytes"])
            if observed_size != expected_size:
                raise ValueError(
                    "Contracted baseline archive size changed: "
                    f"expected={expected_size};observed={observed_size}"
                )
            observed_archive_sha256 = sha256_file(baseline_archive)
            if observed_archive_sha256 != expected_archive["sha256"]:
                raise ValueError("Contracted baseline archive SHA-256 changed")

            expected_report = baseline_contract["assembly_report"]
            if assembly_report.stat().st_size != int(expected_report["size_bytes"]):
                raise ValueError("Contracted baseline assembly report size changed")
            if sha256_file(assembly_report) != expected_report["sha256"]:
                raise ValueError("Contracted baseline assembly report SHA-256 changed")

            expected_index, expected_members, report_summary = build_baseline_index(
                assembly_report, targets, policies
            )

            current_index = baseline_path.read_text(encoding="utf-8")
            rows = list(csv.DictReader(current_index.splitlines(), delimiter="\t"))
            required = {"protein_id", "modality", "archive_member"}
            fieldnames = set(rows[0]) if rows else set()
            if not rows or not required.issubset(fieldnames):
                raise ValueError(
                    "Baseline index must contain protein_id, modality and archive_member"
                )
            current_mapping = [
                (row["protein_id"], row["modality"], row["archive_member"])
                for row in rows
            ]
            expected_rows = list(
                csv.DictReader(expected_index.splitlines(), delimiter="\t")
            )
            expected_mapping = [
                (row["protein_id"], row["modality"], row["archive_member"])
                for row in expected_rows
            ]
            if current_mapping != expected_mapping:
                raise ValueError(
                    "Baseline index mapping differs from the contracted assembly report"
                )
            baseline_pairs = len(rows)
            baseline_archive_summary, member_hashes = verify_baseline_archive_members(
                baseline_archive, expected_members
            )
            baseline_archive_summary["assembly_report"] = report_summary
            upgraded_index = bind_baseline_index_hashes(
                expected_index, member_hashes
            )
            baseline_index_changed = upgraded_index != current_index
            if baseline_index_changed:
                atomic_write_text(baseline_path, upgraded_index)
        elif baseline_archive is not None:
            raise ValueError(
                "State contract has a baseline archive but baseline_accepted.tsv is missing"
            )

        refresh_reports(state_root)
        membership_after = accepted_membership(state_root, targets, policies)
        counts_after = {
            modality: len(protein_ids)
            for modality, protein_ids in membership_after.items()
        }
        if membership_after != membership_before:
            raise ValueError(
                "Evidence upgrade changed accepted membership: "
                f"before={counts_before};after={counts_after}"
            )

        pair_status_path = state_root / "pair_status.tsv"
        result = {
            "schema_version": 1,
            "state_root": str(state_root.resolve()),
            "contract_sha256": contract["contract_sha256"],
            "upgraded_at": utc_now(),
            "accepted_counts_before": counts_before,
            "accepted_counts_after": counts_after,
            "accepted_membership_unchanged": True,
            "baseline_pairs_hashed": baseline_pairs,
            "baseline_index_changed": baseline_index_changed,
            "baseline_archive": baseline_archive_summary,
            "pair_status_sha256": sha256_file(pair_status_path),
            "coverage_sha256": sha256_file(state_root / "coverage.json"),
        }
        atomic_write_json(state_root / "EVIDENCE_HASHES_COMPLETE.json", result)
        if args.report:
            atomic_write_json(Path(args.report), result)
        return result


def add_state_root(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--state-root", type=Path, required=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    initialize = subparsers.add_parser("initialize")
    add_state_root(initialize)
    initialize.add_argument("--benchmark-id", required=True)
    initialize.add_argument("--benchmark-dir", type=Path, required=True)
    target_source = initialize.add_mutually_exclusive_group(required=True)
    target_source.add_argument("--data-dir", type=Path)
    target_source.add_argument(
        "--target-table", type=Path, action="append", default=[]
    )
    initialize.add_argument("--policy", type=Path, required=True)
    initialize.add_argument("--pfp-commit", required=True)
    initialize.add_argument("--framework-commit", required=True)
    initialize.add_argument("--environment-report", type=Path)
    initialize.add_argument("--source-file", action="append", default=[])
    initialize.add_argument("--runtime-value", action="append", default=[])
    initialize.add_argument("--baseline-archive", type=Path)
    initialize.add_argument("--baseline-assembly-report", type=Path)
    framework_commit_mode = initialize.add_mutually_exclusive_group()
    framework_commit_mode.add_argument(
        "--strict-framework-commit",
        action="store_true",
        help="require the current framework commit to equal an existing state contract",
    )
    framework_commit_mode.add_argument(
        "--allow-framework-commit-drift",
        action="store_false",
        dest="strict_framework_commit",
        help=argparse.SUPPRESS,
    )
    initialize.set_defaults(strict_framework_commit=False)

    merge = subparsers.add_parser("merge")
    add_state_root(merge)
    merge.add_argument("--generated-cache-root", type=Path, required=True)
    merge.add_argument("--attempt-id", required=True)
    merge.add_argument("--requested-pairs", type=Path)
    merge.add_argument("--allowed-extra-pairs", type=Path)
    merge.add_argument("--modality-status", type=Path)
    merge.add_argument("--alphafold-report", type=Path)
    merge.add_argument("--alphafold-prefetch-report", type=Path)
    merge.add_argument("--report-dir", type=Path)

    pending = subparsers.add_parser("pending")
    add_state_root(pending)
    pending.add_argument("--modality", choices=("sequence", "text", "structure", "ppi"))
    pending.add_argument("--output", type=Path, required=True)

    controls = subparsers.add_parser("controls")
    add_state_root(controls)
    controls.add_argument("--modality", choices=("sequence", "text", "structure", "ppi"), required=True)
    controls.add_argument("--count", type=int, default=20)
    controls.add_argument("--plan-dir", type=Path)
    controls.add_argument("--balance-global-splits", action="store_true")
    controls.add_argument("--output", type=Path, required=True)

    materialize = subparsers.add_parser("materialize")
    add_state_root(materialize)
    materialize.add_argument("--pairs", type=Path, required=True)
    materialize.add_argument("--output-cache-root", type=Path, required=True)
    materialize.add_argument("--report", type=Path)

    hydrate = subparsers.add_parser("hydrate")
    add_state_root(hydrate)
    hydrate.add_argument("--output-cache-root", type=Path, required=True)
    hydrate.add_argument("--report", type=Path)
    hydrate.add_argument(
        "--preserve-evidence",
        action="store_true",
        help="materialize without refreshing or invalidating existing state evidence",
    )

    summary = subparsers.add_parser("summary")
    add_state_root(summary)
    summary.add_argument("--report-dir", type=Path)

    upgrade_evidence = subparsers.add_parser("upgrade-evidence-hashes")
    add_state_root(upgrade_evidence)
    upgrade_evidence.add_argument("--report", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    commands = {
        "initialize": command_initialize,
        "merge": command_merge,
        "pending": command_pending,
        "controls": command_controls,
        "materialize": command_materialize,
        "hydrate": command_hydrate,
        "summary": command_summary,
        "upgrade-evidence-hashes": command_upgrade_evidence_hashes,
    }
    try:
        result = commands[args.command](args)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
