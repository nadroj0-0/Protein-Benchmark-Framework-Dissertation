#!/usr/bin/env python3
"""Publish a PPI-only delta from an audited unambiguous STRING mapping."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import tarfile
import tempfile
import uuid
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np


POLICY = "widened_unambiguous"
EXPECTED_SOURCES = {
    "Ensembl_HGNC_uniprot_ids",
    "Ensembl_UniProt",
    "Ensembl_flybase_gene_id",
    "Ensembl_gene",
    "UniProt_AC",
    "UniProt_DR_FlyBase",
    "UniProt_ID",
}
SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy-details", type=Path, required=True)
    parser.add_argument("--audit-summary", type=Path, required=True)
    parser.add_argument("--base-pair-status", type=Path, required=True)
    parser.add_argument("--string-h5", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--expected-policy-details-sha256", required=True)
    parser.add_argument("--expected-audit-summary-sha256", required=True)
    parser.add_argument("--expected-base-pair-status-sha256", required=True)
    parser.add_argument("--expected-string-h5-sha256", required=True)
    parser.add_argument("--expected-target-count", type=int, required=True)
    parser.add_argument("--expected-base-count", type=int, required=True)
    parser.add_argument("--expected-delta-count", type=int, required=True)
    parser.add_argument("--expected-final-count", type=int, required=True)
    parser.add_argument("--protein-chunk-size", type=int, default=250_000)
    parser.add_argument("--embedding-batch-size", type=int, default=4_096)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_hash(path: Path, expected: str, label: str) -> str:
    if not path.is_file():
        raise ValueError(f"Missing {label}: {path}")
    if not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError(f"Invalid expected SHA-256 for {label}")
    observed = sha256_file(path)
    if observed != expected:
        raise ValueError(f"{label} SHA-256 mismatch: {observed} != {expected}")
    return observed


def bool_field(value: str, label: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise ValueError(f"Invalid boolean for {label}: {value!r}")


def load_base_status(path: Path) -> tuple[set[str], set[str]]:
    targets: set[str] = set()
    accepted: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "modality", "state"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Base pair-status table has an unexpected schema")
        for row in reader:
            protein_id = row["protein_id"]
            if not SAFE_ID.fullmatch(protein_id):
                raise ValueError(f"Unsafe protein ID in pair-status table: {protein_id!r}")
            targets.add(protein_id)
            if row["modality"] == "ppi" and row["state"] == "accepted":
                accepted.add(protein_id)
    return targets, accepted


def load_policy_details(
    path: Path, targets: set[str]
) -> tuple[dict[str, tuple[str, str]], int]:
    mappings: dict[str, tuple[str, str]] = {}
    seen: set[str] = set()
    ambiguous = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "protein_id",
            "policy",
            "covered",
            "ambiguous",
            "selected_string_id",
            "selected_source",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Policy detail table has an unexpected schema")
        for row in reader:
            if row["policy"] != POLICY or row["protein_id"] not in targets:
                continue
            protein_id = row["protein_id"]
            if protein_id in seen:
                raise ValueError(f"Duplicate {POLICY} row for {protein_id}")
            seen.add(protein_id)
            covered = bool_field(row["covered"], f"{protein_id}.covered")
            is_ambiguous = bool_field(row["ambiguous"], f"{protein_id}.ambiguous")
            if is_ambiguous:
                ambiguous += 1
            if not covered:
                continue
            if is_ambiguous:
                raise ValueError(f"Ambiguous protein was marked covered: {protein_id}")
            string_id = row["selected_string_id"]
            source = row["selected_source"]
            if not string_id or "." not in string_id or not source:
                raise ValueError(f"Incomplete selected mapping for {protein_id}")
            selected_sources = source.split(";")
            if any(
                not selected_source or selected_source not in EXPECTED_SOURCES
                for selected_source in selected_sources
            ):
                raise ValueError(
                    f"Selected source is outside the frozen policy for {protein_id}: "
                    f"{source}"
                )
            mappings[protein_id] = (string_id, source)
    missing = targets - seen
    if missing:
        sample = ", ".join(sorted(missing)[:5])
        raise ValueError(
            f"Policy details omit {len(missing)} contemporary targets; sample: {sample}"
        )
    return mappings, ambiguous


def validate_audit_summary(path: Path, h5_sha256: str) -> dict:
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("schema_name") != "string-alias-policy-coverage-audit":
        raise ValueError("Unexpected alias-policy audit schema")
    policy = summary.get("policies", {}).get(POLICY, {})
    if set(policy.get("source_tokens", [])) != EXPECTED_SOURCES:
        raise ValueError("Audit widened-policy source tokens changed")
    recorded_h5 = summary.get("input_files", {}).get("string_h5", {}).get("sha256")
    if recorded_h5 != h5_sha256:
        raise ValueError(
            f"Audit summary STRING H5 differs: {recorded_h5} != {h5_sha256}"
        )
    return summary


def decode_string_id(value: object) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def extract_vectors(
    h5_path: Path,
    mappings: dict[str, tuple[str, str]],
    output_dir: Path,
    *,
    protein_chunk_size: int,
    embedding_batch_size: int,
) -> dict[str, str]:
    by_string: dict[str, list[str]] = defaultdict(list)
    for protein_id, (string_id, _source) in mappings.items():
        by_string[string_id].append(protein_id)
    by_species: dict[str, set[str]] = defaultdict(set)
    for string_id in by_string:
        by_species[string_id.split(".", 1)[0]].add(string_id)

    output_dir.mkdir(parents=True)
    vector_hashes: dict[str, str] = {}
    found_string_ids: set[str] = set()
    with h5py.File(h5_path, "r") as handle:
        species_root = handle.get("species")
        if species_root is None:
            raise ValueError("STRING H5 has no species group")
        for species_id in sorted(by_species):
            if species_id not in species_root:
                raise ValueError(f"STRING H5 has no species group {species_id}")
            group = species_root[species_id]
            proteins = group["proteins"]
            embeddings = group["embeddings"]
            wanted = by_species[species_id]
            positions: dict[str, int] = {}
            for start in range(0, len(proteins), protein_chunk_size):
                stop = min(len(proteins), start + protein_chunk_size)
                for offset, value in enumerate(proteins[start:stop]):
                    string_id = decode_string_id(value)
                    if string_id not in wanted:
                        continue
                    if string_id in positions:
                        raise ValueError(f"Duplicate STRING ID in H5: {string_id}")
                    positions[string_id] = start + offset
            missing = wanted - positions.keys()
            if missing:
                sample = ", ".join(sorted(missing)[:5])
                raise ValueError(
                    f"Audit-selected STRING IDs are absent from H5; sample: {sample}"
                )

            ordered = sorted((index, string_id) for string_id, index in positions.items())
            for start in range(0, len(ordered), embedding_batch_size):
                batch = ordered[start : start + embedding_batch_size]
                indices = np.asarray([index for index, _ in batch], dtype=np.int64)
                vectors = np.asarray(embeddings[indices])
                if vectors.ndim != 2 or vectors.shape != (len(batch), 512):
                    raise ValueError(
                        f"Unexpected STRING embedding batch shape: {vectors.shape}"
                    )
                if not np.issubdtype(vectors.dtype, np.number) or not np.isfinite(
                    vectors
                ).all():
                    raise ValueError("STRING embedding batch is non-numeric or non-finite")
                for vector, (_index, string_id) in zip(vectors, batch):
                    for protein_id in sorted(by_string[string_id]):
                        destination = output_dir / f"{protein_id}.npy"
                        np.save(destination, vector, allow_pickle=False)
                        vector_hashes[protein_id] = sha256_file(destination)
                    found_string_ids.add(string_id)
    if found_string_ids != set(by_string):
        raise ValueError("Not every selected STRING ID was extracted")
    if set(vector_hashes) != set(mappings):
        raise ValueError("Not every selected protein received one PPI vector")
    return vector_hashes


def write_mapping_report(
    path: Path,
    mappings: dict[str, tuple[str, str]],
    vector_hashes: dict[str, str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    ["protein_id", "string_id", "selected_source", "npy_sha256"]
                )
                for protein_id in sorted(mappings):
                    string_id, source = mappings[protein_id]
                    writer.writerow(
                        [protein_id, string_id, source, vector_hashes[protein_id]]
                    )


def create_deterministic_archive(source_dir: Path, archive: Path) -> None:
    archive.parent.mkdir(parents=True, exist_ok=True)
    with archive.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as tar:
                for source in sorted(source_dir.glob("*.npy")):
                    info = tar.gettarinfo(str(source), arcname=f"ppi/{source.name}")
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mtime = 0
                    info.mode = 0o644
                    with source.open("rb") as handle:
                        tar.addfile(info, handle)


def safe_extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True)
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar:
            if not member.isfile() or not member.name.startswith("ppi/"):
                raise ValueError(f"Unexpected delta archive member: {member.name}")
            name = Path(member.name)
            if len(name.parts) != 2 or not SAFE_ID.fullmatch(name.stem):
                raise ValueError(f"Unsafe delta archive member: {member.name}")
            source = tar.extractfile(member)
            if source is None:
                raise ValueError(f"Could not read archive member: {member.name}")
            output = destination / name
            output.parent.mkdir(parents=True, exist_ok=True)
            with output.open("wb") as handle:
                shutil.copyfileobj(source, handle)


def validate_roundtrip(directory: Path, expected_hashes: dict[str, str]) -> None:
    files = sorted((directory / "ppi").glob("*.npy"))
    if {path.stem for path in files} != set(expected_hashes):
        raise ValueError("Round-trip archive membership differs")
    for path in files:
        if sha256_file(path) != expected_hashes[path.stem]:
            raise ValueError(f"Round-trip hash differs for {path.stem}")
        vector = np.load(path, allow_pickle=False)
        if vector.shape != (512,) or not np.isfinite(vector).all():
            raise ValueError(f"Invalid round-trip PPI array: {path.name}")


def member_content_sha256(vector_hashes: dict[str, str]) -> str:
    digest = hashlib.sha256()
    for protein_id in sorted(vector_hashes):
        digest.update(
            f"{protein_id}\0{vector_hashes[protein_id]}\n".encode("utf-8")
        )
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main() -> int:
    args = parse_args()
    if args.work_dir.exists():
        raise SystemExit(f"Refusing existing work directory: {args.work_dir}")
    if args.output_root.exists():
        raise SystemExit(f"Refusing existing output root: {args.output_root}")
    if args.protein_chunk_size <= 0 or args.embedding_batch_size <= 0:
        raise SystemExit("Chunk and batch sizes must be positive")

    try:
        policy_sha = require_hash(
            args.policy_details, args.expected_policy_details_sha256, "policy details"
        )
        summary_sha = require_hash(
            args.audit_summary, args.expected_audit_summary_sha256, "audit summary"
        )
        pair_status_sha = require_hash(
            args.base_pair_status,
            args.expected_base_pair_status_sha256,
            "base pair status",
        )
        h5_sha = require_hash(
            args.string_h5, args.expected_string_h5_sha256, "STRING H5"
        )
        audit_summary = validate_audit_summary(args.audit_summary, h5_sha)
        targets, accepted = load_base_status(args.base_pair_status)
        if len(targets) != args.expected_target_count:
            raise ValueError(
                f"Target count differs: {len(targets)} != {args.expected_target_count}"
            )
        if len(accepted) != args.expected_base_count:
            raise ValueError(
                f"Accepted PPI count differs: {len(accepted)} != "
                f"{args.expected_base_count}"
            )
        direct, ambiguous_count = load_policy_details(args.policy_details, targets)
        missing_base = accepted - direct.keys()
        if missing_base:
            sample = ", ".join(sorted(missing_base)[:5])
            raise ValueError(
                f"Widened policy loses {len(missing_base)} accepted baseline PPI "
                f"proteins; sample: {sample}"
            )
        if len(direct) != args.expected_final_count:
            raise ValueError(
                f"Widened policy coverage differs: {len(direct)} != "
                f"{args.expected_final_count}"
            )
        delta = {
            protein_id: mapping
            for protein_id, mapping in direct.items()
            if protein_id not in accepted
        }
        if set(delta) & accepted:
            raise ValueError("Delta selection overlaps accepted PPI membership")
        if len(delta) != args.expected_delta_count:
            raise ValueError(
                f"Delta count differs: {len(delta)} != {args.expected_delta_count}"
            )
        if len(accepted | set(delta)) != args.expected_final_count:
            raise ValueError("Final PPI union count differs")

        args.work_dir.mkdir(parents=True)
        generated = args.work_dir / "generated_ppi"
        roundtrip = args.work_dir / "roundtrip"
        vector_hashes = extract_vectors(
            args.string_h5,
            delta,
            generated,
            protein_chunk_size=args.protein_chunk_size,
            embedding_batch_size=args.embedding_batch_size,
        )

        args.output_root.parent.mkdir(parents=True, exist_ok=True)
        staging = args.output_root.parent / (
            f".{args.output_root.name}.staging-{uuid.uuid4().hex}"
        )
        staging.mkdir()
        try:
            archive = staging / "ppi_delta.tar.gz"
            mapping_report = staging / "mapping.tsv.gz"
            create_deterministic_archive(generated, archive)
            write_mapping_report(mapping_report, delta, vector_hashes)
            safe_extract(archive, roundtrip)
            validate_roundtrip(roundtrip, vector_hashes)
            payload = {
                "schema_name": "validated-unambiguous-ppi-delta",
                "schema_version": 1,
                "complete": True,
                "created_at": datetime.now(timezone.utc).isoformat(),
                "policy": POLICY,
                "policy_description": audit_summary["policies"][POLICY][
                    "description"
                ],
                "target_count": len(targets),
                "base_accepted_count": len(accepted),
                "direct_policy_covered_count": len(direct),
                "ambiguous_rejected_count": ambiguous_count,
                "delta_count": len(delta),
                "final_union_count": len(accepted | set(delta)),
                "base_overlap_count": 0,
                "dimension": 512,
                "archive": archive.name,
                "archive_sha256": sha256_file(archive),
                "mapping_report": mapping_report.name,
                "mapping_report_sha256": sha256_file(mapping_report),
                "member_content_sha256": member_content_sha256(vector_hashes),
                "roundtrip_validated": True,
                "inputs": {
                    "policy_details": str(args.policy_details.resolve()),
                    "policy_details_sha256": policy_sha,
                    "audit_summary": str(args.audit_summary.resolve()),
                    "audit_summary_sha256": summary_sha,
                    "base_pair_status": str(args.base_pair_status.resolve()),
                    "base_pair_status_sha256": pair_status_sha,
                    "string_h5": str(args.string_h5.resolve()),
                    "string_h5_sha256": h5_sha,
                },
            }
            atomic_json(staging / "DELTA_COMPLETE.json", payload)
            os.replace(staging, args.output_root)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    raise SystemExit(main())
