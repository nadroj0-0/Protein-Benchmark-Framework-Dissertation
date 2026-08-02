#!/usr/bin/env python3
"""Assemble a validated PFP cache from a pair ledger and modality deltas."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import tarfile
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import BinaryIO

import numpy as np

from resolve_embedding_reuse_sources import (
    EXPECTED_DIMENSIONS,
    MODALITIES,
    canonical_array_sha256,
    normalized_member,
)


MODALITY_DIRECTORIES = {
    "sequence": "prott5",
    "text": "exp_text_embeddings_temporal",
    "structure": "IF1",
    "ppi": "ppi",
}
PAIR_KEY = tuple[str, str]


class AssemblyError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_array(data: bytes, modality: str) -> tuple[str, str, list[int]]:
    try:
        array = np.load(io.BytesIO(data), allow_pickle=False)
    except Exception as error:
        raise AssemblyError(f"Unreadable {modality} NPY: {error}") from error
    expected = (EXPECTED_DIMENSIONS[modality],)
    if array.shape != expected:
        raise AssemblyError(f"Wrong {modality} shape: {array.shape} != {expected}")
    if array.dtype.kind not in "fc" or not np.isfinite(array).all():
        raise AssemblyError(f"Invalid {modality} array dtype or values")
    return canonical_array_sha256(array), str(array.dtype), list(array.shape)


def verify_ledger(ledger_dir: Path) -> tuple[dict, dict[str, str]]:
    manifest_path = ledger_dir / "output_manifest.json"
    complete_path = ledger_dir / "RUN_COMPLETE.json"
    summary_path = ledger_dir / "summary.json"
    if not all(path.is_file() for path in (manifest_path, complete_path, summary_path)):
        raise AssemblyError("Source-resolved ledger is incomplete")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for record in manifest.get("files", []):
        path = ledger_dir / record["path"]
        if not path.is_file() or path.stat().st_size != int(record["bytes"]):
            raise AssemblyError(f"Ledger payload missing or resized: {record['path']}")
        if sha256_file(path) != record["sha256"]:
            raise AssemblyError(f"Ledger payload hash mismatch: {record['path']}")
    complete = json.loads(complete_path.read_text(encoding="utf-8"))
    if complete.get("complete") is not True:
        raise AssemblyError("Ledger completion marker is false")
    if complete.get("output_manifest_sha256") != sha256_file(manifest_path):
        raise AssemblyError("Ledger completion marker has the wrong manifest hash")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    source_hashes = {
        str(Path(source["archive"]).resolve()): source["archive_sha256"]
        for source in summary["sources"]
    }
    return summary, source_hashes


def load_pairs(ledger_dir: Path) -> dict[PAIR_KEY, dict[str, str]]:
    path = ledger_dir / "resolved_embedding_pairs.tsv.gz"
    pairs: dict[PAIR_KEY, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "protein_id",
            "modality",
            "action",
            "reason",
            "selected_archive",
            "selected_member",
            "array_sha256",
            "file_sha256",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise AssemblyError("Resolved pair ledger has an incomplete schema")
        for row in reader:
            key = (row["protein_id"], row["modality"])
            if key in pairs:
                raise AssemblyError(f"Pair occurs twice in ledger: {key}")
            if row["modality"] not in MODALITIES or row["action"] not in {
                "reuse",
                "regenerate",
            }:
                raise AssemblyError(f"Invalid resolved pair: {key}")
            if row["action"] == "reuse" and not all(
                row[field]
                for field in (
                    "selected_archive",
                    "selected_member",
                    "array_sha256",
                    "file_sha256",
                )
            ):
                raise AssemblyError(f"Reusable pair lacks source evidence: {key}")
            pairs[key] = row
    proteins = {protein_id for protein_id, _ in pairs}
    expected = {(protein_id, modality) for protein_id in proteins for modality in MODALITIES}
    if set(pairs) != expected:
        raise AssemblyError("Resolved ledger is not a complete protein/modality matrix")
    return pairs


def load_policy(path: Path) -> dict:
    policy = json.loads(path.read_text(encoding="utf-8"))
    if policy.get("schema_version") != 1 or set(policy.get("modalities", {})) != set(
        MODALITIES
    ):
        raise AssemblyError("Embedding policy must define exactly four modalities")
    for modality, specification in policy["modalities"].items():
        if specification.get("cache_directory") != MODALITY_DIRECTORIES[modality]:
            raise AssemblyError(f"Wrong cache directory in policy for {modality}")
        if int(specification.get("dimension", 0)) != EXPECTED_DIMENSIONS[modality]:
            raise AssemblyError(f"Wrong dimension in policy for {modality}")
    return policy


def parse_generated(values: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        modality, separator, raw_path = value.partition("=")
        if not separator or modality not in MODALITIES or modality in result:
            raise AssemblyError(f"Invalid --generated-archive value: {value}")
        path = Path(raw_path).expanduser().resolve()
        if not path.is_file() or path.is_symlink():
            raise AssemblyError(f"Generated archive is missing or unsafe: {path}")
        result[modality] = path
    if set(result) != set(MODALITIES):
        raise AssemblyError("One generated archive is required for every modality")
    return result


def deterministic_tar_info(name: str, size: int) -> tarfile.TarInfo:
    info = tarfile.TarInfo(name)
    info.size = size
    info.mode = 0o644
    info.mtime = 0
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    return info


def open_output_archive(path: Path) -> tuple[BinaryIO, gzip.GzipFile, tarfile.TarFile]:
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    archive = tarfile.open(fileobj=compressed, mode="w")
    return raw, compressed, archive


def close_output_archive(handles: tuple[BinaryIO, gzip.GzipFile, tarfile.TarFile]) -> None:
    raw, compressed, archive = handles
    archive.close()
    compressed.close()
    raw.flush()
    os.fsync(raw.fileno())
    raw.close()


def add_array(
    output: tarfile.TarFile, protein_id: str, modality: str, data: bytes
) -> str:
    name = f"data/embedding_cache/{MODALITY_DIRECTORIES[modality]}/{protein_id}.npy"
    output.addfile(deterministic_tar_info(name, len(data)), io.BytesIO(data))
    return name


def assembly_writer(path: Path):
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    fields = (
        "protein_id",
        "modality",
        "ledger_action",
        "ledger_reason",
        "status",
        "source_archive",
        "source_member",
        "output_member",
        "file_sha256",
        "array_sha256",
        "dtype",
        "dimension",
    )
    writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    return raw, compressed, text, writer


def close_writer(handles) -> None:
    raw, compressed, text, _writer = handles
    text.flush()
    text.detach()
    compressed.close()
    raw.flush()
    os.fsync(raw.fileno())
    raw.close()


def publish_cache(
    ledger_dir: Path,
    generated_archives: dict[str, Path],
    policy_path: Path,
    output_archive: Path,
    report_dir: Path,
) -> dict:
    if output_archive.exists() or report_dir.exists():
        raise AssemblyError("Output archive and report directory must not already exist")
    summary, source_hashes = verify_ledger(ledger_dir)
    pairs = load_pairs(ledger_dir)
    policy = load_policy(policy_path)
    proteins = {protein_id for protein_id, _ in pairs}
    output_archive.parent.mkdir(parents=True, exist_ok=True)
    report_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{report_dir.name}.staging-", dir=report_dir.parent))
    archive_tmp = output_archive.parent / f".{output_archive.name}.building"
    archive_tmp.unlink(missing_ok=True)
    details_path = stage / "embedding_assembly.tsv.gz"
    output_handles = open_output_archive(archive_tmp)
    report_handles = assembly_writer(details_path)
    counts: Counter[tuple[str, str, str]] = Counter()
    accepted: set[PAIR_KEY] = set()
    inputs = []

    try:
        reuse_by_archive: dict[str, dict[str, tuple[PAIR_KEY, dict[str, str]]]] = defaultdict(dict)
        for key, row in pairs.items():
            if row["action"] != "reuse":
                continue
            archive_path = str(Path(row["selected_archive"]).resolve())
            member = row["selected_member"].removeprefix("./")
            if member in reuse_by_archive[archive_path]:
                raise AssemblyError(f"Source member selected twice: {archive_path}:{member}")
            reuse_by_archive[archive_path][member] = (key, row)

        for archive_name in sorted(reuse_by_archive):
            source = Path(archive_name)
            expected_sha = source_hashes.get(archive_name)
            if expected_sha is None:
                raise AssemblyError(f"Ledger has no bound hash for source archive: {source}")
            observed_sha = sha256_file(source)
            if observed_sha != expected_sha:
                raise AssemblyError(f"Source archive hash mismatch: {source}")
            inputs.append({"role": "reuse", "path": archive_name, "sha256": observed_sha})
            wanted = reuse_by_archive[archive_name]
            seen_members: set[str] = set()
            with tarfile.open(source, mode="r:gz") as archive:
                for member in archive:
                    normalized_name = member.name.removeprefix("./")
                    selection = wanted.get(normalized_name)
                    if selection is None:
                        continue
                    normalized = normalized_member(member)
                    if normalized is None:
                        raise AssemblyError(f"Selected source member is not a file: {member.name}")
                    key, row = selection
                    modality, protein_id = normalized
                    if key != (protein_id, modality):
                        raise AssemblyError(f"Selected source member resolves to the wrong pair: {key}")
                    extracted = archive.extractfile(member)
                    if extracted is None:
                        raise AssemblyError(f"Cannot read source member: {member.name}")
                    with extracted:
                        data = extracted.read()
                    file_sha = hashlib.sha256(data).hexdigest()
                    array_sha, dtype, shape = canonical_array(data, modality)
                    if file_sha != row["file_sha256"] or array_sha != row["array_sha256"]:
                        raise AssemblyError(f"Reusable array changed since resolution: {key}")
                    output_member = add_array(output_handles[2], protein_id, modality, data)
                    report_handles[3].writerow(
                        {
                            "protein_id": protein_id,
                            "modality": modality,
                            "ledger_action": "reuse",
                            "ledger_reason": row["reason"],
                            "status": "available",
                            "source_archive": archive_name,
                            "source_member": member.name,
                            "output_member": output_member,
                            "file_sha256": file_sha,
                            "array_sha256": array_sha,
                            "dtype": dtype,
                            "dimension": shape[0],
                        }
                    )
                    accepted.add(key)
                    seen_members.add(normalized_name)
                    counts[(modality, "reuse", "available")] += 1
            missing = set(wanted) - seen_members
            if missing:
                raise AssemblyError(f"Source archive lacks {len(missing)} selected members: {source}")

        for modality in MODALITIES:
            source = generated_archives[modality]
            observed_sha = sha256_file(source)
            inputs.append(
                {"role": f"generated-{modality}", "path": str(source), "sha256": observed_sha}
            )
            seen_pairs: set[PAIR_KEY] = set()
            with tarfile.open(source, mode="r:gz") as archive:
                for member in archive:
                    normalized = normalized_member(member)
                    if normalized is None:
                        continue
                    observed_modality, protein_id = normalized
                    key = (protein_id, observed_modality)
                    if observed_modality != modality:
                        raise AssemblyError(
                            f"Generated {modality} archive contains {observed_modality}: {member.name}"
                        )
                    row = pairs.get(key)
                    if row is None or row["action"] != "regenerate":
                        raise AssemblyError(f"Generated archive contains an unrequested pair: {key}")
                    if key in seen_pairs or key in accepted:
                        raise AssemblyError(f"Generated pair occurs more than once: {key}")
                    extracted = archive.extractfile(member)
                    if extracted is None:
                        raise AssemblyError(f"Cannot read generated member: {member.name}")
                    with extracted:
                        data = extracted.read()
                    file_sha = hashlib.sha256(data).hexdigest()
                    array_sha, dtype, shape = canonical_array(data, modality)
                    output_member = add_array(output_handles[2], protein_id, modality, data)
                    report_handles[3].writerow(
                        {
                            "protein_id": protein_id,
                            "modality": modality,
                            "ledger_action": "regenerate",
                            "ledger_reason": row["reason"],
                            "status": "available",
                            "source_archive": str(source),
                            "source_member": member.name,
                            "output_member": output_member,
                            "file_sha256": file_sha,
                            "array_sha256": array_sha,
                            "dtype": dtype,
                            "dimension": shape[0],
                        }
                    )
                    accepted.add(key)
                    seen_pairs.add(key)
                    counts[(modality, "regenerate", "available")] += 1

        for key in sorted(set(pairs) - accepted):
            protein_id, modality = key
            row = pairs[key]
            if row["action"] != "regenerate":
                raise AssemblyError(f"Reusable pair was not materialized: {key}")
            report_handles[3].writerow(
                {
                    "protein_id": protein_id,
                    "modality": modality,
                    "ledger_action": "regenerate",
                    "ledger_reason": row["reason"],
                    "status": "missing-after-regeneration",
                    "source_archive": str(generated_archives[modality]),
                    "source_member": "",
                    "output_member": "",
                    "file_sha256": "",
                    "array_sha256": "",
                    "dtype": "",
                    "dimension": EXPECTED_DIMENSIONS[modality],
                }
            )
            counts[(modality, "regenerate", "missing")] += 1
    finally:
        close_writer(report_handles)
        close_output_archive(output_handles)

    modality_summary = {}
    for modality in MODALITIES:
        available = sum(
            counts[(modality, action, "available")] for action in ("reuse", "regenerate")
        )
        missing = len(proteins) - available
        specification = policy["modalities"][modality]
        minimum = specification.get("min_accepted_count")
        if minimum is None:
            minimum = int(np.ceil(float(specification["min_accepted_fraction"]) * len(proteins)))
        if available < int(minimum):
            raise AssemblyError(
                f"{modality} availability {available} is below required minimum {minimum}"
            )
        modality_summary[modality] = {
            "target_proteins": len(proteins),
            "available": available,
            "missing": missing,
            "coverage_fraction": available / len(proteins),
            "reuse_available": counts[(modality, "reuse", "available")],
            "regenerated_available": counts[(modality, "regenerate", "available")],
            "regeneration_missing": counts[(modality, "regenerate", "missing")],
            "minimum_required": int(minimum),
        }

    archive_sha = sha256_file(archive_tmp)
    result = {
        "schema_version": 1,
        "analysis_kind": "pair_resolved_homology_embedding_cache",
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "target_proteins": len(proteins),
        "target_pairs": len(pairs),
        "available_pairs": len(accepted),
        "missing_pairs": len(pairs) - len(accepted),
        "modalities": modality_summary,
        "ledger": {
            "path": str(ledger_dir.resolve()),
            "output_manifest_sha256": sha256_file(ledger_dir / "output_manifest.json"),
            "policy": summary["comparison_policy"],
        },
        "policy": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "inputs": inputs,
        "cache_archive": {
            "path": str(output_archive.resolve()),
            "sha256": archive_sha,
            "size_bytes": archive_tmp.stat().st_size,
        },
    }
    (stage / "assembly_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (stage / "RUN_COMPLETE.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "complete": True,
                "analysis_kind": result["analysis_kind"],
                "cache_archive_sha256": archive_sha,
                "assembly_summary_sha256": sha256_file(stage / "assembly_summary.json"),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    os.replace(archive_tmp, output_archive)
    os.replace(stage, report_dir)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger-dir", type=Path, required=True)
    parser.add_argument(
        "--generated-archive",
        action="append",
        required=True,
        metavar="MODALITY=PATH",
    )
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--output-archive", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    args = parser.parse_args()
    result = publish_cache(
        args.ledger_dir.resolve(),
        parse_generated(args.generated_archive),
        args.policy.resolve(),
        args.output_archive.resolve(),
        args.report_dir.resolve(),
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, json.JSONDecodeError, AssemblyError) as error:
        raise SystemExit(f"ERROR: {error}") from error
