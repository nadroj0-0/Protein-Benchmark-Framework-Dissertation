#!/usr/bin/env python3
"""Resolve a benchmark-level reuse plan to exact PFP cache members.

The benchmark reuse planner proves that a target protein has the same ID and
sequence in one or more earlier benchmark populations.  This tool performs the
next, deliberately separate step: it inventories authenticated cache archives,
compares duplicate arrays exactly, and emits one reuse/regenerate decision for
every target protein/modality pair.
"""

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
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

import numpy as np


MODALITIES: Tuple[str, ...] = ("sequence", "text", "structure", "ppi")
DIRECTORY_TO_MODALITY = {
    "prott5": "sequence",
    "exp_text_embeddings_temporal": "text",
    "IF1": "structure",
    "ppi": "ppi",
}
EXPECTED_DIMENSIONS = {"sequence": 1024, "text": 768, "structure": 512, "ppi": 512}
ARCHIVE_PREFIX = PurePosixPath("data/embedding_cache")
SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SAFE_ID = re.compile(r"^[^\s/\\]+$")
PLAN_COLUMNS = {
    "protein_id",
    "sequence",
    "sequence_sha256",
    "action",
    "reason",
    "matching_embedded_benchmarks",
    "embedded_benchmark_memberships",
    "target_memberships",
}
TEXT_REUSE_POLICIES = {"never", "same-role", "source-current"}
PAIR_COLUMNS = (
    "protein_id",
    "modality",
    "action",
    "reason",
    "sequence_sha256",
    "coarse_action",
    "matching_embedded_benchmarks",
    "eligible_sources",
    "present_valid_sources",
    "selected_source",
    "selected_archive",
    "selected_member",
    "array_sha256",
    "file_sha256",
    "agreeing_sources",
    "conflicting_array_sha256s",
    "source_issues",
)
CANDIDATE_COLUMNS = (
    "protein_id",
    "modality",
    "source_name",
    "embedded_benchmark",
    "archive",
    "archive_sha256",
    "member",
    "status",
    "array_sha256",
    "file_sha256",
    "shape",
    "dtype",
    "detail",
)


class ResolutionError(RuntimeError):
    pass


@dataclass(frozen=True)
class SourceSpec:
    name: str
    embedded_benchmark: str
    archive: Path
    archive_sha256: str
    priority: int
    text_reuse_policy: str = "never"


@dataclass(frozen=True)
class PlanProtein:
    protein_id: str
    sequence: str
    sequence_sha256: str
    coarse_action: str
    coarse_reason: str
    matching_benchmarks: Tuple[str, ...]
    embedded_memberships: Tuple[str, ...]
    target_memberships: Tuple[str, ...]


@dataclass(frozen=True)
class Candidate:
    protein_id: str
    modality: str
    source: SourceSpec
    member: str
    status: str
    array_sha256: str = ""
    file_sha256: str = ""
    shape: str = ""
    dtype: str = ""
    detail: str = ""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def json_list(value: str, *, field: str) -> Tuple[str, ...]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise ResolutionError(f"Invalid JSON in {field}: {error}") from error
    if not isinstance(parsed, list) or not all(isinstance(item, str) for item in parsed):
        raise ResolutionError(f"{field} must be a JSON list of strings")
    if parsed != sorted(set(parsed)):
        raise ResolutionError(f"{field} must be sorted and unique")
    return tuple(parsed)


def parse_source(
    value: str,
    priority: int,
    text_reuse_policy: str = "never",
) -> SourceSpec:
    fields = value.split("=", 3)
    if len(fields) != 4 or any(not field for field in fields):
        raise ResolutionError(
            "--cache-source must use NAME=EMBEDDED_BENCHMARK=ARCHIVE=SHA256"
        )
    name, benchmark, archive_value, digest = fields
    if not SAFE_NAME.fullmatch(name) or not SAFE_NAME.fullmatch(benchmark):
        raise ResolutionError(f"Unsafe source or embedded benchmark name: {value}")
    if not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
        raise ResolutionError(f"Invalid source archive SHA-256: {digest}")
    if text_reuse_policy not in TEXT_REUSE_POLICIES:
        raise ResolutionError(
            f"Unsupported text reuse policy for {name}: {text_reuse_policy}"
        )
    archive = Path(archive_value).expanduser().resolve()
    if not archive.is_file() or archive.is_symlink():
        raise ResolutionError(f"Source archive is missing or unsafe: {archive}")
    return SourceSpec(
        name,
        benchmark,
        archive,
        digest.lower(),
        priority,
        text_reuse_policy,
    )


def parse_text_policies(values: Sequence[str]) -> Dict[str, str]:
    policies: Dict[str, str] = {}
    for value in values:
        fields = value.split("=", 1)
        if len(fields) != 2 or any(not field for field in fields):
            raise ResolutionError("--source-text-policy must use SOURCE=POLICY")
        name, policy = fields
        if not SAFE_NAME.fullmatch(name):
            raise ResolutionError(f"Unsafe source name in text policy: {name}")
        if policy not in TEXT_REUSE_POLICIES:
            raise ResolutionError(
                f"Unsupported text reuse policy for {name}: {policy}; expected one of "
                + ", ".join(sorted(TEXT_REUSE_POLICIES))
            )
        if name in policies:
            raise ResolutionError(f"Repeated text reuse policy for source: {name}")
        policies[name] = policy
    return policies


def verify_plan_manifest(plan_dir: Path) -> str:
    manifest_path = plan_dir / "output_manifest.json"
    complete_path = plan_dir / "RUN_COMPLETE.json"
    if not manifest_path.is_file() or not complete_path.is_file():
        raise ResolutionError("Coarse plan lacks output_manifest.json or RUN_COMPLETE.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    files = manifest.get("files")
    if not isinstance(files, list) or not files:
        raise ResolutionError("Coarse-plan output manifest has no files")
    for item in files:
        if not isinstance(item, dict):
            raise ResolutionError("Coarse-plan manifest contains a non-object file entry")
        relative = item.get("relative_path", item.get("path"))
        if not isinstance(relative, str):
            raise ResolutionError("Coarse-plan manifest contains an invalid path")
        path = plan_dir / relative
        if not path.is_file():
            raise ResolutionError(f"Coarse-plan payload is missing: {relative}")
        expected_size = item.get("size_bytes", item.get("bytes"))
        if path.stat().st_size != expected_size:
            raise ResolutionError(f"Coarse-plan payload size mismatch: {relative}")
        if sha256_file(path) != item.get("sha256"):
            raise ResolutionError(f"Coarse-plan payload hash mismatch: {relative}")
    complete = json.loads(complete_path.read_text(encoding="utf-8"))
    if complete.get("complete") is not True:
        raise ResolutionError("Coarse-plan completion marker is not complete")
    manifest_identity = complete.get("output_manifest")
    if isinstance(manifest_identity, dict):
        if manifest_identity.get("sha256") != sha256_file(manifest_path):
            raise ResolutionError("Coarse-plan completion marker has the wrong manifest hash")
        if manifest_identity.get("size_bytes") != manifest_path.stat().st_size:
            raise ResolutionError("Coarse-plan completion marker has the wrong manifest size")
    return sha256_file(manifest_path)


def read_plan_file(path: Path, expected_action: str) -> Iterator[PlanProtein]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not PLAN_COLUMNS.issubset(reader.fieldnames):
            raise ResolutionError(f"Coarse-plan columns are incomplete: {path}")
        for line_number, row in enumerate(reader, start=2):
            if row["action"] != expected_action:
                raise ResolutionError(f"Unexpected action at {path}:{line_number}")
            protein_id = row["protein_id"]
            sequence = row["sequence"]
            if not SAFE_ID.fullmatch(protein_id) or not sequence:
                raise ResolutionError(f"Unsafe or empty protein at {path}:{line_number}")
            observed_sequence_hash = sequence_sha256(sequence)
            if observed_sequence_hash != row["sequence_sha256"]:
                raise ResolutionError(f"Sequence hash mismatch at {path}:{line_number}")
            yield PlanProtein(
                protein_id=protein_id,
                sequence=sequence,
                sequence_sha256=observed_sequence_hash,
                coarse_action=expected_action,
                coarse_reason=row["reason"],
                matching_benchmarks=json_list(
                    row["matching_embedded_benchmarks"],
                    field="matching_embedded_benchmarks",
                ),
                embedded_memberships=json_list(
                    row["embedded_benchmark_memberships"],
                    field="embedded_benchmark_memberships",
                ),
                target_memberships=json_list(
                    row["target_memberships"], field="target_memberships"
                ),
            )


def load_plan(plan_dir: Path) -> Dict[str, PlanProtein]:
    proteins: Dict[str, PlanProtein] = {}
    for filename, action in (
        ("reuse_proteins.tsv", "reuse"),
        ("regenerate_proteins.tsv", "regenerate"),
    ):
        path = plan_dir / filename
        if not path.is_file():
            raise ResolutionError(f"Coarse-plan payload is missing: {filename}")
        for protein in read_plan_file(path, action):
            if protein.protein_id in proteins:
                raise ResolutionError(f"Coarse plan repeats protein ID: {protein.protein_id}")
            proteins[protein.protein_id] = protein
    if not proteins:
        raise ResolutionError("Coarse plan contains no proteins")
    return proteins


def normalized_member(member: tarfile.TarInfo) -> Optional[Tuple[str, str]]:
    name = member.name.removeprefix("./")
    path = PurePosixPath(name)
    allowed_directories = {
        PurePosixPath("data"),
        ARCHIVE_PREFIX,
        *(ARCHIVE_PREFIX / directory for directory in DIRECTORY_TO_MODALITY),
    }
    if member.isdir() and path in allowed_directories:
        return None
    if not member.isfile():
        raise ResolutionError(f"Archive member is not a regular file: {member.name}")
    if len(path.parts) != 4 or PurePosixPath(*path.parts[:2]) != ARCHIVE_PREFIX:
        raise ResolutionError(f"Archive member is outside the PFP cache: {member.name}")
    directory, filename = path.parts[2], path.parts[3]
    if directory not in DIRECTORY_TO_MODALITY:
        raise ResolutionError(f"Unknown PFP cache directory in archive: {member.name}")
    if (
        PurePosixPath(filename).name != filename
        or not filename.endswith(".npy")
        or not SAFE_ID.fullmatch(filename[:-4])
    ):
        raise ResolutionError(f"Unsafe embedding archive member: {member.name}")
    return DIRECTORY_TO_MODALITY[directory], filename[:-4]


def temporal_text_role(memberships: Sequence[str]) -> str:
    return "historical" if any(item.endswith("-test.csv") for item in memberships) else "current"


def source_memberships(protein: PlanProtein, source: SourceSpec) -> Tuple[str, ...]:
    prefix = source.embedded_benchmark + ":"
    return tuple(
        item[len(prefix):]
        for item in protein.embedded_memberships
        if item.startswith(prefix)
    )


def source_is_eligible(
    protein: PlanProtein,
    modality: str,
    source: SourceSpec,
) -> bool:
    if (
        protein.coarse_action != "reuse"
        or source.embedded_benchmark not in protein.matching_benchmarks
    ):
        return False
    if modality != "text":
        return True
    if source.text_reuse_policy == "never":
        return False
    memberships = source_memberships(protein, source)
    if not memberships:
        raise ResolutionError(
            f"Reusable protein {protein.protein_id} has no recorded memberships for "
            f"source benchmark {source.embedded_benchmark}"
        )
    source_role = temporal_text_role(memberships)
    if source.text_reuse_policy == "source-current":
        return source_role == "current"
    if source.text_reuse_policy == "same-role":
        if not protein.target_memberships:
            raise ResolutionError(
                f"Reusable protein {protein.protein_id} has no target memberships"
            )
        return source_role == temporal_text_role(protein.target_memberships)
    raise ResolutionError(
        f"Unsupported text reuse policy for {source.name}: {source.text_reuse_policy}"
    )


def canonical_array_sha256(array: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(contiguous.dtype.str.encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(list(contiguous.shape), separators=(",", ":")).encode("ascii"))
    digest.update(b"\0")
    digest.update(contiguous.view(np.uint8).tobytes())
    return digest.hexdigest()


def inspect_array(data: bytes, modality: str) -> Tuple[str, str, str]:
    try:
        array = np.load(io.BytesIO(data), allow_pickle=False)
    except Exception as error:  # numpy emits several exception types for malformed NPY
        raise ResolutionError(f"unreadable-npy: {type(error).__name__}: {error}") from error
    expected_shape = (EXPECTED_DIMENSIONS[modality],)
    if array.shape != expected_shape:
        raise ResolutionError(f"wrong-shape: expected {expected_shape}, observed {array.shape}")
    if array.dtype.kind not in "fc":
        raise ResolutionError(f"wrong-dtype: expected floating point, observed {array.dtype}")
    if not np.isfinite(array).all():
        raise ResolutionError("non-finite-array")
    return canonical_array_sha256(array), json.dumps(list(array.shape)), str(array.dtype)


def scan_source(
    source: SourceSpec,
    proteins: Mapping[str, PlanProtein],
) -> List[Candidate]:
    observed_archive_sha = sha256_file(source.archive)
    if observed_archive_sha != source.archive_sha256:
        raise ResolutionError(
            f"Source archive hash mismatch for {source.name}: "
            f"{observed_archive_sha} != {source.archive_sha256}"
        )
    candidates: List[Candidate] = []
    seen: set[Tuple[str, str]] = set()
    with tarfile.open(source.archive, mode="r:gz") as archive:
        for member in archive:
            normalized = normalized_member(member)
            if normalized is None:
                continue
            modality, protein_id = normalized
            key = (protein_id, modality)
            if key in seen:
                raise ResolutionError(
                    f"Source archive repeats protein/modality member: {source.name} {key}"
                )
            seen.add(key)
            protein = proteins.get(protein_id)
            if protein is None or not source_is_eligible(protein, modality, source):
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ResolutionError(f"Cannot read archive member: {member.name}")
            with extracted:
                data = extracted.read()
            file_hash = hashlib.sha256(data).hexdigest()
            try:
                array_hash, shape, dtype = inspect_array(data, modality)
            except ResolutionError as error:
                candidates.append(
                    Candidate(
                        protein_id,
                        modality,
                        source,
                        member.name,
                        "invalid",
                        file_sha256=file_hash,
                        detail=str(error),
                    )
                )
            else:
                candidates.append(
                    Candidate(
                        protein_id,
                        modality,
                        source,
                        member.name,
                        "valid",
                        array_sha256=array_hash,
                        file_sha256=file_hash,
                        shape=shape,
                        dtype=dtype,
                    )
                )
    return candidates


def compact_json(values: Iterable[str]) -> str:
    return json.dumps(sorted(set(values)), separators=(",", ":"))


def candidate_row(candidate: Candidate) -> Dict[str, str]:
    return {
        "protein_id": candidate.protein_id,
        "modality": candidate.modality,
        "source_name": candidate.source.name,
        "embedded_benchmark": candidate.source.embedded_benchmark,
        "archive": str(candidate.source.archive),
        "archive_sha256": candidate.source.archive_sha256,
        "member": candidate.member,
        "status": candidate.status,
        "array_sha256": candidate.array_sha256,
        "file_sha256": candidate.file_sha256,
        "shape": candidate.shape,
        "dtype": candidate.dtype,
        "detail": candidate.detail,
    }


def resolve_pair(
    protein: PlanProtein,
    modality: str,
    sources: Sequence[SourceSpec],
    candidates: Sequence[Candidate],
) -> Dict[str, str]:
    eligible_sources = [
        source.name
        for source in sources
        if source_is_eligible(protein, modality, source)
    ]
    base = {
        "protein_id": protein.protein_id,
        "modality": modality,
        "sequence_sha256": protein.sequence_sha256,
        "coarse_action": protein.coarse_action,
        "matching_embedded_benchmarks": compact_json(protein.matching_benchmarks),
        "eligible_sources": compact_json(eligible_sources),
        "present_valid_sources": "[]",
        "selected_source": "",
        "selected_archive": "",
        "selected_member": "",
        "array_sha256": "",
        "file_sha256": "",
        "agreeing_sources": "[]",
        "conflicting_array_sha256s": "[]",
        "source_issues": "[]",
    }
    if protein.coarse_action == "regenerate":
        return {
            **base,
            "action": "regenerate",
            "reason": f"coarse-plan-{protein.coarse_reason}",
        }
    invalid = [candidate for candidate in candidates if candidate.status != "valid"]
    valid = [candidate for candidate in candidates if candidate.status == "valid"]
    base["present_valid_sources"] = compact_json(candidate.source.name for candidate in valid)
    base["source_issues"] = compact_json(
        f"{candidate.source.name}:{candidate.detail}" for candidate in invalid
    )
    if invalid:
        return {**base, "action": "regenerate", "reason": "invalid-source-array"}
    if not valid:
        reason = "no-valid-source-array" if eligible_sources else "no-compatible-source"
        return {**base, "action": "regenerate", "reason": reason}
    hashes = sorted({candidate.array_sha256 for candidate in valid})
    if len(hashes) != 1:
        base["conflicting_array_sha256s"] = compact_json(hashes)
        return {
            **base,
            "action": "regenerate",
            "reason": "conflicting-source-arrays",
        }
    selected = min(valid, key=lambda item: item.source.priority)
    agreeing = [candidate.source.name for candidate in valid]
    return {
        **base,
        "action": "reuse",
        "reason": "single-valid-source" if len(valid) == 1 else "identical-source-arrays",
        "selected_source": selected.source.name,
        "selected_archive": str(selected.source.archive),
        "selected_member": selected.member,
        "array_sha256": selected.array_sha256,
        "file_sha256": selected.file_sha256,
        "agreeing_sources": compact_json(agreeing),
    }


def open_gzip_tsv(path: Path, fieldnames: Sequence[str]):
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    writer = csv.DictWriter(text, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    return raw, compressed, text, writer


def close_gzip_tsv(handles: Tuple[Any, ...]) -> None:
    raw, compressed, text, _writer = handles
    text.flush()
    text.detach()
    compressed.close()
    raw.flush()
    os.fsync(raw.fileno())
    raw.close()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def output_manifest(root: Path) -> Dict[str, Any]:
    files = []
    for path in sorted(root.iterdir(), key=lambda item: item.name):
        if not path.is_file() or path.name in {"output_manifest.json", "RUN_COMPLETE.json"}:
            continue
        files.append(
            {
                "path": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return {
        "schema_version": 1,
        "payload_file_count": len(files),
        "payload_bytes": sum(item["bytes"] for item in files),
        "files": files,
    }


def publish_resolution(
    plan_dir: Path,
    output_dir: Path,
    sources: Sequence[SourceSpec],
) -> Path:
    if output_dir.exists() or output_dir.is_symlink():
        raise ResolutionError(f"Refusing to overwrite output: {output_dir}")
    plan_manifest_sha = verify_plan_manifest(plan_dir)
    proteins = load_plan(plan_dir)
    if len({source.name for source in sources}) != len(sources):
        raise ResolutionError("Cache source names must be unique")
    reusable_benchmark_names = {
        benchmark
        for protein in proteins.values()
        if protein.coarse_action == "reuse"
        for benchmark in protein.matching_benchmarks
    }
    unused_sources = [
        source
        for source in sources
        if source.embedded_benchmark not in reusable_benchmark_names
    ]
    if unused_sources:
        configured = ", ".join(
            f"{source.name}={source.embedded_benchmark}" for source in unused_sources
        )
        available = ", ".join(sorted(reusable_benchmark_names)) or "<none>"
        raise ResolutionError(
            "Cache source embedded-benchmark label is absent from every reusable "
            f"coarse-plan row: {configured}; available labels: {available}"
        )

    candidates: List[Candidate] = []
    for source in sources:
        print(f"Authenticating and scanning source: {source.name}", flush=True)
        candidates.extend(scan_source(source, proteins))
    candidate_index: Dict[Tuple[str, str], List[Candidate]] = defaultdict(list)
    for candidate in candidates:
        candidate_index[(candidate.protein_id, candidate.modality)].append(candidate)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.staging-", dir=str(output_dir.parent)))
    try:
        candidate_handles = open_gzip_tsv(
            stage / "source_candidates.tsv.gz", CANDIDATE_COLUMNS
        )
        try:
            for candidate in sorted(
                candidates,
                key=lambda item: (
                    item.protein_id,
                    MODALITIES.index(item.modality),
                    item.source.priority,
                ),
            ):
                candidate_handles[3].writerow(candidate_row(candidate))
        finally:
            close_gzip_tsv(candidate_handles)

        pair_handles = open_gzip_tsv(
            stage / "resolved_embedding_pairs.tsv.gz", PAIR_COLUMNS
        )
        conflict_handles = open_gzip_tsv(
            stage / "conflicting_embedding_pairs.tsv.gz", PAIR_COLUMNS
        )
        action_counts: Counter[str] = Counter()
        reason_counts: Counter[str] = Counter()
        modality_action_counts: Counter[Tuple[str, str]] = Counter()
        selected_source_counts: Counter[str] = Counter()
        protein_rows: List[Tuple[PlanProtein, List[Dict[str, str]]]] = []
        try:
            for protein_id in sorted(proteins):
                protein = proteins[protein_id]
                rows = []
                for modality in MODALITIES:
                    row = resolve_pair(
                        protein,
                        modality,
                        sources,
                        candidate_index.get((protein_id, modality), ()),
                    )
                    pair_handles[3].writerow(row)
                    if row["reason"] == "conflicting-source-arrays":
                        conflict_handles[3].writerow(row)
                    rows.append(row)
                    action_counts[row["action"]] += 1
                    reason_counts[row["reason"]] += 1
                    modality_action_counts[(modality, row["action"])] += 1
                    if row["selected_source"]:
                        selected_source_counts[row["selected_source"]] += 1
                protein_rows.append((protein, rows))
        finally:
            close_gzip_tsv(pair_handles)
            close_gzip_tsv(conflict_handles)

        protein_columns = (
            "protein_id",
            "sequence",
            "sequence_sha256",
            "action",
            "reason",
            "reuse_modalities",
            "regenerate_modalities",
            "coarse_action",
            "target_memberships",
        )
        all_reuse = 0
        any_regenerate = 0
        all_regenerate = 0
        with (stage / "reuse_proteins.tsv").open("w", encoding="utf-8", newline="") as reuse_file, (
            stage / "regenerate_proteins.tsv"
        ).open("w", encoding="utf-8", newline="") as regenerate_file, (
            stage / "regenerate_proteins.txt"
        ).open("w", encoding="utf-8") as regenerate_ids, (
            stage / "regenerate_proteins.fasta"
        ).open("w", encoding="utf-8") as regenerate_fasta:
            reuse_writer = csv.DictWriter(
                reuse_file,
                fieldnames=protein_columns,
                delimiter="\t",
                lineterminator="\n",
            )
            regenerate_writer = csv.DictWriter(
                regenerate_file,
                fieldnames=protein_columns,
                delimiter="\t",
                lineterminator="\n",
            )
            reuse_writer.writeheader()
            regenerate_writer.writeheader()
            for protein, rows in protein_rows:
                reuse_modalities = [
                    row["modality"] for row in rows if row["action"] == "reuse"
                ]
                regenerate_modalities = [
                    row["modality"] for row in rows if row["action"] == "regenerate"
                ]
                row = {
                    "protein_id": protein.protein_id,
                    "sequence": protein.sequence,
                    "sequence_sha256": protein.sequence_sha256,
                    "action": "reuse" if not regenerate_modalities else "regenerate",
                    "reason": (
                        "all-modalities-reusable"
                        if not regenerate_modalities
                        else "one-or-more-modalities-require-regeneration"
                    ),
                    "reuse_modalities": compact_json(reuse_modalities),
                    "regenerate_modalities": compact_json(regenerate_modalities),
                    "coarse_action": protein.coarse_action,
                    "target_memberships": compact_json(protein.target_memberships),
                }
                if regenerate_modalities:
                    regenerate_writer.writerow(row)
                    regenerate_ids.write(protein.protein_id + "\n")
                    regenerate_fasta.write(f">{protein.protein_id}\n{protein.sequence}\n")
                    any_regenerate += 1
                    if len(regenerate_modalities) == len(MODALITIES):
                        all_regenerate += 1
                else:
                    reuse_writer.writerow(row)
                    all_reuse += 1

        summary = {
            "schema_version": 1,
            "analysis_kind": "embedding_source_resolved_reuse_ledger",
            "comparison_policy": {
                "granularity": "protein-modality",
                "array_equality": "exact dtype, shape, and numeric bytes",
                "conflicting_valid_sources": "regenerate affected modality",
                "invalid_candidate_source": "regenerate affected modality",
                "identical_duplicates": "reuse deterministic first source and record all copies",
                "source_priority": [source.name for source in sources],
                "text_reuse_policy": {
                    source.name: source.text_reuse_policy for source in sources
                },
            },
            "counts": {
                "target_proteins": len(proteins),
                "target_pairs": len(proteins) * len(MODALITIES),
                "all_modalities_reusable_proteins": all_reuse,
                "proteins_requiring_any_regeneration": any_regenerate,
                "proteins_requiring_all_modalities_regeneration": all_regenerate,
                "pair_actions": dict(sorted(action_counts.items())),
                "pair_reasons": dict(sorted(reason_counts.items())),
                "modality_actions": {
                    modality: {
                        action: modality_action_counts[(modality, action)]
                        for action in ("reuse", "regenerate")
                    }
                    for modality in MODALITIES
                },
                "selected_source_pairs": dict(sorted(selected_source_counts.items())),
                "source_candidate_records": len(candidates),
            },
            "coarse_plan": {
                "path": str(plan_dir),
                "output_manifest_sha256": plan_manifest_sha,
            },
            "sources": [
                {
                    "name": source.name,
                    "embedded_benchmark": source.embedded_benchmark,
                    "archive": str(source.archive),
                    "archive_sha256": source.archive_sha256,
                    "priority": source.priority,
                    "text_reuse_policy": source.text_reuse_policy,
                }
                for source in sources
            ],
        }
        write_json(stage / "summary.json", summary)
        lines = [
            "# Source-resolved embedding reuse ledger",
            "",
            "Every target protein/modality pair is tied to an exact authenticated archive member.",
            "Exact duplicate arrays are reusable; conflicting or invalid source arrays are regenerated.",
            "",
            "## Counts",
            "",
            f"- Target proteins: {len(proteins):,}",
            f"- Target protein/modality pairs: {len(proteins) * len(MODALITIES):,}",
            f"- Reusable pairs: {action_counts['reuse']:,}",
            f"- Regenerate pairs: {action_counts['regenerate']:,}",
            f"- Fully reusable proteins: {all_reuse:,}",
            f"- Proteins requiring at least one regenerated modality: {any_regenerate:,}",
            "",
            "## Per modality",
            "",
            "| Modality | Reuse | Regenerate |",
            "|---|---:|---:|",
        ]
        for modality in MODALITIES:
            lines.append(
                f"| {modality} | {modality_action_counts[(modality, 'reuse')]:,} | "
                f"{modality_action_counts[(modality, 'regenerate')]:,} |"
            )
        lines.extend(
            [
                "",
                "`resolved_embedding_pairs.tsv.gz` is the authoritative pair-level ledger.",
                "`regenerate_proteins.fasta` contains every protein requiring one or more fresh modalities.",
                "",
            ]
        )
        (stage / "summary.md").write_text("\n".join(lines), encoding="utf-8")
        write_json(
            stage / "run_manifest.json",
            {
                "schema_version": 1,
                "created_at": datetime.now(timezone.utc).isoformat(),
                "tool": "resolve_embedding_reuse_sources.py",
                "coarse_plan": str(plan_dir),
                "coarse_plan_output_manifest_sha256": plan_manifest_sha,
                "sources": summary["sources"],
                "output_dir": str(output_dir),
            },
        )
        manifest = output_manifest(stage)
        write_json(stage / "output_manifest.json", manifest)
        manifest_sha = sha256_file(stage / "output_manifest.json")
        write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "analysis_kind": "embedding_source_resolved_reuse_ledger",
                "output_manifest": "output_manifest.json",
                "output_manifest_sha256": manifest_sha,
                "counts": summary["counts"],
            },
        )
        if output_dir.exists():
            raise ResolutionError(f"Output appeared during publication: {output_dir}")
        os.rename(stage, output_dir)
    finally:
        if stage.exists():
            shutil.rmtree(stage, ignore_errors=True)
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse-plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--cache-source",
        action="append",
        required=True,
        metavar="NAME=EMBEDDED_BENCHMARK=ARCHIVE=SHA256",
    )
    parser.add_argument(
        "--source-text-policy",
        action="append",
        required=True,
        metavar="SOURCE=never|same-role|source-current",
        help=(
            "Declare how each source's mixed temporal text cache may be reused. "
            "One policy is required for every --cache-source."
        ),
    )
    args = parser.parse_args()
    plan_dir = args.coarse_plan_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not plan_dir.is_dir() or plan_dir.is_symlink():
        raise ResolutionError(f"Coarse plan is missing or unsafe: {plan_dir}")
    text_policies = parse_text_policies(args.source_text_policy)
    source_names = [value.split("=", 1)[0] for value in args.cache_source]
    missing_policies = sorted(set(source_names) - set(text_policies))
    unknown_policies = sorted(set(text_policies) - set(source_names))
    if missing_policies or unknown_policies:
        raise ResolutionError(
            "Text reuse policies must exactly match cache sources; "
            f"missing={missing_policies}, unknown={unknown_policies}"
        )
    sources = [
        parse_source(value, index, text_policies[source_names[index]])
        for index, value in enumerate(args.cache_source)
    ]
    result = publish_resolution(plan_dir, output_dir, sources)
    print(f"Published source-resolved ledger: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
