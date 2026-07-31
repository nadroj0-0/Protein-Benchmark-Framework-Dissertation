from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .inputs import sha256_file


SCHEMA_NAME = "homology-external-cluster-assignments"
SCHEMA_VERSION = 1


def _mapping(payload: dict[str, Any], key: str) -> dict[str, Any]:
    value = payload.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"External cluster provenance lacks a {key!r} object")
    return value


def _float_equal(observed: object, expected: float) -> bool:
    try:
        return abs(float(observed) - expected) < 1e-12
    except (TypeError, ValueError):
        return False


def load_external_cluster_provenance(
    provenance_path: Path,
    assignments_path: Path,
    *,
    identity: float,
    coverage: float,
    cov_mode: int,
    cluster_mode: int,
    sensitivity: float,
    uniref_level: int,
    uniref_release: str,
    uniref_sha256: str,
    uniref_records: int,
) -> dict[str, Any]:
    provenance = provenance_path.expanduser().resolve()
    assignments = assignments_path.expanduser().resolve()
    if not provenance.is_file() or not assignments.is_file():
        raise FileNotFoundError(provenance if not provenance.is_file() else assignments)
    try:
        payload = json.loads(provenance.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ValueError(f"External cluster provenance is not valid UTF-8 JSON: {provenance}") from exc
    if not isinstance(payload, dict):
        raise ValueError("External cluster provenance root must be a JSON object")
    if payload.get("schema_name") != SCHEMA_NAME or payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("External cluster provenance has an unsupported schema")

    artifact = _mapping(payload, "artifact")
    expected_assignment_sha = str(artifact.get("sha256", ""))
    observed_assignment_sha = sha256_file(assignments)
    if observed_assignment_sha != expected_assignment_sha:
        raise ValueError(
            "External cluster assignment SHA-256 mismatch: "
            f"expected={expected_assignment_sha} observed={observed_assignment_sha}"
        )
    if artifact.get("size_bytes") != assignments.stat().st_size:
        raise ValueError("External cluster assignment byte size does not match provenance")

    input_record = _mapping(payload, "input")
    expected_input = {
        "uniref_level": uniref_level,
        "release": uniref_release,
        "expected_fasta_sha256": uniref_sha256,
        "expected_records": uniref_records,
    }
    for key, expected in expected_input.items():
        if input_record.get(key) != expected:
            raise ValueError(
                f"External cluster provenance input mismatch for {key}: "
                f"expected={expected!r} observed={input_record.get(key)!r}"
            )

    method = _mapping(payload, "method")
    numeric_method = {
        "identity_percent": identity * 100,
        "coverage": coverage,
        "coverage_mode": float(cov_mode),
        "cluster_mode": float(cluster_mode),
        "sensitivity": sensitivity,
    }
    for key, expected in numeric_method.items():
        if not _float_equal(method.get(key), expected):
            raise ValueError(
                f"External cluster provenance method mismatch for {key}: "
                f"expected={expected!r} observed={method.get(key)!r}"
            )
    for key in ("evalue", "createdb_shuffle", "cluster_reassignment"):
        if method.get(key) != "MMseqs2 default":
            raise ValueError(
                f"External Daniel-aligned artifact must declare {key}=MMseqs2 default"
            )

    usage = _mapping(payload, "usage_policy")
    if usage.get("lineage") != "supervisor-generated":
        raise ValueError("External cluster provenance must identify supervisor-generated lineage")
    if usage.get("do_not_merge_with_framework_generated_cluster_cache") is not True:
        raise ValueError("External cluster provenance must prohibit framework-cache merging")
    return payload


def validate_external_cluster_counts(
    payload: dict[str, Any], *, members: int, clusters: int
) -> None:
    artifact = _mapping(payload, "artifact")
    if artifact.get("members") != members or artifact.get("clusters") != clusters:
        raise ValueError(
            "External cluster assignment counts do not match provenance: "
            f"expected_members={artifact.get('members')} observed_members={members}, "
            f"expected_clusters={artifact.get('clusters')} observed_clusters={clusters}"
        )
