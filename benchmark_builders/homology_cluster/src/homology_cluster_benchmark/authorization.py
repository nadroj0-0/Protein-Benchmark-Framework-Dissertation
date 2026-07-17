from __future__ import annotations

from datetime import date
import json
import math
from pathlib import Path
import re
from typing import Any

from .attrition import (
    METRIC_DEFINITIONS,
    evaluate_attrition,
    observation,
    require_reviewed_text,
)
from .inputs import sha256_file


APPROVAL_SCHEMA = "homology-cluster-pilot-approval"
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"{label} is not readable valid JSON: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{label} root must be an object")
    return payload


def _pilot_observations(attrition: dict[str, Any]) -> dict[str, dict[str, Any]]:
    metrics = attrition.get("metrics")
    if not isinstance(metrics, list):
        raise ValueError("Pilot attrition report metrics must be a list")
    by_name: dict[str, dict[str, Any]] = {}
    for item in metrics:
        if not isinstance(item, dict) or not isinstance(item.get("name"), str):
            raise ValueError("Every pilot attrition metric must be a named object")
        name = item["name"]
        if name not in METRIC_DEFINITIONS or name in by_name:
            raise ValueError(f"Pilot attrition report has an unknown or duplicate metric: {name}")
        numerator = item.get("numerator")
        denominator = item.get("denominator")
        if (
            isinstance(numerator, bool)
            or isinstance(denominator, bool)
            or not isinstance(numerator, (int, float))
            or not isinstance(denominator, (int, float))
            or not math.isfinite(numerator)
            or not math.isfinite(denominator)
            or numerator < 0
            or denominator <= 0
            or numerator > denominator
        ):
            raise ValueError(f"Pilot attrition metric {name} has invalid counts")
        rebuilt = observation(name, numerator, denominator)
        for key in (
            "ratio", "numerator_definition", "denominator_definition", "bound_type",
        ):
            observed = item.get(key)
            expected = rebuilt[key]
            if key == "ratio":
                if (
                    isinstance(observed, bool)
                    or not isinstance(observed, (int, float))
                    or not math.isfinite(observed)
                    or not math.isclose(float(observed), float(expected), rel_tol=0, abs_tol=1e-12)
                ):
                    raise ValueError(f"Pilot attrition metric {name} ratio is not reproducible")
            elif observed != expected:
                raise ValueError(f"Pilot attrition metric {name} {key} is inconsistent")
        by_name[name] = rebuilt
    if set(by_name) != set(METRIC_DEFINITIONS):
        missing = sorted(set(METRIC_DEFINITIONS) - set(by_name))
        raise ValueError("Pilot attrition report is missing registered metrics: " + ", ".join(missing))
    return by_name


def validate_pilot_attrition_against_reviewed_policy(
    attrition: dict[str, Any],
    *,
    reviewed_policy: dict[str, Any],
    reviewed_policy_sha256: str,
    source_scope: str,
    framework_commit: str,
) -> dict[str, Any]:
    input_manifest_sha256 = str(attrition.get("input_manifest_sha256", ""))
    if SHA256_RE.fullmatch(input_manifest_sha256) is None:
        raise ValueError("Pilot attrition report input_manifest_sha256 must be one SHA-256")
    if SHA256_RE.fullmatch(reviewed_policy_sha256) is None:
        raise ValueError("Reviewed attrition policy hash must be one SHA-256")
    reviewed_result = evaluate_attrition(
        reviewed_policy,
        reviewed_policy_sha256,
        _pilot_observations(attrition),
        source_scope=source_scope,
        framework_commit=framework_commit,
        input_manifest_sha256=input_manifest_sha256,
        diagnostic=True,
    )
    if reviewed_result["policy_passed"] is not True:
        failed = ", ".join(item["name"] for item in reviewed_result["failed_metrics"])
        raise ValueError(
            "Reviewed attrition policy does not accept the reviewed 30% pilot observations: "
            + failed
        )
    if reviewed_result["production_authorized"] is not False:
        raise ValueError("Diagnostic pilot policy review unexpectedly authorized production")
    return reviewed_result


def validate_pilot_approval(
    approval_path: Path,
    *,
    completion_marker_path: Path,
    attrition_report_path: Path,
    task_context_path: Path,
    measurement_evidence_path: Path,
    framework_commit: str,
    frozen_input_manifest_sha256: str,
    source_scope: str,
    split_policy: str,
    training_population: str,
    mmseqs_version: str,
    reviewed_attrition_policy: dict[str, Any],
    reviewed_attrition_policy_sha256: str,
) -> dict[str, Any]:
    approval = _json_object(approval_path, "Pilot approval")
    if approval.get("schema_name") != APPROVAL_SCHEMA or approval.get("schema_version") != 1:
        raise ValueError("Pilot approval has an unsupported schema_name/schema_version")
    if approval.get("approved") is not True:
        raise ValueError("Pilot approval must be manually changed to approved=true after review")
    marker = _json_object(completion_marker_path, "Pilot completion marker")
    attrition = _json_object(attrition_report_path, "Pilot attrition report")
    task_context = _json_object(task_context_path, "Pilot task context")
    measurement = _json_object(measurement_evidence_path, "Pilot measurement evidence")
    expected = {
        "pilot_task_id": 1,
        "pilot_identity_percent": 30,
        "successful_completion_marker_sha256": sha256_file(completion_marker_path),
        "framework_commit": framework_commit,
        "frozen_input_manifest_sha256": frozen_input_manifest_sha256,
        "uniprot_source_scope": source_scope,
        "split_policy": split_policy,
        "training_population": training_population,
        "mmseqs_version": mmseqs_version,
        "attrition_report_sha256": sha256_file(attrition_report_path),
        "reviewed_attrition_policy_sha256": reviewed_attrition_policy_sha256,
        "pilot_task_context_sha256": sha256_file(task_context_path),
        "pilot_measurement_evidence_sha256": sha256_file(measurement_evidence_path),
        "validation_outcome": "pass",
    }
    for key, value in expected.items():
        if approval.get(key) != value:
            raise ValueError(
                f"Pilot approval {key} mismatch: expected {value!r}, "
                f"observed {approval.get(key)!r}"
            )
    require_reviewed_text(approval.get("pilot_job_id"), "Pilot approval pilot_job_id")
    if (
        not isinstance(approval.get("pilot_run_id"), str)
        or not approval["pilot_run_id"].strip()
        or re.fullmatch(r"[A-Za-z0-9._-]+", approval["pilot_run_id"]) is None
        or re.search(r"[A-Za-z0-9]", approval["pilot_run_id"]) is None
    ):
        raise ValueError("Pilot approval pilot_run_id must be a safe nonempty identifier")
    for key in ("runtime_seconds", "peak_memory_bytes", "scratch_peak_bytes", "output_size_bytes"):
        value = approval.get(key)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError(f"Pilot approval {key} must be a finite positive measured value")
    require_reviewed_text(approval.get("pilot_run_id"), "Pilot approval pilot_run_id")
    for key in ("reviewer", "evidence_notes"):
        require_reviewed_text(approval.get(key), f"Pilot approval {key}")
    try:
        date.fromisoformat(str(approval.get("review_date", "")))
    except ValueError as exc:
        raise ValueError("Pilot approval review_date must be an ISO date") from exc
    if COMMIT_RE.fullmatch(framework_commit) is None:
        raise ValueError("Expected framework commit must be a full lowercase SHA")
    if SHA256_RE.fullmatch(frozen_input_manifest_sha256) is None:
        raise ValueError("Expected frozen-input manifest hash must be one SHA-256")
    marker_expectations = {
        "complete": True,
        "benchmark_scope": "diagnostic-pilot",
        "production_eligible": False,
        "identity_percent": 30,
        "framework_revision": framework_commit,
        "repository_commit": framework_commit,
        "frozen_input_manifest_sha256": frozen_input_manifest_sha256,
        "uniprot_source_scope": source_scope,
        "split_policy": split_policy,
        "training_population": training_population,
        "observed_mmseqs_version": mmseqs_version,
        "run_id": approval.get("pilot_run_id"),
    }
    for key, value in marker_expectations.items():
        if marker.get(key) != value:
            raise ValueError(f"Pilot completion marker failed required binding {key}")
    if (
        attrition.get("schema_name") != "homology-cluster-attrition-report"
        or attrition.get("schema_version") != 1
    ):
        raise ValueError("Pilot attrition report has the wrong schema")
    if attrition.get("diagnostic") is not True:
        raise ValueError("Pilot attrition report must declare diagnostic=true")
    if attrition.get("production_authorized") is not False:
        raise ValueError("Pilot attrition report must declare production_authorized=false")
    attrition_expectations = {
        "uniprot_source_scope": source_scope,
        "framework_commit": framework_commit,
        "input_manifest_sha256": marker.get("run_input_manifest_sha256"),
    }
    for key, value in attrition_expectations.items():
        if attrition.get(key) != value:
            raise ValueError(f"Pilot attrition report failed required binding {key}")
    for key in ("input_manifest_sha256", "policy_sha256"):
        if not SHA256_RE.fullmatch(str(attrition.get(key, ""))):
            raise ValueError(f"Pilot attrition report {key} must be one SHA-256")
    validate_pilot_attrition_against_reviewed_policy(
        attrition,
        reviewed_policy=reviewed_attrition_policy,
        reviewed_policy_sha256=reviewed_attrition_policy_sha256,
        source_scope=source_scope,
        framework_commit=framework_commit,
    )

    task_expectations = {
        "job_id": approval.get("pilot_job_id"),
        "sge_task_id": 1,
        "identity_percent": 30,
        "uniprot_source_scope": source_scope,
        "run_id": approval.get("pilot_run_id"),
        "framework_revision": framework_commit,
        "requested_smp_slots": 2,
        "nslots": 2,
        "mmseqs_threads": 2,
    }
    for key, value in task_expectations.items():
        if task_context.get(key) != value:
            raise ValueError(f"Pilot task context failed required binding {key}")

    if (
        measurement.get("schema_name") != "homology-cluster-pilot-measurement-evidence"
        or measurement.get("schema_version") != 1
    ):
        raise ValueError("Pilot measurement evidence has the wrong schema")
    measurement_expectations = {
        "pilot_job_id": approval.get("pilot_job_id"),
        "pilot_task_id": 1,
        "pilot_identity_percent": 30,
        "run_id": approval.get("pilot_run_id"),
        "framework_commit": framework_commit,
        "uniprot_source_scope": source_scope,
        "successful_completion_marker_sha256": sha256_file(completion_marker_path),
        "runtime_seconds": approval.get("runtime_seconds"),
        "peak_memory_bytes": approval.get("peak_memory_bytes"),
        "scratch_peak_bytes": approval.get("scratch_peak_bytes"),
        "output_size_bytes": approval.get("output_size_bytes"),
    }
    for key, value in measurement_expectations.items():
        if measurement.get(key) != value:
            raise ValueError(f"Pilot measurement evidence failed required binding {key}")
    sources = measurement.get("measurement_sources")
    required_sources = {"runtime", "peak_memory", "scratch_peak", "output_size"}
    if not isinstance(sources, dict) or set(sources) != required_sources:
        raise ValueError("Pilot measurement evidence must name every measurement source")
    for key, value in sources.items():
        require_reviewed_text(value, f"Pilot measurement source {key}")
    for key in ("reviewer", "evidence_notes"):
        require_reviewed_text(
            measurement.get(key), f"Pilot measurement evidence {key}"
        )
    _review_date = str(measurement.get("review_date", ""))
    try:
        date.fromisoformat(_review_date)
    except ValueError as exc:
        raise ValueError("Pilot measurement evidence review_date must be an ISO date") from exc
    return approval
