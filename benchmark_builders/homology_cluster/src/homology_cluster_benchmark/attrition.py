from __future__ import annotations

from dataclasses import dataclass
from datetime import date
import json
from pathlib import Path
import re
from typing import Any

from .inputs import sha256_file


POLICY_SCHEMA = "homology-cluster-attrition-policy"
OVERRIDE_SCHEMA = "homology-cluster-attrition-override"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
PLACEHOLDER_PREFIXES = ("REPLACE", "TODO", "TBD", "PLACEHOLDER")


@dataclass(frozen=True)
class MetricDefinition:
    numerator: str
    denominator: str
    bound: str


METRIC_DEFINITIONS: dict[str, MetricDefinition] = {
    "goa_to_selected_uniprot_mapping_ratio": MetricDefinition(
        "qualifying GOA accessions resolving to one sequence in the selected UniProt scope",
        "all qualifying GOA accessions after evidence, NOT, object, and ontology filtering",
        "minimum",
    ),
    "selected_uniprot_to_uniref90_mapping_ratio": MetricDefinition(
        "selected-UniProt accessions mapped uniquely to a UniRef90 FASTA member",
        "qualifying GOA accessions resolving to one selected-UniProt sequence",
        "minimum",
    ),
    "qualifying_annotation_retention_ratio": MetricDefinition(
        "qualifying annotation rows whose own accession completes the mapping and retained-cluster chain",
        "qualifying annotation rows entering mapping and label construction",
        "minimum",
    ),
    "retained_cluster_member_ratio": MetricDefinition(
        "UniRef90 members in clusters retained by a qualifying selected-UniProt annotation",
        "all frozen UniRef90 FASTA members clustered by MMseqs2",
        "minimum",
    ),
    "evaluable_protein_ratio": MetricDefinition(
        "qualifying raw GOA accessions producing at least one evaluable propagated PFP term",
        "qualifying raw GOA protein accessions entering the ontology and mapping chain",
        "minimum",
    ),
    "propagated_term_evaluable_ratio": MetricDefinition(
        "mapped labelled protein rows retaining at least one development-universe term",
        "mapped labelled protein rows with propagated GO annotations",
        "minimum",
    ),
    "bp_evaluable_ratio": MetricDefinition(
        "mapped labelled protein rows with at least one evaluable biological-process term",
        "mapped labelled protein rows with propagated GO annotations",
        "minimum",
    ),
    "cc_evaluable_ratio": MetricDefinition(
        "mapped labelled protein rows with at least one evaluable cellular-component term",
        "mapped labelled protein rows with propagated GO annotations",
        "minimum",
    ),
    "mf_evaluable_ratio": MetricDefinition(
        "mapped labelled protein rows with at least one evaluable molecular-function term",
        "mapped labelled protein rows with propagated GO annotations",
        "minimum",
    ),
    "development_split_deviation": MetricDefinition(
        "absolute achieved-minus-requested development member fraction",
        "one whole ratio unit",
        "maximum",
    ),
    "training_split_deviation": MetricDefinition(
        "absolute achieved-minus-requested training-within-development member fraction",
        "one whole ratio unit",
        "maximum",
    ),
}


def observation(name: str, numerator: int | float, denominator: int | float) -> dict[str, Any]:
    definition = METRIC_DEFINITIONS[name]
    numerator_value = float(numerator)
    denominator_value = float(denominator)
    if denominator_value <= 0:
        raise ValueError(f"Attrition metric {name} has a non-positive denominator")
    return {
        "name": name,
        "numerator": numerator_value,
        "denominator": denominator_value,
        "ratio": numerator_value / denominator_value,
        "numerator_definition": definition.numerator,
        "denominator_definition": definition.denominator,
        "bound_type": definition.bound,
    }


def _load_json(path: Path, label: str) -> tuple[dict[str, Any], str]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"{label} is not readable valid JSON: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{label} root must be an object")
    return payload, sha256_file(path)


def _review_date(value: object, label: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be an ISO date")
    try:
        date.fromisoformat(value)
    except ValueError as exc:
        raise ValueError(f"{label} must be an ISO date") from exc
    return value


def require_reviewed_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty string")
    normalized = value.strip()
    if normalized.upper().startswith(PLACEHOLDER_PREFIXES):
        raise ValueError(f"{label} still contains a template placeholder")
    return normalized


def load_attrition_policy(
    path: Path,
    *,
    source_scope: str,
    expected_releases: dict[str, str],
    framework_commit: str,
    frozen_input_manifest_sha256: str,
) -> tuple[dict[str, Any], str]:
    payload, digest = _load_json(path, "Attrition policy")
    if payload.get("schema_name") != POLICY_SCHEMA or payload.get("schema_version") != 1:
        raise ValueError("Attrition policy has an unsupported schema_name/schema_version")
    if payload.get("uniprot_source_scope") != source_scope:
        raise ValueError("Attrition policy is bound to the wrong UniProt source scope")
    if payload.get("expected_releases") != expected_releases:
        raise ValueError("Attrition policy is bound to the wrong frozen releases")
    if payload.get("framework_commit") != framework_commit or not COMMIT_RE.fullmatch(
        str(payload.get("framework_commit", ""))
    ):
        raise ValueError("Attrition policy is bound to the wrong framework commit")
    if (
        payload.get("frozen_input_manifest_sha256") != frozen_input_manifest_sha256
        or not SHA256_RE.fullmatch(str(payload.get("frozen_input_manifest_sha256", "")))
    ):
        raise ValueError("Attrition policy is bound to the wrong frozen-input manifest")
    for key in ("rationale", "evidence_source", "author", "reviewer"):
        require_reviewed_text(payload.get(key), f"Attrition policy field {key}")
    for key, value in payload["expected_releases"].items():
        require_reviewed_text(value, f"Attrition policy expected release {key}")
    _review_date(payload.get("review_date"), "Attrition policy review_date")
    metrics = payload.get("metrics")
    if not isinstance(metrics, dict) or set(metrics) != set(METRIC_DEFINITIONS):
        raise ValueError(
            "Attrition policy metrics must exactly match the registered production metric set"
        )
    for name, definition in METRIC_DEFINITIONS.items():
        entry = metrics[name]
        if not isinstance(entry, dict):
            raise ValueError(f"Attrition policy metric {name} must be an object")
        if entry.get("numerator_definition") != definition.numerator:
            raise ValueError(f"Attrition policy numerator definition mismatch for {name}")
        if entry.get("denominator_definition") != definition.denominator:
            raise ValueError(f"Attrition policy denominator definition mismatch for {name}")
        expected_key = f"allowed_{definition.bound}_ratio"
        other_key = "allowed_maximum_ratio" if definition.bound == "minimum" else "allowed_minimum_ratio"
        if expected_key not in entry or other_key in entry:
            raise ValueError(f"Attrition policy metric {name} has the wrong bound type")
        bound = entry[expected_key]
        if isinstance(bound, bool) or not isinstance(bound, (int, float)) or not 0 <= bound <= 1:
            raise ValueError(f"Attrition policy bound for {name} must be between 0 and 1")
        for key in ("rationale", "evidence_source"):
            require_reviewed_text(
                entry.get(key), f"Attrition policy metric {name} field {key}"
            )
    return payload, digest


def _load_override(
    path: Path,
    *,
    failures: list[dict[str, Any]],
    source_scope: str,
    framework_commit: str,
    input_manifest_sha256: str,
) -> tuple[dict[str, Any], str]:
    payload, digest = _load_json(path, "Attrition override")
    if payload.get("schema_name") != OVERRIDE_SCHEMA or payload.get("schema_version") != 1:
        raise ValueError("Attrition override has an unsupported schema_name/schema_version")
    if payload.get("uniprot_source_scope") != source_scope:
        raise ValueError("Attrition override is bound to the wrong source scope")
    if payload.get("framework_commit") != framework_commit:
        raise ValueError("Attrition override is bound to the wrong framework commit")
    if payload.get("input_manifest_sha256") != input_manifest_sha256:
        raise ValueError("Attrition override is bound to the wrong input-manifest hash")
    for key in ("justification", "reviewer", "pilot_or_run_identifier"):
        require_reviewed_text(payload.get(key), f"Attrition override field {key}")
    _review_date(payload.get("review_date"), "Attrition override review_date")
    failed_metrics = payload.get("failed_metrics")
    expected = [
        {"name": item["name"], "observed_ratio": item["observed_ratio"]}
        for item in failures
    ]
    if failed_metrics != expected:
        raise ValueError("Attrition override failed_metrics do not exactly match observed failures")
    return payload, digest


def evaluate_attrition(
    policy: dict[str, Any],
    policy_sha256: str,
    observations: dict[str, dict[str, Any]],
    *,
    source_scope: str,
    framework_commit: str,
    input_manifest_sha256: str,
    override_path: Path | None = None,
    diagnostic: bool = False,
) -> dict[str, Any]:
    if set(observations) != set(METRIC_DEFINITIONS):
        raise ValueError("Observed attrition metrics do not exactly match the registered set")
    evaluations = []
    failures = []
    for name in METRIC_DEFINITIONS:
        observed = observations[name]
        definition = METRIC_DEFINITIONS[name]
        bound_key = f"allowed_{definition.bound}_ratio"
        bound = float(policy["metrics"][name][bound_key])
        ratio = float(observed["ratio"])
        passed = ratio >= bound if definition.bound == "minimum" else ratio <= bound
        item = {
            **observed,
            "bound_key": bound_key,
            "allowed_ratio": bound,
            "passed": passed,
        }
        evaluations.append(item)
        if not passed:
            failures.append({
                "name": name,
                "observed_ratio": ratio,
                "allowed_ratio": bound,
                "bound_type": definition.bound,
            })
    override_sha256 = None
    override_valid = False
    if failures and override_path is not None:
        _, override_sha256 = _load_override(
            override_path,
            failures=failures,
            source_scope=source_scope,
            framework_commit=framework_commit,
            input_manifest_sha256=input_manifest_sha256,
        )
        override_valid = True
    passed = not failures
    return {
        "schema_name": "homology-cluster-attrition-report",
        "schema_version": 1,
        "diagnostic": diagnostic,
        "production_authorized": (not diagnostic) and (passed or override_valid),
        "policy_passed": passed,
        "override_valid": override_valid,
        "uniprot_source_scope": source_scope,
        "framework_commit": framework_commit,
        "input_manifest_sha256": input_manifest_sha256,
        "policy_sha256": policy_sha256,
        "override_sha256": override_sha256,
        "metrics": evaluations,
        "failed_metrics": failures,
    }
