import json
import re
from pathlib import Path
from typing import Any, Dict

from .models import (
    BenchmarkContract,
    MODALITIES,
    ModalitySpec,
    PlannerConfig,
    ProvenanceSpec,
)


class ConfigError(ValueError):
    pass


_OVERLAP_POLICIES = {"allow", "global-disjoint", "per-ontology-disjoint"}
_COMPATIBILITY = {"compatible", "artifact-scoped", "unknown", "incompatible"}
_PLAN_ACTIONS = {"generate", "unavailable", "leave-masked", "manual-review"}
_EXPECTED_DIMS = {"prott5": 1024, "text": 768, "structure": 512, "ppi": 512}


def load_config(path: Path) -> PlannerConfig:
    try:
        raw = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigError("Cannot read configuration %s: %s" % (path, exc)) from exc

    if raw.get("schema_version") != 1:
        raise ConfigError("schema_version must be 1")

    contract_raw = _mapping(raw, "benchmark_contract")
    id_overlap = str(contract_raw.get("id_overlap", "per-ontology-disjoint"))
    seq_overlap = str(contract_raw.get("sequence_overlap", "per-ontology-disjoint"))
    if id_overlap not in _OVERLAP_POLICIES or seq_overlap not in _OVERLAP_POLICIES:
        raise ConfigError("Unsupported benchmark overlap policy")
    id_pattern = str(contract_raw.get("protein_id_pattern", r"^[^\s/\\]+$"))
    seq_pattern = str(contract_raw.get("sequence_pattern", r"^[A-Za-z*.-]+$"))
    try:
        re.compile(id_pattern)
        re.compile(seq_pattern)
    except re.error as exc:
        raise ConfigError("Invalid benchmark validation regex: %s" % exc) from exc

    modalities_raw = _mapping(raw, "modalities")
    if set(modalities_raw) != set(MODALITIES):
        raise ConfigError("modalities must be exactly: %s" % ", ".join(MODALITIES))

    modalities: Dict[str, ModalitySpec] = {}
    for name in MODALITIES:
        item = _mapping(modalities_raw, name)
        directory_value = item.get("directory", "")
        if isinstance(directory_value, list):
            raise ConfigError("modality %s directory must be singular, not a fallback list" % name)
        directory = str(directory_value).strip()
        if not directory:
            raise ConfigError("modality %s requires one explicit directory" % name)
        directory_path = Path(directory)
        if directory_path.is_absolute() or ".." in directory_path.parts or directory in {".", ".."}:
            raise ConfigError(
                "modality %s directory must be a cache-relative path without traversal" % name
            )
        expected_dim = item.get("expected_dim")
        if type(expected_dim) is not int or expected_dim != _EXPECTED_DIMS[name]:
            raise ConfigError(
                "modality %s expected_dim must be the PFP dimension %d"
                % (name, _EXPECTED_DIMS[name])
            )
        missing_action = str(item.get("missing_action", "manual-review"))
        invalid_action = str(item.get("invalid_action", "manual-review"))
        if missing_action not in _PLAN_ACTIONS or invalid_action not in _PLAN_ACTIONS:
            raise ConfigError("modality %s has an unsupported planning action" % name)

        provenance_raw = _mapping(item, "provenance")
        compatibility = str(provenance_raw.get("compatibility", "unknown"))
        if compatibility not in _COMPATIBILITY:
            raise ConfigError("modality %s has invalid provenance compatibility" % name)
        source_identity = str(provenance_raw.get("source_identity", "")).strip()
        target_identity = str(provenance_raw.get("target_identity", "")).strip()
        evidence = str(provenance_raw.get("evidence", "")).strip()
        label = str(provenance_raw.get("label", compatibility)).strip()
        requires_mapping_evidence = _required_bool(
            provenance_raw, "requires_mapping_evidence", name
        )
        expected_mapping_evidence = name != "prott5"
        if requires_mapping_evidence != expected_mapping_evidence:
            raise ConfigError(
                "modality %s requires_mapping_evidence must be %s"
                % (name, str(expected_mapping_evidence).lower())
            )
        if compatibility in {"compatible", "artifact-scoped"}:
            if not source_identity or not target_identity or not evidence:
                raise ConfigError(
                    "modality %s compatible provenance requires source_identity, "
                    "target_identity, and evidence" % name
                )
            if source_identity != target_identity:
                raise ConfigError(
                    "modality %s cannot be compatible when source and target identities differ" % name
                )

        sequence_dependent = _required_bool(item, "sequence_dependent", name)
        allow_sequence_hash_reuse = _required_bool(
            item, "allow_sequence_hash_reuse", name
        )
        expected_sequence_dependent = name in {"prott5", "structure"}
        expected_hash_reuse = name == "prott5"
        if sequence_dependent != expected_sequence_dependent:
            raise ConfigError(
                "modality %s sequence_dependent must be %s"
                % (name, str(expected_sequence_dependent).lower())
            )
        if allow_sequence_hash_reuse != expected_hash_reuse:
            raise ConfigError(
                "modality %s allow_sequence_hash_reuse must be %s"
                % (name, str(expected_hash_reuse).lower())
            )

        modalities[name] = ModalitySpec(
            name=name,
            directory=directory,
            expected_dim=expected_dim,
            sequence_dependent=sequence_dependent,
            allow_sequence_hash_reuse=allow_sequence_hash_reuse,
            missing_action=missing_action,
            invalid_action=invalid_action,
            provenance=ProvenanceSpec(
                compatibility=compatibility,
                label=label,
                source_identity=source_identity,
                target_identity=target_identity,
                evidence=evidence,
                text_role_policy=str(provenance_raw.get("text_role_policy", "none")),
                requires_mapping_evidence=requires_mapping_evidence,
            ),
        )

    return PlannerConfig(
        schema_version=1,
        name=str(raw.get("name", path.stem)),
        benchmark_contract=BenchmarkContract(
            id_overlap=id_overlap,
            sequence_overlap=seq_overlap,
            protein_id_pattern=id_pattern,
            sequence_pattern=seq_pattern,
        ),
        modalities=modalities,
    )


def _mapping(parent: Dict[str, Any], key: str) -> Dict[str, Any]:
    value = parent.get(key)
    if not isinstance(value, dict):
        raise ConfigError("%s must be an object" % key)
    return value


def _required_bool(parent: Dict[str, Any], key: str, modality: str) -> bool:
    value = parent.get(key)
    if type(value) is not bool:
        raise ConfigError("modality %s %s must be a JSON boolean" % (modality, key))
    return value
