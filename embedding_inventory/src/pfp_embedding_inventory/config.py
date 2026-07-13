import json
import re
from pathlib import Path
from typing import Any, Dict

from .models import (
    ArchiveSpec,
    ArtifactScopeSpec,
    BenchmarkContract,
    MODALITIES,
    ModalitySpec,
    PlannerConfig,
    ProvenanceSpec,
    ReferenceFileSpec,
)


class ConfigError(ValueError):
    pass


_OVERLAP_POLICIES = {
    "allow",
    "global-disjoint",
    "global-evaluation-disjoint",
    "per-ontology-disjoint",
}
_COMPATIBILITY = {"compatible", "artifact-scoped", "unknown", "incompatible"}
_PLAN_ACTIONS = {"generate", "unavailable", "leave-masked", "manual-review"}
_EXPECTED_DIMS = {"prott5": 1024, "text": 768, "structure": 512, "ppi": 512}
_ARTIFACT_MODES = {"none", "verified-published-cache"}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_GIT_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")


def load_config(path: Path) -> PlannerConfig:
    try:
        raw = json.loads(path.read_text())
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ConfigError("Cannot read configuration %s: %s" % (path, exc)) from exc

    if not isinstance(raw, dict):
        raise ConfigError("configuration root must be a JSON object")
    version = raw.get("schema_version")
    if version != 2:
        if version == 1:
            raise ConfigError(
                "schema_version 1 is not accepted: migrate benchmark_contract to explicit "
                "target_benchmark_contract and source_benchmark_contract"
            )
        raise ConfigError("schema_version must be 2")
    if "benchmark_contract" in raw:
        raise ConfigError("schema v2 forbids the ambiguous benchmark_contract key")
    _exact_keys(
        raw,
        {
            "schema_version",
            "name",
            "target_benchmark_contract",
            "source_benchmark_contract",
            "artifact_scope",
            "modalities",
        },
        "top-level configuration",
    )

    target_contract = _load_contract(_mapping(raw, "target_benchmark_contract"), "target")
    source_contract = _load_contract(_mapping(raw, "source_benchmark_contract"), "source")
    artifact_scope = _load_artifact_scope(_mapping(raw, "artifact_scope"))

    modalities_raw = _mapping(raw, "modalities")
    if set(modalities_raw) != set(MODALITIES):
        raise ConfigError("modalities must be exactly: %s" % ", ".join(MODALITIES))

    modalities: Dict[str, ModalitySpec] = {}
    for name in MODALITIES:
        item = _mapping(modalities_raw, name)
        _exact_keys(
            item,
            {
                "directory",
                "expected_dim",
                "sequence_dependent",
                "allow_sequence_hash_reuse",
                "missing_action",
                "invalid_action",
                "provenance",
            },
            "modality %s" % name,
        )
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
        _allowed_keys(
            provenance_raw,
            {
                "compatibility",
                "label",
                "source_identity",
                "target_identity",
                "evidence",
                "text_role_policy",
                "requires_mapping_evidence",
            },
            "modality %s provenance" % name,
        )
        _require_keys(
            provenance_raw,
            {
                "compatibility",
                "label",
                "source_identity",
                "target_identity",
                "evidence",
                "requires_mapping_evidence",
            },
            "modality %s provenance" % name,
        )
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
        if compatibility == "artifact-scoped" and artifact_scope.mode == "none":
            raise ConfigError(
                "modality %s cannot be artifact-scoped without verified artifact metadata" % name
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

    directories = [modalities[name].directory for name in MODALITIES]
    if len(directories) != len(set(directories)):
        raise ConfigError("each modality must use a distinct cache directory")

    return PlannerConfig(
        schema_version=2,
        name=str(raw.get("name", path.stem)),
        target_benchmark_contract=target_contract,
        source_benchmark_contract=source_contract,
        modalities=modalities,
        artifact_scope=artifact_scope,
    )


def _load_contract(raw: Dict[str, Any], label: str) -> BenchmarkContract:
    _exact_keys(
        raw,
        {"id_overlap", "sequence_overlap", "protein_id_pattern", "sequence_pattern"},
        "%s benchmark contract" % label,
    )
    id_overlap = str(raw["id_overlap"])
    sequence_overlap = str(raw["sequence_overlap"])
    if id_overlap not in _OVERLAP_POLICIES or sequence_overlap not in _OVERLAP_POLICIES:
        raise ConfigError("Unsupported %s benchmark overlap policy" % label)
    id_pattern = str(raw["protein_id_pattern"])
    sequence_pattern = str(raw["sequence_pattern"])
    try:
        re.compile(id_pattern)
        re.compile(sequence_pattern)
    except re.error as exc:
        raise ConfigError("Invalid %s benchmark validation regex: %s" % (label, exc)) from exc
    return BenchmarkContract(
        id_overlap=id_overlap,
        sequence_overlap=sequence_overlap,
        protein_id_pattern=id_pattern,
        sequence_pattern=sequence_pattern,
    )


def _load_artifact_scope(raw: Dict[str, Any]) -> ArtifactScopeSpec:
    mode = str(raw.get("mode", "none"))
    if mode not in _ARTIFACT_MODES:
        raise ConfigError("Unsupported artifact_scope mode: %s" % mode)
    if mode == "none":
        _exact_keys(raw, {"mode"}, "artifact_scope")
        return ArtifactScopeSpec(
            mode=mode,
            artifact_id=str(raw.get("artifact_id", "")),
            metadata_url=str(raw.get("metadata_url", "")),
            expected_benchmark_fingerprint="",
            expected_cache_catalog_fingerprint="",
            expected_modality_counts={},
            expected_total_files=0,
            expected_total_bytes=0,
            archives=(),
            expected_reference_commit="",
            reference_files=(),
        )

    _exact_keys(
        raw,
        {
            "mode", "artifact_id", "metadata_url",
            "expected_benchmark_fingerprint", "expected_cache_catalog_fingerprint",
            "expected_modality_counts", "expected_total_files", "expected_total_bytes",
            "archives", "expected_reference_commit", "reference_files",
        },
        "artifact_scope",
    )
    artifact_id = str(raw.get("artifact_id", "")).strip()
    metadata_url = str(raw.get("metadata_url", "")).strip()
    benchmark_fingerprint = _sha256_value(raw, "expected_benchmark_fingerprint")
    cache_fingerprint = _sha256_value(raw, "expected_cache_catalog_fingerprint")
    counts_raw = _mapping(raw, "expected_modality_counts")
    if set(counts_raw) != set(MODALITIES):
        raise ConfigError("artifact expected_modality_counts must cover all modalities")
    counts: Dict[str, int] = {}
    for modality in MODALITIES:
        value = counts_raw[modality]
        if type(value) is not int or value < 0:
            raise ConfigError("artifact modality counts must be non-negative integers")
        counts[modality] = value
    total_files = raw.get("expected_total_files")
    total_bytes = raw.get("expected_total_bytes")
    if type(total_files) is not int or total_files != sum(counts.values()):
        raise ConfigError("artifact expected_total_files must equal modality counts")
    if type(total_bytes) is not int or total_bytes <= 0:
        raise ConfigError("artifact expected_total_bytes must be a positive integer")
    if not artifact_id or not metadata_url:
        raise ConfigError("verified artifact scope requires artifact_id and metadata_url")

    archives_raw = raw.get("archives")
    if not isinstance(archives_raw, list) or not archives_raw:
        raise ConfigError("verified artifact scope requires archive checksums")
    archives = []
    archive_paths = set()
    for item in archives_raw:
        if not isinstance(item, dict):
            raise ConfigError("artifact archive entries must be objects")
        _exact_keys(item, {"path", "sha256"}, "artifact archive entry")
        archive_path = str(item.get("path", "")).strip()
        if not archive_path or Path(archive_path).is_absolute() or ".." in Path(archive_path).parts:
            raise ConfigError("artifact archive paths must be safe relative paths")
        if archive_path in archive_paths:
            raise ConfigError("artifact archive paths must be unique")
        archive_paths.add(archive_path)
        archives.append(ArchiveSpec(path=archive_path, sha256=_sha256_value(item, "sha256")))

    reference_commit = str(raw.get("expected_reference_commit", "")).strip().lower()
    if _GIT_COMMIT_RE.fullmatch(reference_commit) is None:
        raise ConfigError(
            "verified artifact scope requires expected_reference_commit as a 40-character commit"
        )
    reference_files_raw = raw.get("reference_files")
    if not isinstance(reference_files_raw, list) or not reference_files_raw:
        raise ConfigError("verified artifact scope requires reference file checksums")
    reference_files = []
    reference_paths = set()
    for item in reference_files_raw:
        if not isinstance(item, dict):
            raise ConfigError("artifact reference file entries must be objects")
        _exact_keys(item, {"path", "sha256"}, "artifact reference file entry")
        reference_path = str(item.get("path", "")).strip()
        if (
            not reference_path
            or Path(reference_path).is_absolute()
            or ".." in Path(reference_path).parts
        ):
            raise ConfigError("artifact reference file paths must be safe relative paths")
        if reference_path in reference_paths:
            raise ConfigError("artifact reference file paths must be unique")
        reference_paths.add(reference_path)
        reference_files.append(
            ReferenceFileSpec(path=reference_path, sha256=_sha256_value(item, "sha256"))
        )
    return ArtifactScopeSpec(
        mode=mode,
        artifact_id=artifact_id,
        metadata_url=metadata_url,
        expected_benchmark_fingerprint=benchmark_fingerprint,
        expected_cache_catalog_fingerprint=cache_fingerprint,
        expected_modality_counts=counts,
        expected_total_files=total_files,
        expected_total_bytes=total_bytes,
        archives=tuple(archives),
        expected_reference_commit=reference_commit,
        reference_files=tuple(reference_files),
    )


def _sha256_value(parent: Dict[str, Any], key: str) -> str:
    value = str(parent.get(key, "")).lower()
    if _SHA256_RE.fullmatch(value) is None:
        raise ConfigError("%s must be a lowercase SHA-256 value" % key)
    return value


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


def _exact_keys(parent: Dict[str, Any], expected: set, label: str) -> None:
    _require_keys(parent, expected, label)
    _allowed_keys(parent, expected, label)


def _require_keys(parent: Dict[str, Any], required: set, label: str) -> None:
    missing = sorted(required - set(parent))
    if missing:
        raise ConfigError("%s is missing required keys: %s" % (label, ", ".join(missing)))


def _allowed_keys(parent: Dict[str, Any], allowed: set, label: str) -> None:
    unknown = sorted(set(parent) - allowed)
    if unknown:
        raise ConfigError("%s contains unknown keys: %s" % (label, ", ".join(unknown)))
