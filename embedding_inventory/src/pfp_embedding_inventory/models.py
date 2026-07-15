from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Set, Tuple


MODALITIES = ("prott5", "text", "structure", "ppi")
ONTOLOGIES = ("BP", "CC", "MF")
SPLITS = ("training", "validation", "test")
Action = Literal["reuse", "regenerate"]


@dataclass
class ProteinRecord:
    protein_id: str
    sequence: str
    sequence_sha256: str
    sequence_length: int
    ontologies: Set[str] = field(default_factory=set)
    splits: Set[str] = field(default_factory=set)
    source_files: Set[str] = field(default_factory=set)
    memberships: Set[Tuple[str, str]] = field(default_factory=set)


@dataclass
class BenchmarkData:
    directory: Path
    proteins: Dict[str, ProteinRecord]
    file_members: Dict[Tuple[str, str], Set[str]]
    duplicate_rows: int = 0
    fingerprint: str = ""


@dataclass(frozen=True)
class ProvenanceSpec:
    compatibility: str
    label: str
    source_identity: str
    target_identity: str
    evidence: str
    text_role_policy: str = "none"
    requires_mapping_evidence: bool = False
    allow_direct_id_reuse: bool = False


@dataclass(frozen=True)
class ModalitySpec:
    name: str
    directory: str
    expected_dim: int
    sequence_dependent: bool
    allow_sequence_hash_reuse: bool
    provenance: ProvenanceSpec


@dataclass(frozen=True)
class BenchmarkContract:
    id_overlap: str
    sequence_overlap: str
    protein_id_pattern: str
    sequence_pattern: str


@dataclass(frozen=True)
class PlannerConfig:
    schema_version: int
    name: str
    target_benchmark_contract: BenchmarkContract
    source_benchmark_contract: BenchmarkContract
    modalities: Dict[str, ModalitySpec]
    artifact_scope: "ArtifactScopeSpec"


@dataclass(frozen=True)
class ArchiveSpec:
    path: str
    sha256: str


@dataclass(frozen=True)
class ReferenceFileSpec:
    path: str
    sha256: str


@dataclass(frozen=True)
class ArtifactScopeSpec:
    mode: str
    artifact_id: str
    metadata_url: str
    expected_benchmark_fingerprint: str
    expected_cache_catalog_fingerprint: str
    expected_modality_counts: Dict[str, int]
    expected_total_files: int
    expected_total_bytes: int
    archives: Tuple[ArchiveSpec, ...]
    expected_reference_commit: str
    reference_files: Tuple[ReferenceFileSpec, ...]


@dataclass
class CacheCatalog:
    schema: str
    fingerprint: str
    modality_fingerprints: Dict[str, str]
    modality_counts: Dict[str, int]
    modality_bytes: Dict[str, int]
    total_files: int
    total_bytes: int

    def as_dict(self) -> Dict[str, Any]:
        return {
            "schema": self.schema,
            "fingerprint": self.fingerprint,
            "modality_fingerprints": dict(self.modality_fingerprints),
            "modality_counts": dict(self.modality_counts),
            "modality_bytes": dict(self.modality_bytes),
            "total_files": self.total_files,
            "total_bytes": self.total_bytes,
        }


@dataclass
class ArtifactVerification:
    configured: bool
    verified: bool
    artifact_id: str
    checks: Dict[str, bool]
    reasons: List[str]
    expected: Dict[str, Any]
    observed: Dict[str, Any]

    def as_dict(self) -> Dict[str, Any]:
        return {
            "configured": self.configured,
            "verified": self.verified,
            "artifact_id": self.artifact_id,
            "checks": dict(self.checks),
            "reasons": list(self.reasons),
            "expected": self.expected,
            "observed": self.observed,
        }


@dataclass(frozen=True)
class AliasEntry:
    protein_id: str
    source_protein_id: str
    modality: str
    mapping_route: str
    source_identity: str
    mapping_evidence: str


@dataclass
class ArrayInfo:
    exists: bool
    observed_shape: str = ""
    dtype: str = ""
    finite: Optional[bool] = None
    dtype_supported: Optional[bool] = None
    valid: bool = False
    error: str = ""


@dataclass(frozen=True)
class InventoryRecord:
    protein_id: str
    sequence_sha256: str
    modality: str
    source_directory: str
    source_file: str
    exists: bool
    observed_shape: str
    expected_shape: str
    dtype: str
    finite: Optional[bool]
    valid: bool
    scientifically_eligible: bool
    source_protein_id: str
    match_route: str
    sequence_match: str
    provenance: str
    factual_status: str
    requested_action: Action
    reason: str

    def __post_init__(self) -> None:
        if self.requested_action not in {"reuse", "regenerate"}:
            raise ValueError(
                "requested_action must be exactly 'reuse' or 'regenerate'"
            )

    def as_dict(self) -> Dict[str, str]:
        return {
            "protein_id": self.protein_id,
            "sequence_sha256": self.sequence_sha256,
            "modality": self.modality,
            "source_directory": self.source_directory,
            "source_file": self.source_file,
            "exists": _bool_text(self.exists),
            "observed_shape": self.observed_shape,
            "expected_shape": self.expected_shape,
            "dtype": self.dtype,
            "finite": "" if self.finite is None else _bool_text(self.finite),
            "valid": _bool_text(self.valid),
            "scientifically_eligible": _bool_text(self.scientifically_eligible),
            "source_protein_id": self.source_protein_id,
            "match_route": self.match_route,
            "sequence_match": self.sequence_match,
            "provenance": self.provenance,
            "factual_status": self.factual_status,
            "requested_action": self.requested_action,
            "reason": self.reason,
        }


@dataclass
class InventoryResult:
    benchmark: BenchmarkData
    source_benchmark: BenchmarkData
    records: List[InventoryRecord]
    cache_ids: Dict[str, Set[str]]
    used_source_ids: Dict[str, Set[str]]
    policy: str
    config: PlannerConfig
    artifact_verification: ArtifactVerification


def _bool_text(value: bool) -> str:
    return "true" if value else "false"
