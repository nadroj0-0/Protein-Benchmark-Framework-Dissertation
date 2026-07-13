from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


MODALITIES = ("prott5", "text", "structure", "ppi")
ONTOLOGIES = ("BP", "CC", "MF")
SPLITS = ("training", "validation", "test")


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


@dataclass(frozen=True)
class ProvenanceSpec:
    compatibility: str
    label: str
    source_identity: str
    target_identity: str
    evidence: str
    text_role_policy: str = "none"
    requires_mapping_evidence: bool = False


@dataclass(frozen=True)
class ModalitySpec:
    name: str
    directory: str
    expected_dim: int
    sequence_dependent: bool
    allow_sequence_hash_reuse: bool
    missing_action: str
    invalid_action: str
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
    benchmark_contract: BenchmarkContract
    modalities: Dict[str, ModalitySpec]


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


@dataclass
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
    requested_action: str
    reason: str

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


def _bool_text(value: bool) -> str:
    return "true" if value else "false"
