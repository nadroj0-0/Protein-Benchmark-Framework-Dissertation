from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Tuple


ONTOLOGIES: Tuple[str, ...] = ("bp", "cc", "mf")
SPLITS: Tuple[str, ...] = ("training", "validation", "test")
REQUIRED_CSV_NAMES: Tuple[str, ...] = tuple(
    f"{ontology}-{split}.csv" for ontology in ONTOLOGIES for split in SPLITS
)
REGENERATE_MODALITIES: Tuple[str, ...] = tuple(
    sorted(("prott5", "text", "structure", "ppi"))
)


class ReusePlannerError(ValueError):
    """Base class for expected input, planning, and reporting failures."""


@dataclass(frozen=True)
class InputFileIdentity:
    relative_path: str
    resolved_path: Path
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class ProteinRecord:
    protein_id: str
    sequence: str
    sequence_sha256: str
    memberships: Tuple[str, ...]


@dataclass(frozen=True)
class BenchmarkData:
    name: str
    directory: Path
    proteins: Mapping[str, ProteinRecord]
    input_files: Tuple[InputFileIdentity, ...]
    duplicate_occurrences: int


@dataclass(frozen=True)
class EmbeddedProtein:
    protein_id: str
    sequence: str
    sequence_sha256: str
    embedded_benchmarks: Tuple[str, ...]
    embedded_benchmark_memberships: Tuple[str, ...]


@dataclass(frozen=True)
class PlanRecord:
    protein_id: str
    sequence: str
    sequence_sha256: str
    action: str
    reason: str
    matching_embedded_benchmarks: Tuple[str, ...]
    embedded_benchmark_memberships: Tuple[str, ...]
    target_memberships: Tuple[str, ...]
    regenerate_modalities: Tuple[str, ...]


@dataclass(frozen=True)
class ReusePlan:
    embedded_benchmarks: Tuple[BenchmarkData, ...]
    target_benchmark: BenchmarkData
    known_embedded_proteins: Tuple[EmbeddedProtein, ...]
    records: Tuple[PlanRecord, ...]

    @property
    def reuse_records(self) -> Tuple[PlanRecord, ...]:
        return tuple(record for record in self.records if record.action == "reuse")

    @property
    def regenerate_records(self) -> Tuple[PlanRecord, ...]:
        return tuple(record for record in self.records if record.action == "regenerate")
