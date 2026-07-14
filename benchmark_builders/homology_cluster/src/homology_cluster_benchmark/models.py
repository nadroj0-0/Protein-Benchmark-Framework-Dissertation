from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class InputSpec:
    """A reproducible local-or-remote input declaration."""

    name: str
    path: Path | None = None
    url: str | None = None
    expected_sha256: str | None = None
    release: str = ""


@dataclass(frozen=True)
class ResolvedInput:
    name: str
    resolved_path: Path
    source_url: str | None
    release: str
    size_bytes: int
    sha256: str
    expected_sha256: str | None
    acquisition: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "resolved_path": str(self.resolved_path.resolve()),
            "source_url": self.source_url,
            "release": self.release,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
            "expected_sha256": self.expected_sha256,
            "acquisition": self.acquisition,
        }


@dataclass(frozen=True)
class ProteinRecord:
    protein_id: str
    sequence: str
    taxon_id: str = ""
    aliases: tuple[str, ...] = ()
    source_header: str = ""


@dataclass
class ProteinCatalog:
    records: dict[str, ProteinRecord] = field(default_factory=dict)
    alias_to_primary: dict[str, str] = field(default_factory=dict)
    ambiguous_aliases: set[str] = field(default_factory=set)


@dataclass(frozen=True)
class AnnotationRecord:
    database: str
    raw_accession: str
    protein_id: str
    symbol: str
    raw_go_id: str
    go_id: str
    namespace: str
    aspect: str
    evidence: str
    qualifier: str
    reference: str
    with_from: str
    taxon_id: str
    assigned_date: str
    assigned_by: str
    annotation_extension: str
    gene_product_form: str
    line_number: int
    term_action: str
    accession_action: str = "exact"


@dataclass(frozen=True)
class ExcludedAnnotation:
    raw_accession: str
    raw_go_id: str
    evidence: str
    qualifier: str
    aspect: str
    object_type: str
    taxon_id: str
    line_number: int
    rejection_reason: str
    detail: str = ""
    database: str = ""
    symbol: str = ""
    reference: str = ""
    with_from: str = ""
    assigned_date: str = ""
    assigned_by: str = ""
    annotation_extension: str = ""
    gene_product_form: str = ""


@dataclass
class GoaLoadResult:
    records: list[AnnotationRecord]
    excluded: list[ExcludedAnnotation]
    annotations: dict[str, set[str]]
    qualifying_accessions: set[str] = field(default_factory=set)
    counters: Counter = field(default_factory=Counter)
    evidence_counts: Counter = field(default_factory=Counter)
    taxonomy_counts: Counter = field(default_factory=Counter)
    direct_go_counts: Counter = field(default_factory=Counter)
    rejected_evidence_reason_counts: Counter = field(default_factory=Counter)
    rejection_reason_counts: Counter = field(default_factory=Counter)
    rejection_evidence_counts: Counter = field(default_factory=Counter)
    rejection_database_counts: Counter = field(default_factory=Counter)
    rejection_object_type_counts: Counter = field(default_factory=Counter)
    rejection_aspect_counts: Counter = field(default_factory=Counter)
    accepted_database_counts: Counter = field(default_factory=Counter)
    accepted_object_type_counts: Counter = field(default_factory=Counter)
    accepted_aspect_counts: Counter = field(default_factory=Counter)
    annotation_decision_counts: Counter = field(default_factory=Counter)
    excluded_sample_counts: Counter = field(default_factory=Counter)
    candidate_accessions: set[str] = field(default_factory=set)
    headers: dict[str, str] = field(default_factory=dict)
    record_spool: Path | None = None
    excluded_spool: Path | None = None
    accession_map: dict[str, tuple[str, str]] = field(default_factory=dict)


@dataclass
class MappingDecision:
    raw_accession: str
    protein_id: str
    accession_action: str
    uniref90_id: str = ""
    status: str = "unmapped"
    detail: str = ""
    # None means not applicable/ambiguous; true/false is asserted only for one selected ID.
    exists_in_fasta: bool | None = None
    canonical_sequence_available: bool = False
    mmseqs_cluster_id: str = ""
    split: str = ""
    accession_lifecycle_status: str = "not-assessed"


@dataclass(frozen=True)
class ClusterInfo:
    cluster_id: str
    member_count: int
    labelled_protein_count: int = 0


@dataclass(frozen=True)
class SplitAssignment:
    cluster_id: str
    split: str
    member_count: int
    labelled_protein_count: int
    stage: str


@dataclass
class LabelBuildResult:
    frames: dict[str, Any]
    unrestricted_annotations: dict[str, tuple[str, ...]]
    restricted_annotations: dict[str, tuple[str, ...]]
    term_universe: tuple[str, ...]
    removed_term_counts: Counter
    no_evaluable_term: Counter
    annotation_exclusion_counts: Counter = field(default_factory=Counter)
    row_attrition_counts: Counter = field(default_factory=Counter)
    protein_attrition_counts: Counter = field(default_factory=Counter)
    intended_annotation_rows: int = 0
    intended_accessions: int = 0
    cluster_assignments: dict[str, SplitAssignment] = field(default_factory=dict)


@dataclass
class ValidationReport:
    checks: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[dict[str, Any]] = field(default_factory=list)

    def add_check(self, name: str, passed: bool, detail: str, **metrics: Any) -> None:
        self.checks.append({"name": name, "passed": bool(passed), "detail": detail, **metrics})

    def add_warning(self, name: str, detail: str, **metrics: Any) -> None:
        self.warnings.append({"name": name, "detail": detail, **metrics})

    @property
    def valid(self) -> bool:
        return all(item["passed"] for item in self.checks)


@dataclass(frozen=True)
class BuildResult:
    output_dir: Path
    identity: float
    files: tuple[Path, ...]
    valid: bool
