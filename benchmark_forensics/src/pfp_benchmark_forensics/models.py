"""Internal records for benchmark forensic analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Tuple


@dataclass(frozen=True)
class Observation:
    dataset_id: str
    aspect: str
    split: str
    protein_id: str
    sequence_sha256: str
    sequence_length: int
    label_count: int
    non_root_label_count: int
    root_positive: bool
    root_only: bool
    all_zero: bool
    positive_terms: Tuple[str, ...]


@dataclass(frozen=True)
class SourceRecord:
    protein_id: str
    sequence: str
    annotations: Tuple[str, ...]


@dataclass(frozen=True)
class Taxon:
    taxon_id: str
    taxon_name: str
    source: str


@dataclass
class DatasetResult:
    dataset_id: str
    config: object
    observations: Tuple[Observation, ...]
    file_profiles: Tuple[dict, ...]
    term_headers: Mapping[str, Tuple[str, ...]]
    sequences: Mapping[str, str]
    source_by_split: Mapping[str, Mapping[str, SourceRecord]]
    taxonomy: Mapping[str, Taxon]
    modality_states: Mapping[Tuple[str, str], Mapping[str, bool]]
    modalities: Tuple[str, ...]
    category_maps: Mapping[str, Mapping[str, Tuple[Tuple[str, str], ...]]]
    input_paths: Tuple[Path, ...]
    diagnostics: Dict[str, object]


@dataclass(frozen=True)
class AnalysisBundle:
    summary: dict
    tables: Mapping[str, Tuple[dict, ...]]
    input_paths: Tuple[Path, ...]
