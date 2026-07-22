"""Strict configuration loading for benchmark forensics."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple


class ConfigError(ValueError):
    """Raised when a forensics configuration is ambiguous or invalid."""


@dataclass(frozen=True)
class SourceAnnotationSpec:
    type: str
    path: Path
    split_files: Mapping[str, str]
    protein_id_column: str
    sequence_column: str
    annotations_column: str
    split_column: str
    annotation_separator: str
    projection_policy: str


@dataclass(frozen=True)
class TaxonomySourceSpec:
    type: str
    path: Path
    id_columns: Tuple[str, ...]
    taxon_id_column: str
    taxon_name_column: str


@dataclass(frozen=True)
class ModalityStateSpec:
    column: str
    true_values: Tuple[str, ...]


@dataclass(frozen=True)
class ModalityInventorySpec:
    type: str
    path: Path
    protein_id_column: str
    modality_column: str
    states: Mapping[str, ModalityStateSpec]


@dataclass(frozen=True)
class CategorySourceSpec:
    name: str
    path: Path
    protein_id_column: str
    category_id_column: str
    category_name_column: str


@dataclass(frozen=True)
class DatasetConfig:
    id: str
    benchmark_dir: Path
    obo_file: Path
    allow_legacy_singular_protein_header: bool
    allow_all_zero_rows: bool
    split_overlap_policy: str
    source_annotations: Optional[SourceAnnotationSpec]
    taxonomy_sources: Tuple[TaxonomySourceSpec, ...]
    modality_inventory: Optional[ModalityInventorySpec]
    category_sources: Tuple[CategorySourceSpec, ...]


@dataclass(frozen=True)
class RunConfig:
    schema_version: int
    run_name: str
    top_n: int
    datasets: Tuple[DatasetConfig, ...]
    source_path: Path


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, dict):
        raise ConfigError(f"{label} must be a JSON object")
    return value


def _keys(
    value: Mapping[str, Any],
    label: str,
    *,
    required: Sequence[str],
    optional: Sequence[str] = (),
) -> None:
    missing = sorted(set(required) - set(value))
    unknown = sorted(set(value) - set(required) - set(optional))
    if missing:
        raise ConfigError(f"{label} is missing required keys: {', '.join(missing)}")
    if unknown:
        raise ConfigError(f"{label} contains unsupported keys: {', '.join(unknown)}")


def _text(value: Any, label: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str) or (not allow_empty and not value.strip()):
        raise ConfigError(f"{label} must be a non-empty string")
    return value


def _boolean(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigError(f"{label} must be a JSON boolean")
    return value


def _path(value: Any, label: str, base: Path) -> Path:
    text = _text(value, label)
    path = Path(text).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def _string_tuple(value: Any, label: str) -> Tuple[str, ...]:
    if not isinstance(value, list) or not value:
        raise ConfigError(f"{label} must be a non-empty JSON array")
    result = tuple(_text(item, f"{label}[]") for item in value)
    if len(result) != len(set(result)):
        raise ConfigError(f"{label} contains duplicates")
    return result


def _source_annotations(raw: Any, label: str, base: Path) -> SourceAnnotationSpec:
    value = _mapping(raw, label)
    _keys(
        value,
        label,
        required=("type", "path"),
        optional=(
            "split_files",
            "protein_id_column",
            "sequence_column",
            "annotations_column",
            "split_column",
            "annotation_separator",
            "projection_policy",
        ),
    )
    source_type = _text(value["type"], f"{label}.type")
    if source_type not in {"pfp-pickle-directory", "long-tsv"}:
        raise ConfigError(f"{label}.type must be pfp-pickle-directory or long-tsv")
    split_files_raw = value.get(
        "split_files",
        {
            "training": "train_data_train.pkl",
            "validation": "train_data_valid.pkl",
            "test": "test_data.pkl",
        },
    )
    split_files = _mapping(split_files_raw, f"{label}.split_files")
    if source_type == "pfp-pickle-directory":
        _keys(
            split_files,
            f"{label}.split_files",
            required=("training", "validation", "test"),
        )
        split_files = {
            split: _text(name, f"{label}.split_files.{split}")
            for split, name in split_files.items()
        }
        for split, name in split_files.items():
            if Path(name).is_absolute() or Path(name).name != name:
                raise ConfigError(
                    f"{label}.split_files.{split} must be one filename within the source directory"
                )
    elif "split_files" in value:
        raise ConfigError(f"{label}.split_files is only valid for pickle directories")
    return SourceAnnotationSpec(
        type=source_type,
        path=_path(value["path"], f"{label}.path", base),
        split_files=split_files,
        protein_id_column=_text(
            value.get("protein_id_column", "proteins"),
            f"{label}.protein_id_column",
        ),
        sequence_column=_text(
            value.get("sequence_column", "sequences"),
            f"{label}.sequence_column",
        ),
        annotations_column=_text(
            value.get("annotations_column", "annotations"),
            f"{label}.annotations_column",
        ),
        split_column=_text(value.get("split_column", "split"), f"{label}.split_column"),
        annotation_separator=_text(
            value.get("annotation_separator", ";"),
            f"{label}.annotation_separator",
        ),
        projection_policy=_text(
            value.get("projection_policy", "unspecified-label-universe-projection"),
            f"{label}.projection_policy",
        ),
    )


def _taxonomy_source(raw: Any, label: str, base: Path) -> TaxonomySourceSpec:
    value = _mapping(raw, label)
    _keys(
        value,
        label,
        required=("type", "path"),
        optional=("id_columns", "taxon_id_column", "taxon_name_column"),
    )
    source_type = _text(value["type"], f"{label}.type")
    if source_type not in {"tsv", "uniprot-dat"}:
        raise ConfigError(f"{label}.type must be tsv or uniprot-dat")
    if source_type == "uniprot-dat" and any(
        key in value for key in ("id_columns", "taxon_id_column", "taxon_name_column")
    ):
        raise ConfigError(f"{label} column settings are only valid for TSV sources")
    return TaxonomySourceSpec(
        type=source_type,
        path=_path(value["path"], f"{label}.path", base),
        id_columns=_string_tuple(
            value.get("id_columns", ["protein_id"]), f"{label}.id_columns"
        ),
        taxon_id_column=_text(
            value.get("taxon_id_column", "taxon_id"),
            f"{label}.taxon_id_column",
        ),
        taxon_name_column=_text(
            value.get("taxon_name_column", "taxon_name"),
            f"{label}.taxon_name_column",
            allow_empty=True,
        ),
    )


def _modality_inventory(raw: Any, label: str, base: Path) -> ModalityInventorySpec:
    value = _mapping(raw, label)
    _keys(
        value,
        label,
        required=("type", "path"),
        optional=("protein_id_column", "modality_column", "states"),
    )
    inventory_type = _text(value["type"], f"{label}.type")
    if inventory_type not in {"embedding-inventory", "long-table"}:
        raise ConfigError(f"{label}.type must be embedding-inventory or long-table")
    if inventory_type == "embedding-inventory":
        if "states" in value:
            raise ConfigError(
                f"{label}.states is fixed for the embedding-inventory format"
            )
        states_raw: Mapping[str, Any] = {
            "artifact_exists": {"column": "exists", "true_values": ["true"]},
            "artifact_valid": {"column": "valid", "true_values": ["true"]},
            "scientifically_eligible": {
                "column": "scientifically_eligible",
                "true_values": ["true"],
            },
            "planned_reuse": {"column": "requested_action", "true_values": ["reuse"]},
        }
    else:
        states_raw = _mapping(value.get("states"), f"{label}.states")
        if not states_raw:
            raise ConfigError(f"{label}.states must not be empty")
    states = {}
    for state_name, state_raw in states_raw.items():
        _text(state_name, f"{label}.states key")
        state = _mapping(state_raw, f"{label}.states.{state_name}")
        _keys(
            state,
            f"{label}.states.{state_name}",
            required=("column", "true_values"),
        )
        states[state_name] = ModalityStateSpec(
            column=_text(state["column"], f"{label}.states.{state_name}.column"),
            true_values=tuple(
                item.casefold()
                for item in _string_tuple(
                    state["true_values"],
                    f"{label}.states.{state_name}.true_values",
                )
            ),
        )
    return ModalityInventorySpec(
        type=inventory_type,
        path=_path(value["path"], f"{label}.path", base),
        protein_id_column=_text(
            value.get("protein_id_column", "protein_id"),
            f"{label}.protein_id_column",
        ),
        modality_column=_text(
            value.get("modality_column", "modality"),
            f"{label}.modality_column",
        ),
        states=states,
    )


def _category_source(raw: Any, label: str, base: Path) -> CategorySourceSpec:
    value = _mapping(raw, label)
    _keys(
        value,
        label,
        required=("name", "path", "protein_id_column", "category_id_column"),
        optional=("category_name_column",),
    )
    return CategorySourceSpec(
        name=_text(value["name"], f"{label}.name"),
        path=_path(value["path"], f"{label}.path", base),
        protein_id_column=_text(
            value["protein_id_column"], f"{label}.protein_id_column"
        ),
        category_id_column=_text(
            value["category_id_column"], f"{label}.category_id_column"
        ),
        category_name_column=_text(
            value.get("category_name_column", ""),
            f"{label}.category_name_column",
            allow_empty=True,
        ),
    )


def _dataset(raw: Any, index: int, base: Path) -> DatasetConfig:
    label = f"datasets[{index}]"
    value = _mapping(raw, label)
    _keys(
        value,
        label,
        required=("id", "benchmark_dir", "obo_file"),
        optional=(
            "allow_legacy_singular_protein_header",
            "allow_all_zero_rows",
            "split_overlap_policy",
            "source_annotations",
            "taxonomy_sources",
            "modality_inventory",
            "category_sources",
        ),
    )
    overlap = _text(
        value.get("split_overlap_policy", "disallow"),
        f"{label}.split_overlap_policy",
    )
    if overlap not in {"allow", "disallow"}:
        raise ConfigError(f"{label}.split_overlap_policy must be allow or disallow")
    taxonomy_raw = value.get("taxonomy_sources", [])
    category_raw = value.get("category_sources", [])
    if not isinstance(taxonomy_raw, list):
        raise ConfigError(f"{label}.taxonomy_sources must be a JSON array")
    if not isinstance(category_raw, list):
        raise ConfigError(f"{label}.category_sources must be a JSON array")
    category_sources = tuple(
        _category_source(item, f"{label}.category_sources[{number}]", base)
        for number, item in enumerate(category_raw)
    )
    category_names = [source.name for source in category_sources]
    if len(category_names) != len(set(category_names)):
        raise ConfigError(f"{label}.category_sources names must be unique")
    return DatasetConfig(
        id=_text(value["id"], f"{label}.id"),
        benchmark_dir=_path(value["benchmark_dir"], f"{label}.benchmark_dir", base),
        obo_file=_path(value["obo_file"], f"{label}.obo_file", base),
        allow_legacy_singular_protein_header=_boolean(
            value.get("allow_legacy_singular_protein_header", False),
            f"{label}.allow_legacy_singular_protein_header",
        ),
        allow_all_zero_rows=_boolean(
            value.get("allow_all_zero_rows", False),
            f"{label}.allow_all_zero_rows",
        ),
        split_overlap_policy=overlap,
        source_annotations=(
            _source_annotations(
                value["source_annotations"], f"{label}.source_annotations", base
            )
            if "source_annotations" in value
            else None
        ),
        taxonomy_sources=tuple(
            _taxonomy_source(item, f"{label}.taxonomy_sources[{number}]", base)
            for number, item in enumerate(taxonomy_raw)
        ),
        modality_inventory=(
            _modality_inventory(
                value["modality_inventory"], f"{label}.modality_inventory", base
            )
            if "modality_inventory" in value
            else None
        ),
        category_sources=category_sources,
    )


def load_config(path: Path) -> RunConfig:
    source_path = path.expanduser().resolve()
    try:
        raw = json.loads(source_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ConfigError(f"Invalid JSON in {source_path}: {exc}") from exc
    value = _mapping(raw, "configuration root")
    _keys(
        value,
        "configuration root",
        required=("schema_version", "run_name", "datasets"),
        optional=("top_n",),
    )
    if value["schema_version"] != 1:
        raise ConfigError("schema_version must be exactly 1")
    datasets_raw = value["datasets"]
    if not isinstance(datasets_raw, list) or not datasets_raw:
        raise ConfigError("datasets must be a non-empty JSON array")
    datasets = tuple(
        _dataset(item, index, source_path.parent)
        for index, item in enumerate(datasets_raw)
    )
    ids = [dataset.id for dataset in datasets]
    if len(ids) != len(set(ids)):
        raise ConfigError("dataset IDs must be unique")
    top_n = value.get("top_n", 20)
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
        raise ConfigError("top_n must be a positive integer")
    return RunConfig(
        schema_version=1,
        run_name=_text(value["run_name"], "run_name"),
        top_n=top_n,
        datasets=datasets,
        source_path=source_path,
    )
