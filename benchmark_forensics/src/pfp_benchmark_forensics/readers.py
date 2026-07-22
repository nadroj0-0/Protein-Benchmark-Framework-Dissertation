"""Strict readers for benchmark CSVs and optional evidence sources."""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, Iterator, Mapping, MutableMapping, Set, TextIO, Tuple

from .config import (
    CategorySourceSpec,
    DatasetConfig,
    ModalityInventorySpec,
    SourceAnnotationSpec,
    TaxonomySourceSpec,
)
from .models import Observation, SourceRecord, Taxon, TaxonomyConflict


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open("r", encoding="utf-8", newline="") as handle:
            yield handle


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def file_snapshot(path: Path) -> dict:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required {label} does not exist: {path}")


def read_benchmark_csv(
    path: Path,
    *,
    dataset_id: str,
    aspect: str,
    split: str,
    root: str,
    expected_namespace: str,
    ontology: object,
    allow_singular_header: bool,
    allow_all_zero_rows: bool,
) -> Tuple[Tuple[Observation, ...], dict, Tuple[str, ...], Mapping[str, str]]:
    _require_file(path, "benchmark CSV")
    initial = file_snapshot(path)
    observations = []
    sequences = {}
    protein_ids: Set[str] = set()
    term_support: MutableMapping[str, int] = defaultdict(int)
    label_counts = []
    non_root_counts = []
    sequence_lengths = []
    source_header = None
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Benchmark CSV is empty: {path}") from exc
        if len(header) < 3 or header[1] != "sequences":
            raise ValueError(f"{path.name} must begin with proteins,sequences")
        source_header = header[0]
        if header[0] == "protein":
            if not allow_singular_header:
                raise ValueError(f"{path.name} uses forbidden legacy protein header")
        elif header[0] != "proteins":
            raise ValueError(f"{path.name} has unsupported first column {header[0]!r}")
        terms = tuple(header[2:])
        if len(terms) != len(set(terms)):
            raise ValueError(f"{path.name} contains duplicate GO columns")
        if root not in terms:
            raise ValueError(f"{path.name} omits ontology root {root}")
        term_set = set(terms)
        ancestors = {}
        for term in terms:
            canonical = ontology.resolve(term)
            if canonical is None:
                raise ValueError(
                    f"{path.name} contains GO term absent from OBO: {term}"
                )
            if ontology.namespace(term) != expected_namespace:
                raise ValueError(
                    f"{path.name} contains {term} outside {expected_namespace}"
                )
            ancestors[term] = ontology.ancestors(term)
            missing = sorted(ancestors[term] - term_set)
            if missing:
                raise ValueError(
                    f"{path.name} omits ancestor columns required by {term}: {missing[:5]}"
                )
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"{path.name}:{line_number} has {len(row)} columns; expected {len(header)}"
                )
            protein_id, sequence = row[:2]
            if not protein_id or protein_id in protein_ids:
                raise ValueError(
                    f"{path.name}:{line_number} has empty or duplicate protein ID {protein_id!r}"
                )
            if not sequence:
                raise ValueError(f"{path.name}:{line_number} has an empty sequence")
            protein_ids.add(protein_id)
            sequences[protein_id] = sequence
            labels = row[2:]
            if any(value not in {"0", "1"} for value in labels):
                raise ValueError(
                    f"{path.name}:{line_number} contains non-binary labels"
                )
            positive = tuple(term for term, value in zip(terms, labels) if value == "1")
            positive_set = set(positive)
            missing_positive_ancestors = sorted(
                ancestor
                for term in positive
                for ancestor in ancestors[term]
                if ancestor not in positive_set
            )
            if missing_positive_ancestors:
                raise ValueError(
                    f"{path.name}:{line_number} labels are not ancestor-closed: "
                    f"{missing_positive_ancestors[:5]}"
                )
            for term in positive:
                term_support[term] += 1
            non_root = tuple(term for term in positive if term != root)
            all_zero = not positive
            root_positive = root in positive_set
            root_only = root_positive and not non_root
            if all_zero and not allow_all_zero_rows:
                raise ValueError(
                    f"{path.name}:{line_number} is all-zero but policy forbids it"
                )
            sequence_sha = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
            observations.append(
                Observation(
                    dataset_id=dataset_id,
                    aspect=aspect,
                    split=split,
                    protein_id=protein_id,
                    sequence_sha256=sequence_sha,
                    sequence_length=len(sequence),
                    label_count=len(positive),
                    non_root_label_count=len(non_root),
                    root_positive=root_positive,
                    root_only=root_only,
                    all_zero=all_zero,
                    positive_terms=positive,
                )
            )
            label_counts.append(len(positive))
            non_root_counts.append(len(non_root))
            sequence_lengths.append(len(sequence))
    if not observations:
        raise ValueError(f"Benchmark split has no rows: {path}")
    final = file_snapshot(path)
    if initial != final:
        raise ValueError(f"Benchmark CSV changed while being read: {path}")
    profile = {
        "dataset_id": dataset_id,
        "aspect": aspect,
        "split": split,
        "file": path.name,
        "path": str(path.resolve()),
        "bytes": final["bytes"],
        "sha256": final["sha256"],
        "source_protein_header": source_header,
        "proteins": len(observations),
        "terms": len(terms),
        "positive_labels": sum(label_counts),
        "root_positive_proteins": sum(item.root_positive for item in observations),
        "root_only_proteins": sum(item.root_only for item in observations),
        "all_zero_proteins": sum(item.all_zero for item in observations),
        "terms_with_support": sum(term_support[term] > 0 for term in terms),
    }
    return tuple(observations), profile, terms, sequences


def _normalise_annotations(value: object, separator: str) -> Tuple[str, ...]:
    if value is None:
        return tuple()
    if isinstance(value, str):
        return tuple(
            sorted(term.strip() for term in value.split(separator) if term.strip())
        )
    if isinstance(value, (list, tuple, set, frozenset)):
        if not all(isinstance(term, str) for term in value):
            raise ValueError("Source annotation collection contains a non-string term")
        return tuple(sorted(set(value)))
    try:
        if value != value:  # pandas NaN
            return tuple()
    except Exception:
        pass
    raise ValueError(f"Unsupported source annotation value: {type(value).__name__}")


def load_source_annotations(
    spec: SourceAnnotationSpec,
) -> Tuple[Mapping[str, Mapping[str, SourceRecord]], Tuple[Path, ...]]:
    import pandas as pd

    by_split: Dict[str, Dict[str, SourceRecord]] = {
        "training": {},
        "validation": {},
        "test": {},
    }
    input_paths = []
    if spec.type == "pfp-pickle-directory":
        for split, filename in spec.split_files.items():
            path = spec.path / filename
            _require_file(path, f"source annotation pickle for {split}")
            input_paths.append(path)
            frame = pd.read_pickle(path)
            required = {
                spec.protein_id_column,
                spec.sequence_column,
                spec.annotations_column,
            }
            missing = sorted(required - set(frame.columns))
            if missing:
                raise ValueError(f"{path} is missing source columns: {missing}")
            for row in frame[
                [spec.protein_id_column, spec.sequence_column, spec.annotations_column]
            ].itertuples(index=False, name=None):
                protein_id, sequence, annotations = row
                if not isinstance(protein_id, str) or not protein_id:
                    raise ValueError(f"{path} contains an invalid protein ID")
                if protein_id in by_split[split]:
                    raise ValueError(f"{path} duplicates source protein {protein_id}")
                if not isinstance(sequence, str) or not sequence:
                    raise ValueError(
                        f"{path} contains an invalid sequence for {protein_id}"
                    )
                by_split[split][protein_id] = SourceRecord(
                    protein_id,
                    sequence,
                    _normalise_annotations(annotations, spec.annotation_separator),
                )
    else:
        path = spec.path
        _require_file(path, "source annotation TSV")
        input_paths.append(path)
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                spec.protein_id_column,
                spec.split_column,
                spec.sequence_column,
                spec.annotations_column,
            }
            missing = sorted(required - set(reader.fieldnames or ()))
            if missing:
                raise ValueError(f"{path} is missing source columns: {missing}")
            for line_number, row in enumerate(reader, start=2):
                split = row[spec.split_column]
                if split not in by_split:
                    raise ValueError(
                        f"{path}:{line_number} has unsupported split {split!r}"
                    )
                protein_id = row[spec.protein_id_column]
                if not protein_id or protein_id in by_split[split]:
                    raise ValueError(
                        f"{path}:{line_number} has empty or duplicate protein ID"
                    )
                sequence = row[spec.sequence_column]
                if not sequence:
                    raise ValueError(f"{path}:{line_number} has an empty sequence")
                by_split[split][protein_id] = SourceRecord(
                    protein_id,
                    sequence,
                    _normalise_annotations(
                        row[spec.annotations_column], spec.annotation_separator
                    ),
                )
    return by_split, tuple(input_paths)


def _record_taxon(
    result: MutableMapping[str, MutableMapping[str, Taxon]],
    protein_id: str,
    taxon: Taxon,
) -> None:
    if not protein_id:
        return
    by_source = result.setdefault(protein_id, {})
    previous = by_source.get(taxon.source_name)
    if previous is not None and previous.taxon_id != taxon.taxon_id:
        raise ValueError(
            f"Taxonomy source {taxon.source_name!r} assigns conflicting taxa to "
            f"{protein_id}: {previous.taxon_id} vs {taxon.taxon_id}"
        )
    if previous is None or (not previous.taxon_name and taxon.taxon_name):
        by_source[taxon.source_name] = taxon


def _load_taxonomy_tsv(
    spec: TaxonomySourceSpec,
    wanted_ids: Set[str],
    result: MutableMapping[str, MutableMapping[str, Taxon]],
) -> None:
    with open_text(spec.path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = set(spec.id_columns) | {spec.taxon_id_column}
        if spec.taxon_name_column:
            required.add(spec.taxon_name_column)
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{spec.path} is missing taxonomy columns: {missing}")
        for line_number, row in enumerate(reader, start=2):
            taxon_id = row[spec.taxon_id_column].strip()
            if not taxon_id:
                continue
            taxon_name = (
                row[spec.taxon_name_column].strip() if spec.taxon_name_column else ""
            )
            taxon = Taxon(
                taxon_id,
                taxon_name,
                str(spec.path.resolve()),
                spec.name,
                spec.priority,
            )
            for column in spec.id_columns:
                for protein_id in re.split(r"[;,]", row[column]):
                    protein_id = protein_id.strip()
                    if protein_id in wanted_ids:
                        _record_taxon(result, protein_id, taxon)


def _load_uniprot_dat(
    spec: TaxonomySourceSpec,
    wanted_ids: Set[str],
    result: MutableMapping[str, MutableMapping[str, Taxon]],
) -> None:
    accessions = []
    taxon_id = ""
    organism_parts = []

    def publish() -> None:
        if not taxon_id:
            return
        organism = " ".join(organism_parts).strip()
        taxon = Taxon(
            taxon_id,
            organism,
            str(spec.path.resolve()),
            spec.name,
            spec.priority,
        )
        for accession in accessions:
            if accession in wanted_ids:
                _record_taxon(result, accession, taxon)

    with open_text(spec.path) as handle:
        for raw in handle:
            if raw.startswith("AC   "):
                accessions.extend(
                    item.strip() for item in raw[5:].split(";") if item.strip()
                )
            elif raw.startswith("OX   "):
                match = re.search(r"NCBI_TaxID=(\d+)", raw)
                if match:
                    taxon_id = match.group(1)
            elif raw.startswith("OS   "):
                organism_parts.append(raw[5:].strip())
            elif raw.startswith("//"):
                publish()
                accessions = []
                taxon_id = ""
                organism_parts = []
        publish()


def _resolve_taxonomy(
    candidates: Mapping[str, Mapping[str, Taxon]],
) -> Tuple[Mapping[str, Taxon], Tuple[TaxonomyConflict, ...]]:
    resolved: Dict[str, Taxon] = {}
    conflicts = []
    for protein_id in sorted(candidates):
        ranked = sorted(
            candidates[protein_id].values(),
            key=lambda item: (
                -item.source_priority,
                not bool(item.taxon_name),
                item.source_name,
                item.source,
            ),
        )
        highest_priority = ranked[0].source_priority
        highest = [item for item in ranked if item.source_priority == highest_priority]
        highest_ids = sorted({item.taxon_id for item in highest})
        if len(highest_ids) > 1:
            details = ", ".join(
                f"{item.source_name}={item.taxon_id}" for item in highest
            )
            raise ValueError(
                f"Equal-priority taxonomy sources disagree for {protein_id} "
                f"at priority {highest_priority}: {details}"
            )
        selected = highest[0]
        resolved[protein_id] = selected
        for alternative in ranked:
            if alternative.source_name == selected.source_name:
                continue
            if alternative.taxon_id == selected.taxon_id:
                continue
            conflicts.append(
                TaxonomyConflict(
                    protein_id=protein_id,
                    selected=selected,
                    alternative=alternative,
                    resolution="selected explicitly higher-priority taxonomy source",
                )
            )
    return resolved, tuple(conflicts)


def load_taxonomy(
    specs: Iterable[TaxonomySourceSpec], wanted_ids: Set[str]
) -> Tuple[Mapping[str, Taxon], Tuple[Path, ...], Tuple[TaxonomyConflict, ...]]:
    candidates: Dict[str, MutableMapping[str, Taxon]] = {}
    paths = []
    for spec in specs:
        _require_file(spec.path, "taxonomy source")
        paths.append(spec.path)
        if spec.type == "tsv":
            _load_taxonomy_tsv(spec, wanted_ids, candidates)
        else:
            _load_uniprot_dat(spec, wanted_ids, candidates)
    result, conflicts = _resolve_taxonomy(candidates)
    return result, tuple(paths), conflicts


def load_modality_inventory(
    spec: ModalityInventorySpec,
) -> Tuple[Mapping[Tuple[str, str], Mapping[str, bool]], Tuple[str, ...]]:
    _require_file(spec.path, "modality inventory")
    records: Dict[Tuple[str, str], Mapping[str, bool]] = {}
    modalities = set()
    with open_text(spec.path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {spec.protein_id_column, spec.modality_column}
        required.update(state.column for state in spec.states.values())
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{spec.path} is missing modality columns: {missing}")
        for line_number, row in enumerate(reader, start=2):
            protein_id = row[spec.protein_id_column].strip()
            modality = row[spec.modality_column].strip()
            if not protein_id or not modality:
                raise ValueError(f"{spec.path}:{line_number} has an empty pair key")
            key = (protein_id, modality)
            if key in records:
                raise ValueError(f"{spec.path}:{line_number} duplicates {key}")
            records[key] = {
                name: row[state.column].strip().casefold() in state.true_values
                for name, state in spec.states.items()
            }
            modalities.add(modality)
    if not records:
        raise ValueError(f"Modality inventory contains no records: {spec.path}")
    return records, tuple(sorted(modalities))


def load_category_source(
    spec: CategorySourceSpec,
) -> Mapping[str, Tuple[Tuple[str, str], ...]]:
    _require_file(spec.path, f"{spec.name} category source")
    result: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    with open_text(spec.path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {spec.protein_id_column, spec.category_id_column}
        if spec.category_name_column:
            required.add(spec.category_name_column)
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{spec.path} is missing category columns: {missing}")
        for line_number, row in enumerate(reader, start=2):
            protein_id = row[spec.protein_id_column].strip()
            category_id = row[spec.category_id_column].strip()
            category_name = (
                row[spec.category_name_column].strip()
                if spec.category_name_column
                else ""
            )
            if not protein_id or not category_id:
                raise ValueError(f"{spec.path}:{line_number} has an empty category key")
            result[protein_id].add((category_id, category_name))
    return {protein: tuple(sorted(values)) for protein, values in result.items()}


def input_paths_for_dataset(config: DatasetConfig) -> Tuple[Path, ...]:
    paths = [config.obo_file]
    from .ontology import ASPECTS, ASPECT_TO_FILE, SPLITS

    for aspect in ASPECTS:
        for split in SPLITS:
            paths.append(config.benchmark_dir / f"{ASPECT_TO_FILE[aspect]}-{split}.csv")
    if config.source_annotations:
        if config.source_annotations.type == "pfp-pickle-directory":
            paths.extend(
                config.source_annotations.path / name
                for name in config.source_annotations.split_files.values()
            )
        else:
            paths.append(config.source_annotations.path)
    paths.extend(spec.path for spec in config.taxonomy_sources)
    if config.modality_inventory:
        paths.append(config.modality_inventory.path)
    paths.extend(spec.path for spec in config.category_sources)
    return tuple(paths)
