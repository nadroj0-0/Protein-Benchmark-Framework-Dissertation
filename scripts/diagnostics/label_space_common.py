#!/usr/bin/env python3
"""Shared, benchmark-agnostic helpers for PFP label-space diagnostics."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import resource
import sys
import tempfile
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


ASPECTS = ("BPO", "CCO", "MFO")
ASPECT_TO_PREFIX = {"BPO": "bp", "CCO": "cc", "MFO": "mf"}
ASPECT_TO_NAMESPACE = {
    "BPO": "biological_process",
    "CCO": "cellular_component",
    "MFO": "molecular_function",
}
ASPECT_TO_ROOT = {
    "BPO": "GO:0008150",
    "CCO": "GO:0005575",
    "MFO": "GO:0003674",
}
CSV_SPLITS = ("training", "validation", "test")
PFP_SPLITS = {"training": "train", "validation": "valid", "test": "test"}
SCHEMA_VERSION = 1


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_snapshot(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Required file is missing: {path}")
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def require_unchanged(path: Path, snapshot: Mapping[str, Any], context: str) -> None:
    if (
        not path.is_file()
        or path.stat().st_size != snapshot["bytes"]
        or sha256_file(path) != snapshot["sha256"]
    ):
        raise ValueError(f"{context} changed while it was being read: {path}")


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def sha256_json(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def atomic_write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_name, path)
    except BaseException:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def required_csv_names() -> list[str]:
    return [
        f"{ASPECT_TO_PREFIX[aspect]}-{split}.csv"
        for aspect in ASPECTS
        for split in CSV_SPLITS
    ]


def quantile(values: Sequence[int | float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(float(value) for value in values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def describe(values: Sequence[int | float]) -> dict[str, float | int | None]:
    if not values:
        return {
            "count": 0,
            "mean": None,
            "median": None,
            "p90": None,
            "p99": None,
            "maximum": None,
        }
    return {
        "count": len(values),
        "mean": sum(float(value) for value in values) / len(values),
        "median": quantile(values, 0.5),
        "p90": quantile(values, 0.9),
        "p99": quantile(values, 0.99),
        "maximum": float(max(values)),
    }


@dataclass(frozen=True)
class OboTerm:
    namespace: str
    parents: tuple[str, ...]


@dataclass
class OboGraph:
    terms: dict[str, OboTerm]

    def validate_roots(self) -> None:
        for aspect, root in ASPECT_TO_ROOT.items():
            term = self.terms.get(root)
            if term is None:
                raise ValueError(f"Ontology does not contain required {aspect} root {root}")
            expected = ASPECT_TO_NAMESPACE[aspect]
            if term.namespace != expected:
                raise ValueError(
                    f"Ontology root {root} has namespace {term.namespace!r}; expected {expected!r}"
                )

    def depths(self, aspect: str) -> tuple[dict[str, int], dict[str, int]]:
        namespace = ASPECT_TO_NAMESPACE[aspect]
        root = ASPECT_TO_ROOT[aspect]
        shortest: dict[str, int] = {root: 0}
        longest: dict[str, int] = {root: 0}
        visiting: set[str] = set()

        def visit(term_id: str) -> tuple[int, int] | None:
            if term_id in shortest:
                return shortest[term_id], longest[term_id]
            term = self.terms.get(term_id)
            if term is None or term.namespace != namespace:
                return None
            if term_id in visiting:
                raise ValueError(f"Cycle detected in GO parent graph at {term_id}")
            visiting.add(term_id)
            candidates = []
            for parent in term.parents:
                result = visit(parent)
                if result is not None:
                    candidates.append(result)
            visiting.remove(term_id)
            if not candidates:
                return None
            shortest[term_id] = 1 + min(value[0] for value in candidates)
            longest[term_id] = 1 + max(value[1] for value in candidates)
            return shortest[term_id], longest[term_id]

        for term_id in self.terms:
            visit(term_id)
        return shortest, longest

    def ancestor_closure(self, aspect: str) -> dict[str, frozenset[str]]:
        namespace = ASPECT_TO_NAMESPACE[aspect]
        root = ASPECT_TO_ROOT[aspect]
        memo: dict[str, frozenset[str]] = {root: frozenset({root})}
        visiting: set[str] = set()

        def visit(term_id: str) -> frozenset[str]:
            if term_id in memo:
                return memo[term_id]
            term = self.terms.get(term_id)
            if term is None or term.namespace != namespace:
                return frozenset()
            if term_id in visiting:
                raise ValueError(f"Cycle detected in GO parent graph at {term_id}")
            visiting.add(term_id)
            values = {term_id}
            for parent in term.parents:
                parent_term = self.terms.get(parent)
                if parent_term is not None and parent_term.namespace == namespace:
                    values.update(visit(parent))
            visiting.remove(term_id)
            memo[term_id] = frozenset(values)
            return memo[term_id]

        for term_id, term in self.terms.items():
            if term.namespace == namespace:
                visit(term_id)
        return memo


def read_obo(path: Path) -> OboGraph:
    if not path.is_file():
        raise FileNotFoundError(f"GO OBO file is missing: {path}")
    terms: dict[str, OboTerm] = {}
    stanza: dict[str, Any] = {}

    def publish() -> None:
        term_id = stanza.get("id")
        namespace = stanza.get("namespace")
        if (
            stanza.get("type") == "Term"
            and term_id
            and namespace
            and stanza.get("is_obsolete") != "true"
        ):
            if str(term_id) in terms:
                raise ValueError(f"Ontology contains duplicate live term {term_id}")
            terms[str(term_id)] = OboTerm(
                namespace=str(namespace),
                parents=tuple(sorted(set(stanza.get("parents", [])))),
            )

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                publish()
                stanza = {"type": "Term", "parents": []}
            elif line.startswith("["):
                publish()
                stanza = {}
            elif stanza.get("type") == "Term":
                if line.startswith("id: "):
                    stanza["id"] = line[4:].strip()
                elif line.startswith("namespace: "):
                    stanza["namespace"] = line[11:].strip()
                elif line.startswith("is_obsolete: "):
                    stanza["is_obsolete"] = line[13:].strip()
                elif line.startswith("is_a: "):
                    stanza["parents"].append(line[6:].split()[0])
                elif line.startswith("relationship: part_of "):
                    stanza["parents"].append(line[22:].split()[0])
        publish()
    if not terms:
        raise ValueError(f"No live GO terms were parsed from {path}")
    graph = OboGraph(terms)
    graph.validate_roots()
    return graph


def read_ia_file(path: Path | None) -> dict[str, float]:
    if path is None:
        return {}
    values: dict[str, float] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(f"Malformed IA row {path}:{line_number}")
            value = float(fields[1])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid IA value {path}:{line_number}")
            if fields[0] in values:
                raise ValueError(f"Duplicate IA term {path}:{line_number}: {fields[0]}")
            values[fields[0]] = value
    if not values:
        raise ValueError(f"IA file contains no values: {path}")
    return values


def _update_ordered_digests(
    label_digest: Any,
    sequence_digest: Any,
    record_digest: Any,
    protein_id: str,
    sequence: str,
    positive_indices: Sequence[int],
) -> None:
    sequence_sha = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
    indices = ",".join(str(value) for value in positive_indices)
    label_digest.update(f"{protein_id}\t{indices}\n".encode("utf-8"))
    sequence_digest.update(f"{protein_id}\t{sequence_sha}\n".encode("utf-8"))
    record_digest.update(
        f"{protein_id}\t{sequence_sha}\t{indices}\n".encode("utf-8")
    )


def audit_csv(
    path: Path,
    aspect: str,
    split: str,
    graph: OboGraph,
    ia_values: Mapping[str, float],
    allow_singular_header: bool,
    allow_all_zero_rows: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], Counter[str]]:
    initial_bytes = path.stat().st_size
    initial_sha256 = sha256_file(path)
    namespace = ASPECT_TO_NAMESPACE[aspect]
    root = ASPECT_TO_ROOT[aspect]
    shortest, longest = graph.depths(aspect)
    ancestors = graph.ancestor_closure(aspect)
    term_support: Counter[str] = Counter()
    root_only_rows: list[dict[str, Any]] = []
    term_rows: list[dict[str, Any]] = []
    label_counts: list[int] = []
    non_root_counts: list[int] = []
    row_max_shortest_depths: list[int] = []
    row_max_longest_depths: list[int] = []
    row_max_ia: list[float] = []
    root_recall_values: list[float] = []
    depth_histogram: Counter[str] = Counter()
    label_digest = hashlib.sha256()
    sequence_digest = hashlib.sha256()
    record_digest = hashlib.sha256()
    protein_ids: set[str] = set()
    all_zero_rows = 0
    root_positive_rows = 0
    alias = None

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Required CSV is empty: {path}") from exc
        if len(header) < 3 or header[1] != "sequences":
            raise ValueError(f"{path.name} must begin with proteins,sequences")
        if header[0] == "protein":
            if not allow_singular_header:
                raise ValueError(
                    f"{path.name} uses legacy singular protein header but policy forbids it"
                )
            alias = {"file": path.name, "source": "protein", "normalized": "proteins"}
        elif header[0] != "proteins":
            raise ValueError(f"{path.name} has unsupported protein header {header[0]!r}")
        terms = header[2:]
        if len(terms) != len(set(terms)):
            raise ValueError(f"{path.name} contains duplicate GO columns")
        if root not in terms:
            raise ValueError(f"{path.name} does not contain ontology root {root}")
        term_set = set(terms)
        for term in terms:
            ontology_term = graph.terms.get(term)
            if ontology_term is None:
                raise ValueError(f"{path.name} contains GO term absent from OBO: {term}")
            if ontology_term.namespace != namespace:
                raise ValueError(
                    f"{path.name} contains {term} from {ontology_term.namespace}, expected {namespace}"
                )
            if term not in shortest or root not in ancestors.get(term, frozenset()):
                raise ValueError(f"{path.name} contains term disconnected from {root}: {term}")
            missing_columns = sorted(ancestors[term] - term_set)
            if missing_columns:
                raise ValueError(
                    f"{path.name} omits required ancestor GO columns for {term}: "
                    f"{missing_columns[:5]}"
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
            labels = row[2:]
            invalid = [value for value in labels if value not in {"0", "1"}]
            if invalid:
                raise ValueError(f"{path.name}:{line_number} contains a non-binary label")
            positive_indices = [index for index, value in enumerate(labels) if value == "1"]
            _update_ordered_digests(
                label_digest,
                sequence_digest,
                record_digest,
                protein_id,
                sequence,
                positive_indices,
            )
            positive_terms = [terms[index] for index in positive_indices]
            positive_set = set(positive_terms)
            missing_ancestors = sorted(
                {
                    ancestor
                    for term in positive_terms
                    for ancestor in ancestors[term]
                    if ancestor not in positive_set
                }
            )
            if missing_ancestors:
                raise ValueError(
                    f"{path.name}:{line_number} labels are not ancestor-closed; "
                    f"missing positives: {missing_ancestors[:5]}"
                )
            for term in positive_terms:
                term_support[term] += 1
            label_count = len(positive_terms)
            label_counts.append(label_count)
            if label_count == 0:
                all_zero_rows += 1
                non_root_counts.append(0)
                root_recall_values.append(0.0)
                depth_histogram["all-zero"] += 1
                continue
            root_present = root in positive_terms
            if root_present:
                root_positive_rows += 1
            root_recall_values.append((1.0 / label_count) if root_present else 0.0)
            non_root = [term for term in positive_terms if term != root]
            non_root_counts.append(len(non_root))
            if not non_root:
                depth_histogram["root-only"] += 1
                root_only_rows.append(
                    {
                        "aspect": aspect,
                        "split": split,
                        "protein_id": protein_id,
                        "sequence_sha256": hashlib.sha256(
                            sequence.encode("utf-8")
                        ).hexdigest(),
                    }
                )
                continue
            shortest_depths = [shortest[term] for term in non_root if term in shortest]
            longest_depths = [longest[term] for term in non_root if term in longest]
            if len(shortest_depths) != len(non_root):
                missing = sorted(set(non_root) - set(shortest))
                raise ValueError(
                    f"{path.name}:{line_number} contains terms disconnected from {root}: {missing[:5]}"
                )
            row_max_shortest = max(shortest_depths)
            row_max_longest = max(longest_depths)
            row_max_shortest_depths.append(row_max_shortest)
            row_max_longest_depths.append(row_max_longest)
            depth_histogram[str(row_max_shortest)] += 1
            if ia_values:
                present_ia = [ia_values[term] for term in non_root if term in ia_values]
                if present_ia:
                    row_max_ia.append(max(present_ia))

    row_count = len(label_counts)
    if row_count == 0:
        raise ValueError(f"Benchmark split contains no rows: {path.name}")
    if all_zero_rows and not allow_all_zero_rows:
        raise ValueError(
            f"{path.name} contains {all_zero_rows} all-zero rows but policy forbids them"
        )
    macro_precision = root_positive_rows / row_count
    macro_recall = sum(root_recall_values) / row_count
    baseline_f = (
        2.0 * macro_precision * macro_recall / (macro_precision + macro_recall)
        if macro_precision + macro_recall
        else 0.0
    )
    ranked_support = sorted(term_support.values(), reverse=True)
    positive_total = sum(ranked_support)
    for term in terms:
        term_rows.append(
            {
                "aspect": aspect,
                "split": split,
                "term": term,
                "is_root": term == root,
                "support": term_support[term],
                "shortest_depth": shortest.get(term),
                "longest_depth": longest.get(term),
                "ia": ia_values.get(term),
            }
        )
    final_bytes = path.stat().st_size
    final_sha256 = sha256_file(path)
    if final_bytes != initial_bytes or final_sha256 != initial_sha256:
        raise ValueError(f"Benchmark CSV changed while it was being audited: {path}")
    report = {
        "file": path.name,
        "bytes": final_bytes,
        "sha256": final_sha256,
        "source_header": header[:2],
        "header_alias": alias,
        "rows": row_count,
        "terms": len(terms),
        "ordered_terms": terms,
        "ordered_terms_sha256": sha256_json(terms),
        "ordered_labels_sha256": label_digest.hexdigest(),
        "ordered_sequences_sha256": sequence_digest.hexdigest(),
        "ordered_records_sha256": record_digest.hexdigest(),
        "positive_labels": positive_total,
        "all_zero_rows": all_zero_rows,
        "root": root,
        "root_positive_rows": root_positive_rows,
        "root_only_rows": len(root_only_rows),
        "root_only_fraction": len(root_only_rows) / row_count,
        "rows_with_non_root_labels": row_count - len(root_only_rows) - all_zero_rows,
        "labels_per_protein": describe(label_counts),
        "non_root_labels_per_protein": describe(non_root_counts),
        "max_shortest_depth_per_non_root_row": describe(row_max_shortest_depths),
        "max_longest_depth_per_non_root_row": describe(row_max_longest_depths),
        "max_ia_per_non_root_row": describe(row_max_ia) if ia_values else None,
        "ia_terms_missing": (
            sorted(term for term in terms if term not in ia_values)
            if ia_values
            else None
        ),
        "terms_with_support": sum(value > 0 for value in term_support.values()),
        "positive_mass_top10_terms": (
            sum(ranked_support[:10]) / positive_total if positive_total else 0.0
        ),
        "positive_mass_top100_terms": (
            sum(ranked_support[:100]) / positive_total if positive_total else 0.0
        ),
        "root_only_diagnostic_baseline": {
            "macro_precision": macro_precision,
            "macro_recall": macro_recall,
            "macro_f": baseline_f,
            "note": "Arithmetic diagnostic on the retained label matrix; not a cafaeval result.",
        },
    }
    return report, root_only_rows, term_rows, depth_histogram


def verify_prepared_data(
    prepared_dir: Path,
    label: str,
    csv_reports: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    try:
        import numpy as np
        from scipy.sparse import load_npz
    except ImportError as exc:
        raise RuntimeError("Prepared-data verification requires numpy and scipy") from exc

    started = time.perf_counter()
    result: dict[str, Any] = {"label": label, "directory": str(prepared_dir), "files": {}}
    for aspect in ASPECTS:
        terms_path = prepared_dir / f"{aspect}_go_terms.json"
        terms_snapshot = file_snapshot(terms_path)
        terms = json.loads(terms_path.read_text(encoding="utf-8"))
        if (
            not isinstance(terms, list)
            or not all(isinstance(term, str) for term in terms)
            or len(terms) != len(set(terms))
        ):
            raise ValueError(f"Prepared GO terms are invalid for {aspect} in {label}")
        for csv_split, pfp_split in PFP_SPLITS.items():
            key = f"{aspect}:{csv_split}"
            source = csv_reports[key]
            names_path = prepared_dir / f"{aspect}_{pfp_split}_names.npy"
            labels_path = prepared_dir / f"{aspect}_{pfp_split}_labels.npz"
            sequences_path = prepared_dir / f"{aspect}_{pfp_split}_sequences.json"
            snapshots = {
                path: file_snapshot(path)
                for path in (names_path, labels_path, sequences_path)
            }
            if terms != source["ordered_terms"]:
                raise ValueError(f"Prepared GO terms differ from CSV for {key} in {label}")
            names_array = np.load(names_path, allow_pickle=True)
            if names_array.ndim != 1 or names_array.dtype.kind not in {"O", "U"}:
                raise ValueError(f"Prepared names have unsupported shape or dtype for {key}")
            if not all(isinstance(name, (str, np.str_)) for name in names_array.tolist()):
                raise ValueError(f"Prepared names contain non-string values for {key}")
            names = [str(name) for name in names_array.tolist()]
            sequences = json.loads(sequences_path.read_text(encoding="utf-8"))
            if not isinstance(sequences, dict):
                raise ValueError(f"Prepared sequences are not an object for {key} in {label}")
            matrix = load_npz(labels_path).tocsr()
            matrix.sum_duplicates()
            matrix.eliminate_zeros()
            matrix.sort_indices()
            if matrix.shape != (len(names), len(terms)):
                raise ValueError(f"Prepared label shape mismatch for {key} in {label}")
            if not np.isfinite(matrix.data).all() or not np.equal(matrix.data, 1).all():
                raise ValueError(f"Prepared labels are not binary for {key} in {label}")
            if len(set(names)) != len(names):
                raise ValueError(f"Prepared names are duplicated for {key} in {label}")
            if set(sequences) != set(names):
                raise ValueError(
                    f"Prepared sequence membership differs from names for {key} in {label}"
                )
            label_digest = hashlib.sha256()
            sequence_digest = hashlib.sha256()
            record_digest = hashlib.sha256()
            for index, protein_id in enumerate(names):
                sequence = sequences.get(protein_id)
                if not isinstance(sequence, str):
                    raise ValueError(f"Prepared sequence is missing for {protein_id} in {label}")
                start = matrix.indptr[index]
                stop = matrix.indptr[index + 1]
                positive_indices = matrix.indices[start:stop].tolist()
                _update_ordered_digests(
                    label_digest,
                    sequence_digest,
                    record_digest,
                    protein_id,
                    sequence,
                    positive_indices,
                )
            observed = {
                "rows": len(names),
                "terms": len(terms),
                "ordered_labels_sha256": label_digest.hexdigest(),
                "ordered_sequences_sha256": sequence_digest.hexdigest(),
                "ordered_records_sha256": record_digest.hexdigest(),
                "names_sha256": snapshots[names_path]["sha256"],
                "labels_sha256": snapshots[labels_path]["sha256"],
                "sequences_sha256": snapshots[sequences_path]["sha256"],
                "terms_sha256": terms_snapshot["sha256"],
            }
            for path, snapshot in snapshots.items():
                require_unchanged(path, snapshot, "Prepared-data artifact")
            for field in (
                "rows",
                "terms",
                "ordered_labels_sha256",
                "ordered_sequences_sha256",
                "ordered_records_sha256",
            ):
                if observed[field] != source[field]:
                    raise ValueError(
                        f"Prepared-data {field} differs from CSV for {key} in {label}"
                    )
            result["files"][key] = observed
        require_unchanged(terms_path, terms_snapshot, "Prepared GO terms")
    result["passed"] = True
    result["fingerprint"] = sha256_json(result["files"])
    result["wall_seconds"] = time.perf_counter() - started
    return result


def output_manifest(root: Path, exclude: Iterable[str] = ()) -> dict[str, Any]:
    excluded = set(exclude)
    files = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        if relative in excluded:
            continue
        files.append(
            {"path": relative, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        )
    return {
        "schema_version": 1,
        "payload_file_count": len(files),
        "payload_bytes": sum(item["bytes"] for item in files),
        "files": files,
    }
