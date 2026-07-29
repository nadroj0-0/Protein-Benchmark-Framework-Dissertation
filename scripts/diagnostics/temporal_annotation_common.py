#!/usr/bin/env python3
"""Policy-neutral temporal annotation ledgers and cohort masks."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from label_space_common import ASPECT_TO_ROOT, file_snapshot, require_unchanged


ASPECTS = ("BPO", "CCO", "MFO")
GO_ID = re.compile(r"GO:\d{7}")

AnnotationMap = dict[str, dict[str, frozenset[str]]]
EXPOSURE_FIELDS = (
    "protein_id",
    "train_id_member",
    "valid_id_member",
    "train_sequence_member",
    "valid_sequence_member",
    "train_homology_cluster_member",
    "modality_availability",
    "feature_temporal_policy",
)


def read_protein_scope(path: Path) -> tuple[list[str], dict[str, Any]]:
    snapshot = file_snapshot(path)
    protein_ids: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if reader.fieldnames != ["protein_id"]:
            raise ValueError(f"{path} must have the single header protein_id")
        for line_number, row in enumerate(reader, start=2):
            protein_id = row["protein_id"].strip()
            if not protein_id:
                raise ValueError(f"Empty protein ID at {path}:{line_number}")
            if protein_id in seen:
                raise ValueError(
                    f"Duplicate protein ID at {path}:{line_number}: {protein_id}"
                )
            seen.add(protein_id)
            protein_ids.append(protein_id)
    if not protein_ids:
        raise ValueError(f"Protein scope is empty: {path}")
    require_unchanged(path, snapshot, "Protein scope")
    return protein_ids, {**snapshot, "protein_count": len(protein_ids)}


def read_annotation_rows(
    path: Path, protein_scope: set[str] | None = None
) -> tuple[AnnotationMap, dict[str, Any]]:
    snapshot = file_snapshot(path)
    values: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {aspect: set() for aspect in ASPECTS}
    )
    seen: set[tuple[str, str, str]] = set()
    term_aspects: dict[str, str] = {}
    rows = 0
    selected_rows = 0
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if reader.fieldnames != ["protein_id", "aspect", "go_term"]:
            raise ValueError(
                f"{path} must have headers protein_id, aspect, go_term in that order"
            )
        for line_number, row in enumerate(reader, start=2):
            rows += 1
            protein_id = row["protein_id"].strip()
            aspect = row["aspect"].strip()
            term = row["go_term"].strip()
            if not protein_id:
                raise ValueError(f"Empty protein ID at {path}:{line_number}")
            if aspect not in ASPECTS:
                raise ValueError(f"Invalid aspect at {path}:{line_number}: {aspect!r}")
            if GO_ID.fullmatch(term) is None:
                raise ValueError(f"Invalid GO ID at {path}:{line_number}: {term!r}")
            previous_aspect = term_aspects.setdefault(term, aspect)
            if previous_aspect != aspect:
                raise ValueError(
                    f"GO term is assigned to multiple aspects in {path}: {term}"
                )
            key = (protein_id, aspect, term)
            if key in seen:
                raise ValueError(
                    f"Duplicate annotation row at {path}:{line_number}: {key}"
                )
            seen.add(key)
            if protein_scope is None or protein_id in protein_scope:
                values[protein_id][aspect].add(term)
                selected_rows += 1
    require_unchanged(path, snapshot, "Temporal annotation input")
    frozen: AnnotationMap = {
        protein_id: {aspect: frozenset(aspects[aspect]) for aspect in ASPECTS}
        for protein_id, aspects in values.items()
    }
    return frozen, {
        **snapshot,
        "rows": rows,
        "selected_rows": selected_rows,
        "selected_proteins": len(frozen),
        "terms": len(term_aspects),
        "term_aspects": term_aspects,
    }


def read_exposure_rows(
    path: Path, protein_scope: set[str]
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    snapshot = file_snapshot(path)
    rows: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if tuple(reader.fieldnames or ()) != EXPOSURE_FIELDS:
            raise ValueError(
                f"{path} must have exposure headers in the documented order"
            )
        for line_number, row in enumerate(reader, start=2):
            protein_id = row["protein_id"].strip()
            if not protein_id or protein_id in rows:
                raise ValueError(
                    f"Empty or duplicate exposure protein at {path}:{line_number}"
                )
            if protein_id not in protein_scope:
                raise ValueError(
                    f"Exposure protein is outside analysis scope at "
                    f"{path}:{line_number}: {protein_id}"
                )
            value: dict[str, Any] = {}
            for field in (
                "train_id_member",
                "valid_id_member",
                "train_sequence_member",
                "valid_sequence_member",
                "train_homology_cluster_member",
            ):
                observed = row[field].strip().lower()
                if observed not in {"0", "1", "unknown"}:
                    raise ValueError(
                        f"Invalid exposure value at {path}:{line_number}: "
                        f"{field}={observed!r}"
                    )
                value[field] = None if observed == "unknown" else observed == "1"
            for field in ("modality_availability", "feature_temporal_policy"):
                observed = row[field].strip()
                if not observed:
                    raise ValueError(
                        f"Empty exposure provenance at {path}:{line_number}: {field}"
                    )
                value[field] = observed
            rows[protein_id] = value
    require_unchanged(path, snapshot, "Exposure table")
    return rows, {**snapshot, "protein_count": len(rows)}


def validate_term_aspects(
    t0_contract: Mapping[str, Any], t1_contract: Mapping[str, Any]
) -> None:
    t0 = dict(t0_contract["term_aspects"])
    t1 = dict(t1_contract["term_aspects"])
    conflicts = sorted(term for term in set(t0) & set(t1) if t0[term] != t1[term])
    if conflicts:
        raise ValueError(
            "GO terms change aspect between temporal inputs: "
            + ", ".join(conflicts[:5])
        )


def build_ledger_rows(
    protein_ids: Sequence[str],
    t0_annotations: Mapping[str, Mapping[str, frozenset[str]]],
    t1_annotations: Mapping[str, Mapping[str, frozenset[str]]],
) -> list[dict[str, Any]]:
    if len(set(protein_ids)) != len(protein_ids):
        raise ValueError("Temporal ledger protein IDs must be unique")
    rows: list[dict[str, Any]] = []
    for protein_id in protein_ids:
        t0_by_aspect = t0_annotations.get(protein_id, {})
        t1_by_aspect = t1_annotations.get(protein_id, {})
        global_t0 = set().union(
            *(set(t0_by_aspect.get(aspect, frozenset())) for aspect in ASPECTS)
        )
        global_t1 = set().union(
            *(set(t1_by_aspect.get(aspect, frozenset())) for aspect in ASPECTS)
        )
        for aspect in ASPECTS:
            t0_terms = set(t0_by_aspect.get(aspect, frozenset()))
            t1_terms = set(t1_by_aspect.get(aspect, frozenset()))
            gained = t1_terms - t0_terms
            lost = t0_terms - t1_terms
            rows.append(
                {
                    "protein_id": protein_id,
                    "aspect": aspect,
                    "t0_terms": tuple(sorted(t0_terms)),
                    "t1_terms": tuple(sorted(t1_terms)),
                    "gained_terms": tuple(sorted(gained)),
                    "lost_terms": tuple(sorted(lost)),
                    "t0_term_count": len(t0_terms),
                    "t1_term_count": len(t1_terms),
                    "gained_term_count": len(gained),
                    "lost_term_count": len(lost),
                    "global_t0_term_count": len(global_t0),
                    "global_t1_term_count": len(global_t1),
                    "global_t0_empty": not global_t0,
                    "aspect_t0_empty": not t0_terms,
                    "aspect_has_gain": bool(gained),
                    "aspect_has_loss": bool(lost),
                }
            )
    return rows


def build_temporal_state_rows(
    protein_ids: Sequence[str],
    t0_direct: Mapping[str, Mapping[str, frozenset[str]]],
    t1_direct: Mapping[str, Mapping[str, frozenset[str]]],
    t0_closure: Mapping[str, Mapping[str, frozenset[str]]],
    t1_closure: Mapping[str, Mapping[str, frozenset[str]]],
    t0_present: set[str],
    t1_present: set[str],
) -> list[dict[str, Any]]:
    """Build policy-bound direct knowledge states and closure transitions."""
    if len(set(protein_ids)) != len(protein_ids):
        raise ValueError("Temporal state protein IDs must be unique")
    scope = set(protein_ids)
    for label, values, presence in (
        ("t0 direct", t0_direct, t0_present),
        ("t1 direct", t1_direct, t1_present),
        ("t0 closure", t0_closure, t0_present),
        ("t1 closure", t1_closure, t1_present),
    ):
        outside = sorted(set(values) - scope)
        if outside:
            raise ValueError(f"{label} contains proteins outside scope: {outside[:5]}")
        absent = sorted(set(values) - presence)
        if absent:
            raise ValueError(
                f"{label} contains annotations for proteins absent from its "
                f"snapshot: {absent[:5]}"
            )

    rows: list[dict[str, Any]] = []
    for protein_id in protein_ids:
        direct0_by_aspect = t0_direct.get(protein_id, {})
        direct1_by_aspect = t1_direct.get(protein_id, {})
        closure0_by_aspect = t0_closure.get(protein_id, {})
        closure1_by_aspect = t1_closure.get(protein_id, {})
        direct0_nonroot_by_aspect: dict[str, set[str]] = {}
        direct1_nonroot_by_aspect: dict[str, set[str]] = {}
        closure0_nonroot_by_aspect: dict[str, set[str]] = {}
        closure1_nonroot_by_aspect: dict[str, set[str]] = {}
        direct0_roots: set[str] = set()

        for aspect in ASPECTS:
            root = ASPECT_TO_ROOT[aspect]
            direct0 = set(direct0_by_aspect.get(aspect, frozenset()))
            direct1 = set(direct1_by_aspect.get(aspect, frozenset()))
            closure0 = set(closure0_by_aspect.get(aspect, frozenset()))
            closure1 = set(closure1_by_aspect.get(aspect, frozenset()))
            direct0_nonroot = direct0 - {root}
            direct1_nonroot = direct1 - {root}
            closure0_nonroot = closure0 - {root}
            closure1_nonroot = closure1 - {root}
            if not direct0_nonroot.issubset(closure0_nonroot):
                missing = sorted(direct0_nonroot - closure0_nonroot)
                raise ValueError(
                    f"t0 direct terms are absent from t0 closure for "
                    f"{protein_id}/{aspect}: {missing[:5]}"
                )
            if not direct1_nonroot.issubset(closure1_nonroot):
                missing = sorted(direct1_nonroot - closure1_nonroot)
                raise ValueError(
                    f"t1 direct terms are absent from t1 closure for "
                    f"{protein_id}/{aspect}: {missing[:5]}"
                )
            direct0_nonroot_by_aspect[aspect] = direct0_nonroot
            direct1_nonroot_by_aspect[aspect] = direct1_nonroot
            closure0_nonroot_by_aspect[aspect] = closure0_nonroot
            closure1_nonroot_by_aspect[aspect] = closure1_nonroot
            if root in direct0:
                direct0_roots.add(root)

        global_direct0 = set().union(*direct0_nonroot_by_aspect.values())
        global_direct1 = set().union(*direct1_nonroot_by_aspect.values())
        t0_available = protein_id in t0_present
        t1_available = protein_id in t1_present
        if not t0_available:
            global_state = "unknown"
        elif global_direct0:
            global_state = "known_qualifying"
        elif direct0_roots:
            global_state = "root_only"
        else:
            global_state = "no_qualifying"

        for aspect in ASPECTS:
            direct0 = direct0_nonroot_by_aspect[aspect]
            direct1 = direct1_nonroot_by_aspect[aspect]
            closure0 = closure0_nonroot_by_aspect[aspect]
            closure1 = closure1_nonroot_by_aspect[aspect]
            retained = closure0 & closure1
            gained = closure1 - closure0
            lost = closure0 - closure1
            root_only_aspect = (
                ASPECT_TO_ROOT[aspect]
                in set(direct0_by_aspect.get(aspect, frozenset()))
                and not direct0
            )
            if not t0_available:
                aspect_state = "unknown"
            elif direct0:
                aspect_state = "same_aspect_partial"
            elif global_direct0:
                aspect_state = "cross_ontology_known"
            elif root_only_aspect or direct0_roots:
                aspect_state = "root_only"
            else:
                aspect_state = "no_qualifying"

            rows.append(
                {
                    "protein_id": protein_id,
                    "aspect": aspect,
                    "t0_state_available": t0_available,
                    "t1_state_available": t1_available,
                    "global_knowledge_state": global_state,
                    "aspect_knowledge_state": aspect_state,
                    "cross_ontology_known": aspect_state == "cross_ontology_known",
                    "same_aspect_partial": aspect_state == "same_aspect_partial",
                    "root_only_t0": root_only_aspect,
                    "direct_t0_terms": tuple(sorted(direct0)),
                    "direct_t1_terms": tuple(sorted(direct1)),
                    "closure_t0_terms": tuple(sorted(closure0)),
                    "closure_t1_terms": tuple(sorted(closure1)),
                    "retained_terms": tuple(sorted(retained)),
                    "gained_terms": tuple(sorted(gained)),
                    "lost_terms": tuple(sorted(lost)),
                    "direct_t0_count": len(direct0),
                    "direct_t1_count": len(direct1),
                    "closure_t0_count": len(closure0),
                    "closure_t1_count": len(closure1),
                    "retained_count": len(retained),
                    "gained_count": len(gained),
                    "lost_count": len(lost),
                    "global_direct_t0_count": len(global_direct0),
                    "global_direct_t1_count": len(global_direct1),
                    "global_t0_empty": t0_available and not global_direct0,
                    "aspect_t0_empty": t0_available and not direct0,
                    "aspect_has_gain": bool(gained),
                    "aspect_has_loss": bool(lost),
                }
            )
    return rows


def build_term_transition_rows(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    transitions: list[dict[str, Any]] = []
    for row in rows:
        terms = (
            set(row["closure_t0_terms"])
            | set(row["closure_t1_terms"])
            | set(row["direct_t0_terms"])
            | set(row["direct_t1_terms"])
        )
        for term in sorted(terms):
            closure_t0 = term in row["closure_t0_terms"]
            closure_t1 = term in row["closure_t1_terms"]
            if closure_t0 and closure_t1:
                transition = "retained_known"
            elif closure_t1:
                transition = "gained"
            else:
                transition = "lost"
            transitions.append(
                {
                    "protein_id": row["protein_id"],
                    "aspect": row["aspect"],
                    "go_id": term,
                    "direct_t0": term in row["direct_t0_terms"],
                    "direct_t1": term in row["direct_t1_terms"],
                    "closure_t0": closure_t0,
                    "closure_t1": closure_t1,
                    "transition": transition,
                }
            )
    return transitions


def descriptive_cohort_masks(
    rows: Sequence[Mapping[str, Any]], protein_ids: Sequence[str], aspect: str
) -> dict[str, np.ndarray]:
    if aspect not in ASPECTS:
        raise ValueError(f"Unsupported aspect: {aspect}")
    if len(set(protein_ids)) != len(protein_ids):
        raise ValueError("Cohort protein IDs must be unique")
    indexed = {(str(row["protein_id"]), str(row["aspect"])): row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError("Temporal ledger contains duplicate protein/aspect rows")
    selected = []
    for protein_id in protein_ids:
        key = (protein_id, aspect)
        if key not in indexed:
            raise ValueError(f"Temporal ledger lacks {aspect} row for {protein_id}")
        selected.append(indexed[key])
    return {
        "global_t0_empty": np.asarray(
            [bool(row["global_t0_empty"]) for row in selected], dtype=bool
        ),
        "global_t0_nonempty": np.asarray(
            [not bool(row["global_t0_empty"]) for row in selected], dtype=bool
        ),
        "aspect_t0_empty": np.asarray(
            [bool(row["aspect_t0_empty"]) for row in selected], dtype=bool
        ),
        "aspect_t0_nonempty": np.asarray(
            [not bool(row["aspect_t0_empty"]) for row in selected], dtype=bool
        ),
        "aspect_has_gain": np.asarray(
            [bool(row["aspect_has_gain"]) for row in selected], dtype=bool
        ),
        "aspect_has_loss": np.asarray(
            [bool(row["aspect_has_loss"]) for row in selected], dtype=bool
        ),
        "global_no_qualifying": np.asarray(
            [row.get("global_knowledge_state") == "no_qualifying" for row in selected],
            dtype=bool,
        ),
        "global_known_qualifying": np.asarray(
            [
                row.get("global_knowledge_state") == "known_qualifying"
                for row in selected
            ],
            dtype=bool,
        ),
        "cross_ontology_known": np.asarray(
            [
                row.get("aspect_knowledge_state") == "cross_ontology_known"
                for row in selected
            ],
            dtype=bool,
        ),
        "same_aspect_partial": np.asarray(
            [
                row.get("aspect_knowledge_state") == "same_aspect_partial"
                for row in selected
            ],
            dtype=bool,
        ),
        "unknown_t0_state": np.asarray(
            [row.get("global_knowledge_state") == "unknown" for row in selected],
            dtype=bool,
        ),
    }
