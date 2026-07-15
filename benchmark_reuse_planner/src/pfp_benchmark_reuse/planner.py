from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Dict, Sequence, Set, Tuple

from .models import (
    BenchmarkData,
    EmbeddedProtein,
    PlanRecord,
    ProteinRecord,
    REGENERATE_MODALITIES,
    ReusePlan,
    ReusePlannerError,
)


class PlanningError(ReusePlannerError):
    pass


@dataclass
class _KnownProtein:
    sequence: str
    sequence_sha256: str
    benchmarks: Set[str] = field(default_factory=set)
    memberships: Set[str] = field(default_factory=set)


def build_plan(
    embedded_benchmarks: Sequence[BenchmarkData], target_benchmark: BenchmarkData
) -> ReusePlan:
    if not embedded_benchmarks:
        raise PlanningError("At least one embedded benchmark is required")

    ordered_embedded = tuple(sorted(embedded_benchmarks, key=lambda item: item.name))
    _validate_unique_names(ordered_embedded, target_benchmark)
    known = _build_known_union(ordered_embedded)
    known_records = _known_records(known)
    records = tuple(
        _expected_record(protein_id, target, known)
        for protein_id, target in sorted(target_benchmark.proteins.items())
    )

    plan = ReusePlan(
        embedded_benchmarks=ordered_embedded,
        target_benchmark=target_benchmark,
        known_embedded_proteins=known_records,
        records=records,
    )
    validate_plan(plan)
    return plan


def validate_plan(plan: ReusePlan) -> None:
    if not plan.embedded_benchmarks:
        raise PlanningError("At least one embedded benchmark is required")
    ordered_embedded = tuple(
        sorted(plan.embedded_benchmarks, key=lambda item: item.name)
    )
    if plan.embedded_benchmarks != ordered_embedded:
        raise PlanningError("Embedded benchmarks must be sorted by name")
    _validate_unique_names(ordered_embedded, plan.target_benchmark)
    for benchmark in (*ordered_embedded, plan.target_benchmark):
        _validate_benchmark_records(benchmark)

    target_ids = set(plan.target_benchmark.proteins)
    record_ids = [record.protein_id for record in plan.records]
    if len(record_ids) != len(set(record_ids)):
        raise PlanningError("The action partition contains duplicate target protein IDs")
    if set(record_ids) != target_ids:
        raise PlanningError("The action partition does not equal the target protein population")

    reuse_ids = {record.protein_id for record in plan.records if record.action == "reuse"}
    regenerate_ids = {
        record.protein_id for record in plan.records if record.action == "regenerate"
    }
    if any(record.action not in {"reuse", "regenerate"} for record in plan.records):
        raise PlanningError("Every target protein action must be reuse or regenerate")
    if reuse_ids & regenerate_ids:
        raise PlanningError("The reuse and regenerate partitions overlap")
    if reuse_ids | regenerate_ids != target_ids:
        raise PlanningError("The reuse and regenerate partitions are incomplete")

    for record in plan.records:
        expected_modalities = REGENERATE_MODALITIES if record.action == "regenerate" else ()
        if record.regenerate_modalities != expected_modalities:
            raise PlanningError("Regenerate modality scheduling is inconsistent")

    known = _build_known_union(ordered_embedded)
    expected_known = _known_records(known)
    if plan.known_embedded_proteins != expected_known:
        raise PlanningError("Known embedded proteins do not match the embedded union")
    expected_records = tuple(
        _expected_record(protein_id, target, known)
        for protein_id, target in sorted(plan.target_benchmark.proteins.items())
    )
    if plan.records != expected_records:
        raise PlanningError(
            "Action records do not exactly match target sequences, memberships, hashes, and policy"
        )


def _validate_unique_names(
    embedded_benchmarks: Sequence[BenchmarkData], target_benchmark: BenchmarkData
) -> None:
    names = [benchmark.name for benchmark in embedded_benchmarks] + [
        target_benchmark.name
    ]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise PlanningError("Benchmark names must be unique: %s" % ", ".join(duplicates))


def _validate_benchmark_records(benchmark: BenchmarkData) -> None:
    for protein_id, protein in benchmark.proteins.items():
        if protein_id != protein.protein_id:
            raise PlanningError("Benchmark protein mapping key does not match its record ID")
        expected_hash = hashlib.sha256(protein.sequence.encode("utf-8")).hexdigest()
        if protein.sequence_sha256 != expected_hash:
            raise PlanningError(
                "Benchmark protein sequence_sha256 does not match its exact sequence: %s"
                % protein_id
            )
        if protein.memberships != tuple(sorted(set(protein.memberships))):
            raise PlanningError("Benchmark protein memberships are not sorted and unique")


def _build_known_union(
    embedded_benchmarks: Sequence[BenchmarkData],
) -> Dict[str, _KnownProtein]:
    known: Dict[str, _KnownProtein] = {}
    for benchmark in embedded_benchmarks:
        for protein_id, protein in sorted(benchmark.proteins.items()):
            existing = known.get(protein_id)
            if existing is not None and existing.sequence != protein.sequence:
                conflicting = sorted(existing.benchmarks | {benchmark.name})
                raise PlanningError(
                    "Embedded protein ID %s has conflicting sequences across benchmarks: %s"
                    % (protein_id, ", ".join(conflicting))
                )
            if existing is None:
                existing = _KnownProtein(
                    sequence=protein.sequence,
                    sequence_sha256=protein.sequence_sha256,
                )
                known[protein_id] = existing
            existing.benchmarks.add(benchmark.name)
            existing.memberships.update(
                "%s:%s" % (benchmark.name, membership)
                for membership in protein.memberships
            )
    return known


def _known_records(known: Dict[str, _KnownProtein]) -> Tuple[EmbeddedProtein, ...]:
    return tuple(
        EmbeddedProtein(
            protein_id=protein_id,
            sequence=protein.sequence,
            sequence_sha256=protein.sequence_sha256,
            embedded_benchmarks=tuple(sorted(protein.benchmarks)),
            embedded_benchmark_memberships=tuple(sorted(protein.memberships)),
        )
        for protein_id, protein in sorted(known.items())
    )


def _expected_record(
    protein_id: str,
    target: ProteinRecord,
    known: Dict[str, _KnownProtein],
) -> PlanRecord:
    target_sequence = target.sequence
    embedded = known.get(protein_id)
    if embedded is None:
        action = "regenerate"
        reason = "protein-id-absent"
        matches: Tuple[str, ...] = ()
        embedded_memberships: Tuple[str, ...] = ()
    elif embedded.sequence != target_sequence:
        action = "regenerate"
        reason = "sequence-mismatch"
        matches = ()
        embedded_memberships = tuple(sorted(embedded.memberships))
    else:
        action = "reuse"
        reason = "exact-id-sequence-match"
        matches = tuple(sorted(embedded.benchmarks))
        embedded_memberships = tuple(sorted(embedded.memberships))
    return PlanRecord(
        protein_id=protein_id,
        sequence=target_sequence,
        sequence_sha256=target.sequence_sha256,
        action=action,
        reason=reason,
        matching_embedded_benchmarks=matches,
        embedded_benchmark_memberships=embedded_memberships,
        target_memberships=target.memberships,
        regenerate_modalities=(REGENERATE_MODALITIES if action == "regenerate" else ()),
    )
