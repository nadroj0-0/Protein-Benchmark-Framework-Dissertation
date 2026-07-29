#!/usr/bin/env python3
"""Build immutable policy-labelled temporal cohort and transition states."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

from label_space_common import (
    atomic_write_json,
    atomic_write_text,
    output_manifest,
    peak_rss_bytes,
    require_unchanged,
    sha256_file,
)
from temporal_annotation_common import (
    ASPECTS,
    build_temporal_state_rows,
    build_term_transition_rows,
    read_annotation_rows,
    read_exposure_rows,
    read_protein_scope,
    validate_term_aspects,
)


def _cohort_tsv(rows: Sequence[Mapping[str, Any]]) -> str:
    fields = (
        "protein_id",
        "aspect",
        "t0_state_available",
        "t1_state_available",
        "global_knowledge_state",
        "aspect_knowledge_state",
        "cross_ontology_known",
        "same_aspect_partial",
        "root_only_t0",
        "direct_t0_terms",
        "direct_t1_terms",
        "closure_t0_terms",
        "closure_t1_terms",
        "retained_terms",
        "gained_terms",
        "lost_terms",
        "direct_t0_count",
        "direct_t1_count",
        "closure_t0_count",
        "closure_t1_count",
        "retained_count",
        "gained_count",
        "lost_count",
        "global_direct_t0_count",
        "global_direct_t1_count",
        "global_t0_empty",
        "aspect_t0_empty",
        "aspect_has_gain",
        "aspect_has_loss",
        "train_id_member",
        "valid_id_member",
        "train_sequence_member",
        "valid_sequence_member",
        "train_homology_cluster_member",
        "modality_availability",
        "feature_temporal_policy",
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        value = dict(row)
        for field in (
            "direct_t0_terms",
            "direct_t1_terms",
            "closure_t0_terms",
            "closure_t1_terms",
            "retained_terms",
            "gained_terms",
            "lost_terms",
        ):
            value[field] = "|".join(value[field])
        for field in (
            "t0_state_available",
            "t1_state_available",
            "cross_ontology_known",
            "same_aspect_partial",
            "root_only_t0",
            "global_t0_empty",
            "aspect_t0_empty",
            "aspect_has_gain",
            "aspect_has_loss",
        ):
            value[field] = int(bool(value[field]))
        for field in (
            "train_id_member",
            "valid_id_member",
            "train_sequence_member",
            "valid_sequence_member",
            "train_homology_cluster_member",
        ):
            value[field] = (
                "unknown" if value[field] is None else int(bool(value[field]))
            )
        writer.writerow(value)
    return stream.getvalue()


def _transition_tsv(rows: Sequence[Mapping[str, Any]]) -> str:
    fields = (
        "protein_id",
        "aspect",
        "go_id",
        "direct_t0",
        "direct_t1",
        "closure_t0",
        "closure_t1",
        "transition",
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        value = dict(row)
        for field in ("direct_t0", "direct_t1", "closure_t0", "closure_t1"):
            value[field] = int(bool(value[field]))
        writer.writerow(value)
    return stream.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t0-direct-annotations", type=Path, required=True)
    parser.add_argument("--t1-direct-annotations", type=Path, required=True)
    parser.add_argument("--t0-closure-annotations", type=Path, required=True)
    parser.add_argument("--t1-closure-annotations", type=Path, required=True)
    parser.add_argument("--t0-protein-presence", type=Path, required=True)
    parser.add_argument("--t1-protein-presence", type=Path, required=True)
    parser.add_argument(
        "--exposure-table",
        type=Path,
        help=(
            "Required when any analyzed protein has qualifying t0 knowledge; "
            "records train/validation ID, sequence, cluster and feature exposure."
        ),
    )
    parser.add_argument("--protein-scope", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--t0-snapshot", required=True)
    parser.add_argument("--t1-snapshot", required=True)
    parser.add_argument("--evidence-policy-id", required=True)
    parser.add_argument("--graph-policy-id", required=True)
    parser.add_argument("--relationship", action="append", default=[], required=True)
    parser.add_argument(
        "--expected-global-knowledge-state",
        choices=("no_qualifying", "known_qualifying", "root_only", "unknown"),
        help="Optional fail-loud assertion for the complete protein scope.",
    )
    parser.add_argument("--benchmark-id", required=True)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")
    if not args.t0_snapshot.strip() or not args.t1_snapshot.strip():
        raise ValueError("Temporal snapshot labels must be non-empty")
    if args.t0_snapshot == args.t1_snapshot:
        raise ValueError("t0 and t1 snapshot labels must differ")
    relationships = tuple(dict.fromkeys(value.strip() for value in args.relationship))
    if any(not value for value in relationships):
        raise ValueError("--relationship values must be non-empty")
    if len(relationships) != len(args.relationship):
        raise ValueError("--relationship values must be unique")
    if not args.evidence_policy_id.strip() or not args.graph_policy_id.strip():
        raise ValueError("Evidence and graph policy IDs must be non-empty")

    protein_ids, scope_contract = read_protein_scope(args.protein_scope.resolve())
    scope = set(protein_ids)
    t0_presence_ids, t0_presence_contract = read_protein_scope(
        args.t0_protein_presence.resolve()
    )
    t1_presence_ids, t1_presence_contract = read_protein_scope(
        args.t1_protein_presence.resolve()
    )
    t0_direct, t0_direct_contract = read_annotation_rows(
        args.t0_direct_annotations.resolve(), scope
    )
    t1_direct, t1_direct_contract = read_annotation_rows(
        args.t1_direct_annotations.resolve(), scope
    )
    t0_closure, t0_closure_contract = read_annotation_rows(
        args.t0_closure_annotations.resolve(), scope
    )
    t1_closure, t1_closure_contract = read_annotation_rows(
        args.t1_closure_annotations.resolve(), scope
    )
    for left, right in (
        (t0_direct_contract, t1_direct_contract),
        (t0_direct_contract, t0_closure_contract),
        (t1_direct_contract, t1_closure_contract),
        (t0_closure_contract, t1_closure_contract),
    ):
        validate_term_aspects(left, right)
    rows = build_temporal_state_rows(
        protein_ids,
        t0_direct,
        t1_direct,
        t0_closure,
        t1_closure,
        set(t0_presence_ids),
        set(t1_presence_ids),
    )
    exposure_contract: dict[str, Any] | None = None
    exposure_rows: dict[str, dict[str, Any]] = {}
    if args.exposure_table is not None:
        exposure_rows, exposure_contract = read_exposure_rows(
            args.exposure_table.resolve(), scope
        )
    known_proteins = {
        str(row["protein_id"])
        for row in rows
        if row["global_knowledge_state"] == "known_qualifying"
    }
    missing_exposure = sorted(known_proteins - set(exposure_rows))
    if missing_exposure:
        raise ValueError(
            "Known-protein cohort analysis requires exposure rows; missing "
            f"{len(missing_exposure)} proteins: {missing_exposure[:5]}"
        )
    for row in rows:
        exposure = exposure_rows.get(str(row["protein_id"]))
        for field in (
            "train_id_member",
            "valid_id_member",
            "train_sequence_member",
            "valid_sequence_member",
            "train_homology_cluster_member",
        ):
            row[field] = exposure[field] if exposure is not None else None
        row["modality_availability"] = (
            exposure["modality_availability"]
            if exposure is not None
            else "not_required_global_no_qualifying"
        )
        row["feature_temporal_policy"] = (
            exposure["feature_temporal_policy"]
            if exposure is not None
            else "not_required_global_no_qualifying"
        )
    transitions = build_term_transition_rows(rows)
    protein_states: dict[str, str] = {}
    for row in rows:
        previous = protein_states.setdefault(
            str(row["protein_id"]), str(row["global_knowledge_state"])
        )
        if previous != row["global_knowledge_state"]:
            raise ValueError(
                f"Global knowledge state differs by aspect for {row['protein_id']}"
            )
    if args.expected_global_knowledge_state is not None:
        unexpected = sorted(
            protein_id
            for protein_id, state in protein_states.items()
            if state != args.expected_global_knowledge_state
        )
        if unexpected:
            raise ValueError(
                f"{len(unexpected)} proteins violate expected global state "
                f"{args.expected_global_knowledge_state}: {unexpected[:5]}"
            )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.stage-", dir=str(output_dir.parent)
        )
    )
    started = time.perf_counter()
    try:
        summary = {
            "schema_version": 1,
            "status": "complete",
            "analysis_kind": "temporal_annotation_state_ledger",
            "scientific_label": "policy_bound_temporal_state_inventory",
            "benchmark_id": args.benchmark_id,
            "t0_snapshot": args.t0_snapshot,
            "t1_snapshot": args.t1_snapshot,
            "evidence_policy_id": args.evidence_policy_id,
            "graph_policy_id": args.graph_policy_id,
            "relationships": list(relationships),
            "expected_global_knowledge_state": (args.expected_global_knowledge_state),
            "policy_boundary": (
                "this Layer-B tool consumes separately normalized direct and "
                "separately propagated t0/t1 inputs; it performs no GAF parsing, "
                "GO resolution, evidence filtering, or propagation"
            ),
            "protein_count": len(protein_ids),
            "aspects": list(ASPECTS),
            "ledger_rows": len(rows),
            "term_transition_rows": len(transitions),
            "inputs": {
                "protein_scope": scope_contract,
                "t0_protein_presence": t0_presence_contract,
                "t1_protein_presence": t1_presence_contract,
                "exposure_table": exposure_contract,
                "t0_direct_annotations": {
                    key: value
                    for key, value in t0_direct_contract.items()
                    if key != "term_aspects"
                },
                "t1_direct_annotations": {
                    key: value
                    for key, value in t1_direct_contract.items()
                    if key != "term_aspects"
                },
                "t0_closure_annotations": {
                    key: value
                    for key, value in t0_closure_contract.items()
                    if key != "term_aspects"
                },
                "t1_closure_annotations": {
                    key: value
                    for key, value in t1_closure_contract.items()
                    if key != "term_aspects"
                },
            },
            "cohort_counts": {
                aspect: {
                    state: sum(
                        row["aspect_knowledge_state"] == state
                        for row in rows
                        if row["aspect"] == aspect
                    )
                    for state in (
                        "no_qualifying",
                        "cross_ontology_known",
                        "same_aspect_partial",
                        "root_only",
                        "unknown",
                    )
                }
                for aspect in ASPECTS
            },
            "global_knowledge_counts": {
                state: sum(value == state for value in protein_states.values())
                for state in (
                    "no_qualifying",
                    "known_qualifying",
                    "root_only",
                    "unknown",
                )
            },
            "transition_counts": {
                transition: sum(row["transition"] == transition for row in transitions)
                for transition in ("retained_known", "gained", "lost")
            },
        }
        atomic_write_text(stage / "protein_cohorts.tsv", _cohort_tsv(rows))
        atomic_write_text(
            stage / "protein_term_states.tsv", _transition_tsv(transitions)
        )
        summary["resource_usage"] = {
            "wall_seconds": time.perf_counter() - started,
            "peak_rss_bytes": peak_rss_bytes(),
        }
        atomic_write_json(stage / "temporal_annotation_ledger.json", summary)

        for path, contract, label in (
            (args.protein_scope.resolve(), scope_contract, "Protein scope"),
            (
                args.t0_protein_presence.resolve(),
                t0_presence_contract,
                "t0 protein presence",
            ),
            (
                args.t1_protein_presence.resolve(),
                t1_presence_contract,
                "t1 protein presence",
            ),
            (
                args.t0_direct_annotations.resolve(),
                t0_direct_contract,
                "t0 direct annotations",
            ),
            (
                args.t1_direct_annotations.resolve(),
                t1_direct_contract,
                "t1 direct annotations",
            ),
            (
                args.t0_closure_annotations.resolve(),
                t0_closure_contract,
                "t0 closure annotations",
            ),
            (
                args.t1_closure_annotations.resolve(),
                t1_closure_contract,
                "t1 closure annotations",
            ),
        ):
            require_unchanged(path, contract, label)
        if args.exposure_table is not None and exposure_contract is not None:
            require_unchanged(
                args.exposure_table.resolve(),
                exposure_contract,
                "Exposure table",
            )
        artifacts = output_manifest(
            stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"}
        )
        atomic_write_json(stage / "output_manifest.json", artifacts)
        atomic_write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "benchmark_id": args.benchmark_id,
                "output_manifest_sha256": sha256_file(stage / "output_manifest.json"),
            },
        )
        os.replace(stage, output_dir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
