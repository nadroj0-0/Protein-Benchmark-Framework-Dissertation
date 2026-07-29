#!/usr/bin/env python3
"""Evaluate policy-bound PFP knowledge cohorts and retained/gained term partitions."""

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

import numpy as np

from evaluate_pfp_information_content import (
    bootstrap_fixed_threshold,
    read_information_accretion,
    threshold_metrics,
)
from label_space_common import (
    atomic_write_json,
    atomic_write_text,
    file_snapshot,
    output_manifest,
    peak_rss_bytes,
    require_unchanged,
    sha256_file,
)
from pfp_sensitivity_common import (
    load_aspect_bundle,
    require_evaluation_split,
    verify_artifact_manifest,
)


COHORTS = (
    "global_no_qualifying",
    "global_known_qualifying",
    "cross_ontology_known",
    "same_aspect_partial",
    "root_only_t0",
    "unknown_t0_state",
)


def _tsv(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _verify_ledger(root: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    summary_path = root / "temporal_annotation_ledger.json"
    cohorts_path = root / "protein_cohorts.tsv"
    output_path = root / "output_manifest.json"
    complete_path = root / "RUN_COMPLETE.json"
    snapshots = {
        path: file_snapshot(path)
        for path in (summary_path, cohorts_path, output_path, complete_path)
    }
    complete = json.loads(complete_path.read_text(encoding="utf-8"))
    if not complete.get("complete"):
        raise ValueError(f"Temporal ledger is incomplete: {root}")
    if complete.get("output_manifest_sha256") != sha256_file(output_path):
        raise ValueError("Temporal ledger completion marker does not bind its manifest")
    manifest = json.loads(output_path.read_text(encoding="utf-8"))
    listed = {value["path"]: value for value in manifest.get("files", [])}
    for path in (summary_path, cohorts_path):
        value = listed.get(path.name)
        if (
            value is None
            or value.get("bytes") != path.stat().st_size
            or value.get("sha256") != sha256_file(path)
        ):
            raise ValueError(f"Temporal ledger manifest does not bind {path.name}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (
        summary.get("status") != "complete"
        or summary.get("analysis_kind") != "temporal_annotation_state_ledger"
    ):
        raise ValueError("Unsupported temporal ledger summary")
    with cohorts_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t", strict=True))
    if not rows:
        raise ValueError("Temporal cohort table is empty")
    for path, snapshot in snapshots.items():
        require_unchanged(path, snapshot, "Temporal ledger artifact")
    return summary, rows


def _parse_terms(value: str) -> set[str]:
    return {term for term in value.split("|") if term}


def _cohort_mask(rows: Sequence[Mapping[str, str]], cohort: str) -> np.ndarray:
    if cohort == "global_no_qualifying":
        return np.asarray(
            [row["global_knowledge_state"] == "no_qualifying" for row in rows],
            dtype=bool,
        )
    if cohort == "global_known_qualifying":
        return np.asarray(
            [row["global_knowledge_state"] == "known_qualifying" for row in rows],
            dtype=bool,
        )
    if cohort == "cross_ontology_known":
        return np.asarray(
            [row["aspect_knowledge_state"] == "cross_ontology_known" for row in rows],
            dtype=bool,
        )
    if cohort == "same_aspect_partial":
        return np.asarray(
            [row["aspect_knowledge_state"] == "same_aspect_partial" for row in rows],
            dtype=bool,
        )
    if cohort == "root_only_t0":
        return np.asarray(
            [row["aspect_knowledge_state"] == "root_only" for row in rows],
            dtype=bool,
        )
    if cohort == "unknown_t0_state":
        return np.asarray(
            [row["global_knowledge_state"] == "unknown" for row in rows],
            dtype=bool,
        )
    raise ValueError(f"Unknown cohort: {cohort}")


def _matrix_from_terms(
    rows: Sequence[Mapping[str, str]],
    field: str,
    term_index: Mapping[str, int],
) -> tuple[np.ndarray, int]:
    matrix = np.zeros((len(rows), len(term_index)), dtype=np.uint8)
    outside = 0
    for protein_index, row in enumerate(rows):
        for term in _parse_terms(row[field]):
            index = term_index.get(term)
            if index is None:
                outside += 1
            else:
                matrix[protein_index, index] = 1
    return matrix, outside


def _metric_rows(
    *,
    benchmark_id: str,
    mode: str,
    aspect: str,
    cohort: str,
    target_component: str,
    truth: np.ndarray,
    scores: np.ndarray,
    weights: np.ndarray,
    threshold: float,
    threshold_type: str,
    bootstrap_replicates: int,
    bootstrap_seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    population = truth.shape[0]
    if population == 0:
        return (
            [
                {
                    "benchmark_id": benchmark_id,
                    "mode": mode,
                    "aspect": aspect,
                    "cohort": cohort,
                    "target_component": target_component,
                    "status": "not_evaluable_empty_cohort",
                    "threshold_type": threshold_type,
                    "threshold": threshold,
                    "population_proteins": 0,
                }
            ],
            [],
        )
    if not truth.any():
        return (
            [
                {
                    "benchmark_id": benchmark_id,
                    "mode": mode,
                    "aspect": aspect,
                    "cohort": cohort,
                    "target_component": target_component,
                    "status": "not_evaluable_no_positive_targets",
                    "threshold_type": threshold_type,
                    "threshold": threshold,
                    "population_proteins": population,
                }
            ],
            [],
        )
    metrics = threshold_metrics(truth, scores, weights, threshold)
    row = {
        "benchmark_id": benchmark_id,
        "mode": mode,
        "aspect": aspect,
        "cohort": cohort,
        "target_component": target_component,
        "status": "complete",
        "threshold_type": threshold_type,
        "precision_interpretation": (
            "known_candidate_precision"
            if target_component == "retained_known"
            else (
                "novel_precision"
                if target_component.startswith("gained_novel")
                else "protein_centric_precision"
            )
        ),
        **metrics,
    }
    intervals = bootstrap_fixed_threshold(
        truth,
        scores,
        weights,
        threshold,
        bootstrap_replicates,
        bootstrap_seed,
    )
    return [row], [
        {
            "benchmark_id": benchmark_id,
            "mode": mode,
            "aspect": aspect,
            "cohort": cohort,
            "target_component": target_component,
            "threshold_type": threshold_type,
            "threshold": threshold,
            **interval,
        }
        for interval in intervals
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-manifest", type=Path, required=True)
    parser.add_argument("--temporal-ledger-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--aspect", action="append", default=[])
    parser.add_argument(
        "--truth-graph-policy-id",
        required=True,
        help="Must exactly match the temporal ledger graph policy.",
    )
    parser.add_argument(
        "--fixed-threshold",
        action="append",
        default=[],
        help="Optional ASPECT=VALUE validation-fixed threshold.",
    )
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=0)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")
    if args.bootstrap_replicates < 0:
        raise ValueError("--bootstrap-replicates must be non-negative")
    threshold_overrides: dict[str, float] = {}
    for value in args.fixed_threshold:
        fields = value.split("=", 1)
        if len(fields) != 2 or fields[0] in threshold_overrides:
            raise ValueError(f"Invalid or duplicate --fixed-threshold: {value}")
        threshold = float(fields[1])
        if not 0 <= threshold <= 1:
            raise ValueError(f"Threshold is outside [0,1]: {value}")
        threshold_overrides[fields[0]] = threshold

    prediction_path = args.prediction_manifest.resolve()
    prediction_sha = sha256_file(prediction_path)
    manifest, artifact_root = verify_artifact_manifest(prediction_path)
    require_evaluation_split(manifest, "test", "Knowledge-cohort analysis")
    ledger_root = args.temporal_ledger_dir.resolve()
    ledger_summary, ledger_rows = _verify_ledger(ledger_root)
    if ledger_summary["benchmark_id"] != manifest["benchmark_id"]:
        raise ValueError("Temporal ledger and prediction benchmark IDs differ")
    if ledger_summary["graph_policy_id"] != args.truth_graph_policy_id:
        raise ValueError(
            "Declared truth graph policy differs from temporal ledger graph policy"
        )
    available = list(manifest["selected_aspects"])
    aspects = list(args.aspect) or available
    if len(aspects) != len(set(aspects)) or not set(aspects).issubset(available):
        raise ValueError("Selected aspects are invalid or duplicated")
    unknown_thresholds = sorted(set(threshold_overrides) - set(aspects))
    if unknown_thresholds:
        raise ValueError(
            f"Threshold overrides target unselected aspects: {unknown_thresholds}"
        )

    ledger_index = {(row["protein_id"], row["aspect"]): row for row in ledger_rows}
    if len(ledger_index) != len(ledger_rows):
        raise ValueError("Temporal ledger has duplicate protein/aspect rows")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.stage-", dir=str(output_dir.parent)
        )
    )
    metric_output: list[dict[str, Any]] = []
    bootstrap_output: list[dict[str, Any]] = []
    copy_output: list[dict[str, Any]] = []
    aspect_reports: dict[str, Any] = {}
    started = time.perf_counter()
    try:
        for aspect_number, aspect in enumerate(aspects):
            bundle = load_aspect_bundle(manifest, artifact_root, aspect)
            selected_rows = []
            for protein_id in bundle["protein_ids"]:
                key = (protein_id, aspect)
                if key not in ledger_index:
                    raise ValueError(
                        f"Temporal ledger lacks {aspect} state for {protein_id}"
                    )
                selected_rows.append(ledger_index[key])
            term_ids = [
                term
                for index, term in enumerate(bundle["go_terms"])
                if index != bundle["root_index"]
            ]
            term_index = {term: index for index, term in enumerate(term_ids)}
            truth = bundle["truth"][
                :,
                [
                    index
                    for index in range(bundle["truth"].shape[1])
                    if index != bundle["root_index"]
                ],
            ]
            scores = bundle["scores"][
                :,
                [
                    index
                    for index in range(bundle["scores"].shape[1])
                    if index != bundle["root_index"]
                ],
            ]
            ia, ia_contract = read_information_accretion(
                bundle["ia_path"], bundle["go_terms"]
            )
            weights = ia[
                [index for index in range(len(ia)) if index != bundle["root_index"]]
            ]
            closure0, outside_t0 = _matrix_from_terms(
                selected_rows, "closure_t0_terms", term_index
            )
            closure1, outside_t1 = _matrix_from_terms(
                selected_rows, "closure_t1_terms", term_index
            )
            if not np.array_equal(closure1, truth):
                mismatch = int(np.count_nonzero(closure1 != truth))
                raise ValueError(
                    f"Temporal t1 closure differs from prediction truth for "
                    f"{aspect}: {mismatch} cells"
                )
            retained = np.logical_and(closure0, closure1).astype(np.uint8)
            gained = np.logical_and(closure1, ~closure0.astype(bool)).astype(np.uint8)
            lost = np.logical_and(closure0, ~closure1.astype(bool)).astype(np.uint8)

            if aspect in threshold_overrides:
                threshold = threshold_overrides[aspect]
                threshold_type = "validation_fixed_operating_threshold"
            else:
                threshold = float(
                    bundle["specification"]["canonical_cafa_metrics"]["threshold"]
                )
                threshold_type = "descriptive_test_oracle_fixed"

            cohort_counts: dict[str, int] = {}
            for cohort_number, cohort in enumerate(COHORTS):
                mask = _cohort_mask(selected_rows, cohort)
                cohort_counts[cohort] = int(mask.sum())
                rows, intervals = _metric_rows(
                    benchmark_id=manifest["benchmark_id"],
                    mode=manifest["mode"],
                    aspect=aspect,
                    cohort=cohort,
                    target_component="combined_t1_truth",
                    truth=truth[mask],
                    scores=scores[mask],
                    weights=weights,
                    threshold=threshold,
                    threshold_type=threshold_type,
                    bootstrap_replicates=args.bootstrap_replicates,
                    bootstrap_seed=(
                        args.bootstrap_seed + 1000 * aspect_number + 10 * cohort_number
                    ),
                )
                metric_output.extend(rows)
                bootstrap_output.extend(intervals)

            partial_mask = _cohort_mask(selected_rows, "same_aspect_partial")
            known_scores = np.where(
                closure0[partial_mask].astype(bool), scores[partial_mask], -1.0
            )
            gained_scores = np.where(
                ~closure0[partial_mask].astype(bool), scores[partial_mask], -1.0
            )
            for component_number, (
                component,
                component_truth,
                component_scores,
            ) in enumerate(
                (
                    ("retained_known", retained[partial_mask], known_scores),
                    (
                        "gained_novel_acquisition_conditioned",
                        gained[partial_mask],
                        gained_scores,
                    ),
                    (
                        "gained_novel_deployment_like",
                        gained[partial_mask],
                        gained_scores,
                    ),
                )
            ):
                component_mask = (
                    component_truth.any(axis=1)
                    if component == "gained_novel_acquisition_conditioned"
                    else np.ones(component_truth.shape[0], dtype=bool)
                )
                rows, intervals = _metric_rows(
                    benchmark_id=manifest["benchmark_id"],
                    mode=manifest["mode"],
                    aspect=aspect,
                    cohort="same_aspect_partial",
                    target_component=component,
                    truth=component_truth[component_mask],
                    scores=component_scores[component_mask],
                    weights=weights,
                    threshold=threshold,
                    threshold_type=threshold_type,
                    bootstrap_replicates=args.bootstrap_replicates,
                    bootstrap_seed=(
                        args.bootstrap_seed
                        + 1000 * aspect_number
                        + 100
                        + component_number
                    ),
                )
                metric_output.extend(rows)
                bootstrap_output.extend(intervals)

            if partial_mask.any():
                copy_scores = closure0[partial_mask].astype(np.float64)
                copy_metrics = threshold_metrics(
                    truth[partial_mask], copy_scores, weights, 0.5
                )
                copy_output.append(
                    {
                        "benchmark_id": manifest["benchmark_id"],
                        "mode": "t0_annotation_copy_baseline",
                        "aspect": aspect,
                        "cohort": "same_aspect_partial",
                        "threshold": 0.5,
                        **copy_metrics,
                    }
                )

            aspect_reports[aspect] = {
                "prediction_proteins": len(bundle["protein_ids"]),
                "model_terms_nonroot": len(term_ids),
                "cohort_counts": cohort_counts,
                "closure_t0_terms_outside_model_universe": outside_t0,
                "closure_t1_terms_outside_model_universe": outside_t1,
                "retained_positive_cells": int(retained.sum()),
                "gained_positive_cells": int(gained.sum()),
                "lost_positive_cells": int(lost.sum()),
                "threshold": threshold,
                "threshold_type": threshold_type,
                "ia_sha256": ia_contract["sha256"],
                "partial_cohort_exposure": {
                    field: {
                        "true": sum(
                            row[field] == "1"
                            for row, selected in zip(selected_rows, partial_mask)
                            if selected
                        ),
                        "false": sum(
                            row[field] == "0"
                            for row, selected in zip(selected_rows, partial_mask)
                            if selected
                        ),
                        "unknown": sum(
                            row[field] == "unknown"
                            for row, selected in zip(selected_rows, partial_mask)
                            if selected
                        ),
                    }
                    for field in (
                        "train_id_member",
                        "valid_id_member",
                        "train_sequence_member",
                        "valid_sequence_member",
                        "train_homology_cluster_member",
                    )
                },
            }

        report = {
            "schema_version": 1,
            "status": "complete",
            "analysis_kind": "knowledge_cohort_combined_truth",
            "scientific_label": (
                "seen_protein_annotation_extension_sensitivity"
                if any(
                    value["cohort_counts"]["same_aspect_partial"]
                    for value in aspect_reports.values()
                )
                else "accepted_global_no_knowledge_cohort_inventory"
            ),
            "canonicality_label": "noncanonical_flat_cohort_diagnostic",
            "benchmark_id": manifest["benchmark_id"],
            "mode": manifest["mode"],
            "selected_aspects": aspects,
            "truth_graph_policy_id": args.truth_graph_policy_id,
            "evidence_policy_id": ledger_summary["evidence_policy_id"],
            "prediction_manifest": {
                "path": str(prediction_path),
                "sha256": prediction_sha,
            },
            "temporal_ledger": {
                "path": str(ledger_root),
                "summary_sha256": sha256_file(
                    ledger_root / "temporal_annotation_ledger.json"
                ),
                "cohorts_sha256": sha256_file(ledger_root / "protein_cohorts.tsv"),
            },
            "interpretation": (
                "cohort and retained/gained results are flat diagnostics at one "
                "fixed threshold; same-aspect partial results are seen-protein "
                "annotation-extension sensitivity unless an independent exposure "
                "audit proves strict holdout"
            ),
            "aspects": aspect_reports,
            "resource_usage": {
                "wall_seconds": time.perf_counter() - started,
                "peak_rss_bytes": peak_rss_bytes(),
            },
        }
        atomic_write_json(stage / "knowledge_cohort_analysis.json", report)
        metric_fields = (
            "benchmark_id",
            "mode",
            "aspect",
            "cohort",
            "target_component",
            "status",
            "threshold_type",
            "precision_interpretation",
            "threshold",
            "population_proteins",
            "target_positive_proteins",
            "predicted_coverage",
            "precision",
            "recall",
            "f",
            "jaccard_set_agreement",
            "micro_precision",
            "micro_recall",
            "micro_f",
            "weighted_precision",
            "weighted_recall",
            "weighted_f",
            "weighted_jaccard_set_agreement",
            "micro_weighted_precision",
            "micro_weighted_recall",
            "micro_weighted_f",
        )
        atomic_write_text(
            stage / "cohort_metrics.tsv", _tsv(metric_output, metric_fields)
        )
        bootstrap_fields = (
            "benchmark_id",
            "mode",
            "aspect",
            "cohort",
            "target_component",
            "threshold_type",
            "threshold",
            "metric_name",
            "estimate",
            "ci_low",
            "ci_high",
            "bootstrap_replicates",
            "bootstrap_seed",
            "resampling_unit",
            "method",
        )
        atomic_write_text(
            stage / "bootstrap_intervals.tsv",
            _tsv(bootstrap_output, bootstrap_fields),
        )
        copy_fields = (
            "benchmark_id",
            "mode",
            "aspect",
            "cohort",
            "threshold",
            "population_proteins",
            "target_positive_proteins",
            "predicted_coverage",
            "precision",
            "recall",
            "f",
            "jaccard_set_agreement",
            "micro_precision",
            "micro_recall",
            "micro_f",
            "weighted_precision",
            "weighted_recall",
            "weighted_f",
            "weighted_jaccard_set_agreement",
            "micro_weighted_precision",
            "micro_weighted_recall",
            "micro_weighted_f",
        )
        atomic_write_text(
            stage / "copy_baseline_metrics.tsv", _tsv(copy_output, copy_fields)
        )

        if sha256_file(prediction_path) != prediction_sha:
            raise ValueError("Prediction manifest changed during cohort analysis")
        artifacts = output_manifest(
            stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"}
        )
        atomic_write_json(stage / "output_manifest.json", artifacts)
        atomic_write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "analysis_kind": "knowledge_cohort_combined_truth",
                "benchmark_id": manifest["benchmark_id"],
                "mode": manifest["mode"],
                "output_manifest_sha256": sha256_file(stage / "output_manifest.json"),
            },
        )
        os.replace(stage, output_dir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
