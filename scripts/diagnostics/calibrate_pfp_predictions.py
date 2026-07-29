#!/usr/bin/env python3
"""Fit validation-only calibrated PFP correctness probabilities and audit test transport."""

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

from calibration_common import (
    CalibrationPolicy,
    apply_calibrator,
    calibration_intercept_slope,
    calibration_metrics,
    expected_calibration_error,
    fit_monotone_hierarchical_calibrator,
    propagate_scores_max,
    reliability_rows,
)
from label_space_common import (
    atomic_write_json,
    atomic_write_text,
    output_manifest,
    peak_rss_bytes,
    read_obo,
    sha256_file,
)
from pfp_sensitivity_common import (
    load_aspect_bundle,
    require_evaluation_split,
    sha256_array,
    verify_artifact_manifest,
)
from specificity_common import (
    SpecificityMeasure,
    assign_specificity_bins,
    read_nonnegative_term_values,
)


def _tsv(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _atomic_savez(path: Path, **arrays: np.ndarray) -> None:
    fd, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".npz", dir=str(path.parent)
    )
    try:
        with os.fdopen(fd, "wb") as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def _selected_aspects(
    requested: Sequence[str],
    validation_manifest: Mapping[str, Any],
    test_manifest: Mapping[str, Any],
) -> list[str]:
    available = list(validation_manifest["selected_aspects"])
    if set(available) != set(test_manifest["selected_aspects"]):
        raise ValueError("Validation and test prediction aspects differ")
    result = list(requested) or available
    if len(result) != len(set(result)) or not set(result).issubset(available):
        raise ValueError("Selected calibration aspects are invalid or duplicated")
    return result


def _validate_manifest_pair(
    validation: Mapping[str, Any], test: Mapping[str, Any]
) -> None:
    for field in ("benchmark_id", "mode", "seed"):
        if validation[field] != test[field]:
            raise ValueError(f"Validation/test prediction {field} differs")
    for field in ("framework_commit", "pfp_commit", "benchmark_fingerprint"):
        if validation["provenance"][field] != test["provenance"][field]:
            raise ValueError(f"Validation/test provenance {field} differs")
    if validation["config"]["sha256"] != test["config"]["sha256"]:
        raise ValueError("Validation/test run configuration differs")
    if validation["obo"]["sha256"] != test["obo"]["sha256"]:
        raise ValueError("Validation/test prediction ontology differs")


def _calibration_metrics_with_slope(
    probabilities: np.ndarray,
    truth: np.ndarray,
    reliability_bin_count: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = reliability_rows(probabilities, truth, reliability_bin_count)
    metrics: dict[str, Any] = {
        **calibration_metrics(probabilities, truth),
        **calibration_intercept_slope(probabilities, truth),
        "ece_equal_count": expected_calibration_error(rows),
        "reliability_bin_count": reliability_bin_count,
    }
    return metrics, rows


def _hierarchy_audit(
    probabilities: np.ndarray,
    go_terms: Sequence[str],
    graph: Any,
    tolerance: float,
) -> dict[str, Any]:
    term_index = {term: index for index, term in enumerate(go_terms)}
    edges = []
    for child, child_term in graph.terms.items():
        child_index = term_index.get(child)
        if child_index is None:
            continue
        for parent in child_term.parents:
            parent_index = term_index.get(parent)
            if parent_index is not None:
                edges.append((child_index, parent_index))
    violations = 0
    affected_edges = 0
    worst = 0.0
    for child_index, parent_index in edges:
        difference = probabilities[:, child_index] - probabilities[:, parent_index]
        selected = difference > tolerance
        count = int(selected.sum())
        if count:
            violations += count
            affected_edges += 1
            worst = max(worst, float(difference[selected].max()))
    return {
        "edges_audited": len(edges),
        "protein_edge_events": probabilities.shape[0] * len(edges),
        "violation_events": violations,
        "affected_edges": affected_edges,
        "maximum_child_minus_parent": worst,
        "tolerance": tolerance,
        "projection": "audit_only",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-prediction-manifest", type=Path, required=True)
    parser.add_argument("--test-prediction-manifest", type=Path, required=True)
    parser.add_argument("--obo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--aspect", action="append", default=[])
    parser.add_argument("--positive-ia-bins", type=int, default=4)
    parser.add_argument("--reliability-bins", type=int, default=10)
    parser.add_argument("--score-clip-epsilon", type=float, default=1e-6)
    parser.add_argument("--bin-l2", type=float, default=5.0)
    parser.add_argument("--term-l2", type=float, default=20.0)
    parser.add_argument("--minimum-bin-positives", type=int, default=20)
    parser.add_argument("--minimum-bin-negatives", type=int, default=20)
    parser.add_argument("--minimum-term-positives", type=int, default=20)
    parser.add_argument("--minimum-term-negatives", type=int, default=20)
    parser.add_argument("--maximum-iterations", type=int, default=200)
    parser.add_argument("--optimizer-tolerance", type=float, default=1e-7)
    parser.add_argument("--protein-chunk-size", type=int, default=256)
    parser.add_argument("--hierarchy-tolerance", type=float, default=1e-12)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")
    if not 1 <= args.positive_ia_bins <= 20:
        raise ValueError("--positive-ia-bins must be between 1 and 20")
    if not 2 <= args.reliability_bins <= 100:
        raise ValueError("--reliability-bins must be between 2 and 100")
    if args.hierarchy_tolerance < 0:
        raise ValueError("--hierarchy-tolerance must be non-negative")

    validation_path = args.validation_prediction_manifest.resolve()
    test_path = args.test_prediction_manifest.resolve()
    validation_sha = sha256_file(validation_path)
    test_sha = sha256_file(test_path)
    validation, validation_root = verify_artifact_manifest(validation_path)
    test, test_root = verify_artifact_manifest(test_path)
    require_evaluation_split(validation, "valid", "Calibration fitting")
    require_evaluation_split(test, "test", "Calibration transport evaluation")
    _validate_manifest_pair(validation, test)
    aspects = _selected_aspects(args.aspect, validation, test)

    obo_path = args.obo.resolve()
    obo_sha = sha256_file(obo_path)
    if obo_sha != validation["obo"]["sha256"]:
        raise ValueError("Calibration OBO hash differs from prediction artifacts")
    graph = read_obo(obo_path)
    policy = CalibrationPolicy(
        score_clip_epsilon=args.score_clip_epsilon,
        bin_l2=args.bin_l2,
        term_l2=args.term_l2,
        minimum_bin_positives=args.minimum_bin_positives,
        minimum_bin_negatives=args.minimum_bin_negatives,
        minimum_term_positives=args.minimum_term_positives,
        minimum_term_negatives=args.minimum_term_negatives,
        maximum_iterations=args.maximum_iterations,
        optimizer_tolerance=args.optimizer_tolerance,
        protein_chunk_size=args.protein_chunk_size,
    )
    policy.validate()

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.stage-", dir=str(output_dir.parent)
        )
    )
    started = time.perf_counter()
    reliability_output: list[dict[str, Any]] = []
    shift_output: list[dict[str, Any]] = []
    hierarchy_output: list[dict[str, Any]] = []
    model_aspects: dict[str, Any] = {}
    prediction_outputs: dict[str, Any] = {}
    try:
        for aspect in aspects:
            valid_bundle = load_aspect_bundle(validation, validation_root, aspect)
            test_bundle = load_aspect_bundle(test, test_root, aspect)
            if valid_bundle["go_terms"] != test_bundle["go_terms"]:
                raise ValueError(f"Validation/test GO-term order differs for {aspect}")
            if (
                valid_bundle["specification"]["checkpoint_sha256"]
                != test_bundle["specification"]["checkpoint_sha256"]
            ):
                raise ValueError(f"Validation/test checkpoint differs for {aspect}")
            if (
                valid_bundle["specification"]["ia_file_sha256"]
                != test_bundle["specification"]["ia_file_sha256"]
            ):
                raise ValueError(f"Validation/test IA differs for {aspect}")
            overlap = sorted(
                set(valid_bundle["protein_ids"]) & set(test_bundle["protein_ids"])
            )
            if overlap:
                raise ValueError(
                    f"Validation/test protein IDs overlap for {aspect}: {overlap[:5]}"
                )

            ia_measure = read_nonnegative_term_values(
                valid_bundle["ia_path"],
                valid_bundle["go_terms"],
                measure_name="information_accretion",
                higher_is_more_specific=True,
                zero_bin_label="zero_ia",
            )
            _, all_assignments = assign_specificity_bins(
                valid_bundle["go_terms"],
                SpecificityMeasure(
                    name="information_accretion",
                    values=ia_measure.values,
                    higher_is_more_specific=True,
                    zero_bin_label="zero_ia",
                    source=ia_measure.source,
                ),
                args.positive_ia_bins,
                bin_prefix="positive_ia_q",
                excluded_indices=(valid_bundle["root_index"],),
                excluded_label="root_excluded",
            )
            nonroot_indices = np.asarray(
                [
                    index
                    for index in range(len(valid_bundle["go_terms"]))
                    if index != valid_bundle["root_index"]
                ],
                dtype=np.int64,
            )
            term_ids = [valid_bundle["go_terms"][index] for index in nonroot_indices]
            term_bins = [all_assignments[index] for index in nonroot_indices]
            if "root_excluded" in term_bins:
                raise RuntimeError("Root entered the calibration candidate universe")

            valid_postprop_all = propagate_scores_max(
                valid_bundle["scores"],
                valid_bundle["go_terms"],
                graph,
                aspect,
            )
            test_postprop_all = propagate_scores_max(
                test_bundle["scores"],
                test_bundle["go_terms"],
                graph,
                aspect,
            )
            valid_postprop = valid_postprop_all[:, nonroot_indices]
            test_postprop = test_postprop_all[:, nonroot_indices]
            valid_truth = valid_bundle["truth"][:, nonroot_indices]
            test_truth = test_bundle["truth"][:, nonroot_indices]
            fitted = fit_monotone_hierarchical_calibrator(
                valid_postprop,
                valid_truth,
                term_ids,
                term_bins,
                policy,
            )
            valid_q, fallback = apply_calibrator(
                valid_postprop,
                term_ids,
                term_bins,
                fitted,
                protein_chunk_size=policy.protein_chunk_size,
            )
            test_q, observed_fallback = apply_calibrator(
                test_postprop,
                term_ids,
                term_bins,
                fitted,
                protein_chunk_size=policy.protein_chunk_size,
            )
            if fallback != observed_fallback:
                raise RuntimeError("Calibration fallback changed between splits")
            calibrated_available = (
                fitted["family"] != "uncalibrated_insufficient_support"
            )

            prediction_path = stage / f"{aspect}_calibration_predictions.npz"
            _atomic_savez(
                prediction_path,
                raw_scores=np.asarray(
                    test_bundle["scores"][:, nonroot_indices], dtype=np.float32
                ),
                postprop_scores=np.asarray(test_postprop, dtype=np.float32),
                calibrated_q=test_q,
                truth=np.asarray(test_truth, dtype=np.uint8),
                protein_ids=np.asarray(test_bundle["protein_ids"], dtype=str),
                go_terms=np.asarray(term_ids, dtype=str),
                information_accretion=np.asarray(
                    ia_measure.values[nonroot_indices], dtype=np.float64
                ),
                ia_bins=np.asarray(term_bins, dtype=str),
                fallback_level=np.asarray(fallback, dtype=str),
            )
            prediction_outputs[aspect] = {
                "path": prediction_path.name,
                "bytes": prediction_path.stat().st_size,
                "sha256": sha256_file(prediction_path),
                "calibrated_q_content_sha256": sha256_array(test_q),
                "shape": list(test_q.shape),
                "prediction_interval_status": (
                    "not_computed_exploratory_point_calibration"
                    if calibrated_available
                    else "not_applicable_uncalibrated_insufficient_support"
                ),
            }

            split_values = [
                ("valid", "raw_postprop", valid_postprop, valid_truth),
                ("test", "raw_postprop", test_postprop, test_truth),
            ]
            if calibrated_available:
                split_values.extend(
                    [
                        ("valid", "calibrated", valid_q, valid_truth),
                        ("test", "calibrated_no_refit", test_q, test_truth),
                    ]
                )
            aspect_metrics: dict[str, Any] = {}
            for split, score_kind, probabilities, target in split_values:
                metrics, reliability = _calibration_metrics_with_slope(
                    probabilities, target, args.reliability_bins
                )
                aspect_metrics[f"{split}:{score_kind}"] = metrics
                shift_output.extend(
                    {
                        "benchmark_id": validation["benchmark_id"],
                        "mode": validation["mode"],
                        "aspect": aspect,
                        "split": split,
                        "score_kind": score_kind,
                        "metric_name": name,
                        "estimate": value,
                    }
                    for name, value in metrics.items()
                )
                reliability_output.extend(
                    {
                        "benchmark_id": validation["benchmark_id"],
                        "mode": validation["mode"],
                        "aspect": aspect,
                        "split": split,
                        "score_kind": score_kind,
                        **row,
                    }
                    for row in reliability
                )

            hierarchy_values = [
                ("valid", "raw_postprop", valid_postprop),
                ("test", "raw_postprop", test_postprop),
            ]
            if calibrated_available:
                hierarchy_values.extend(
                    [
                        ("valid", "calibrated", valid_q),
                        ("test", "calibrated_no_refit", test_q),
                    ]
                )
            for split, score_kind, probabilities in hierarchy_values:
                hierarchy_output.append(
                    {
                        "benchmark_id": validation["benchmark_id"],
                        "mode": validation["mode"],
                        "aspect": aspect,
                        "split": split,
                        "score_kind": score_kind,
                        **_hierarchy_audit(
                            probabilities,
                            term_ids,
                            graph,
                            args.hierarchy_tolerance,
                        ),
                    }
                )

            model_aspects[aspect] = {
                "checkpoint_sha256": valid_bundle["specification"]["checkpoint_sha256"],
                "ia_sha256": valid_bundle["specification"]["ia_file_sha256"],
                "go_terms_sha256": valid_bundle["specification"]["go_terms_sha256"],
                "validation_truth_sha256": valid_bundle["specification"][
                    "truth_content_sha256"
                ],
                "test_truth_sha256": test_bundle["specification"][
                    "truth_content_sha256"
                ],
                "model": fitted,
                "metrics": aspect_metrics,
            }

        model = {
            "schema_version": 1,
            "status": "complete",
            "calibrator_id": (
                f"{validation['benchmark_id']}__{validation['mode']}__"
                "post_selection_validation_calibration"
            ),
            "analysis_label": "post_selection_validation_calibration",
            "target_definition": "benchmark_observed_qualifying_propagated_t1",
            "probability_interpretation": (
                "estimated probability that the term is present in the "
                "benchmark-observed qualifying propagated t1 label set; not "
                "biological truth and not a p-value"
            ),
            "checkpoint_selection_overlap_disclosure": (
                "the same validation population previously influenced checkpoint "
                "selection and early stopping; this is exploratory post-selection "
                "calibration, not an independent calibration holdout"
            ),
            "benchmark_id": validation["benchmark_id"],
            "mode": validation["mode"],
            "selected_aspects": aspects,
            "score_stage": "postprop_is_a_part_of_max",
            "candidate_universe": {
                "terms": "all non-root model-output terms",
                "reporting_floor": 0.0,
                "negative_sampling": "none",
            },
            "provenance": {
                "validation_prediction_manifest": {
                    "path": str(validation_path),
                    "sha256": validation_sha,
                },
                "test_prediction_manifest": {
                    "path": str(test_path),
                    "sha256": test_sha,
                },
                "obo": {"path": str(obo_path), "sha256": obo_sha},
                "framework_commit": validation["provenance"]["framework_commit"],
                "pfp_commit": validation["provenance"]["pfp_commit"],
                "benchmark_fingerprint": validation["provenance"][
                    "benchmark_fingerprint"
                ],
            },
            "policy": {
                **policy.as_dict(),
                "positive_ia_bins": args.positive_ia_bins,
                "reliability_bins": args.reliability_bins,
                "fallback_order": [
                    "term_shrinkage",
                    "aspect_mode_ia_bin",
                    "aspect_mode_platt",
                    "uncalibrated",
                ],
                "p_values": "prohibited",
                "hierarchy_projection": "audit_only",
            },
            "aspects": model_aspects,
            "prediction_outputs": prediction_outputs,
            "uncertainty": {
                "per_prediction_interval": (
                    "not_computed_in_first_exploratory_implementation"
                ),
                "required_for_final_user_facing_tool": (
                    "protein-cluster bootstrap or separately validated interval method"
                ),
            },
        }
        atomic_write_json(stage / "calibration_model.json", model)
        atomic_write_text(
            stage / "calibration_reliability.tsv",
            _tsv(
                reliability_output,
                (
                    "benchmark_id",
                    "mode",
                    "aspect",
                    "split",
                    "score_kind",
                    "reliability_bin",
                    "events",
                    "mean_probability",
                    "observed_fraction",
                    "minimum_probability",
                    "maximum_probability",
                ),
            ),
        )
        atomic_write_text(
            stage / "calibration_shift.tsv",
            _tsv(
                shift_output,
                (
                    "benchmark_id",
                    "mode",
                    "aspect",
                    "split",
                    "score_kind",
                    "metric_name",
                    "estimate",
                ),
            ),
        )
        atomic_write_text(
            stage / "calibration_hierarchy_audit.tsv",
            _tsv(
                hierarchy_output,
                (
                    "benchmark_id",
                    "mode",
                    "aspect",
                    "split",
                    "score_kind",
                    "edges_audited",
                    "protein_edge_events",
                    "violation_events",
                    "affected_edges",
                    "maximum_child_minus_parent",
                    "tolerance",
                    "projection",
                ),
            ),
        )
        analysis = {
            "schema_version": 1,
            "status": "complete",
            "analysis_kind": "post_selection_calibration",
            "scientific_label": "post_selection_validation_calibration",
            "benchmark_id": validation["benchmark_id"],
            "mode": validation["mode"],
            "selected_aspects": aspects,
            "calibration_model_sha256": sha256_file(stage / "calibration_model.json"),
            "resource_usage": {
                "wall_seconds": time.perf_counter() - started,
                "peak_rss_bytes": peak_rss_bytes(),
            },
        }
        atomic_write_json(stage / "calibration_analysis.json", analysis)

        if sha256_file(validation_path) != validation_sha:
            raise ValueError(
                "Validation prediction manifest changed during calibration"
            )
        if sha256_file(test_path) != test_sha:
            raise ValueError("Test prediction manifest changed during calibration")
        if sha256_file(obo_path) != obo_sha:
            raise ValueError("Calibration ontology changed during calibration")

        artifacts = output_manifest(
            stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"}
        )
        atomic_write_json(stage / "output_manifest.json", artifacts)
        atomic_write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "analysis_kind": "post_selection_calibration",
                "benchmark_id": validation["benchmark_id"],
                "mode": validation["mode"],
                "output_manifest_sha256": sha256_file(stage / "output_manifest.json"),
            },
        )
        os.replace(stage, output_dir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps(analysis, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
