#!/usr/bin/env python3
"""Build the compact provenance and results report for a generic PFP run."""

from __future__ import annotations

import argparse
import json
import math
import platform
from pathlib import Path
from typing import Any, Dict, Optional

from importlib.metadata import PackageNotFoundError, version

from common import (
    atomic_write_json,
    atomic_write_text,
    expected_result_dir,
    load_json,
    selected_aspects,
    sha256_file,
    validate_mandatory_metrics,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument("--execution-mode", choices=("prepare-only", "eval-only", "train-eval"), required=True)
    parser.add_argument("--modality-mode", choices=("full", "sequence-only"), required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--framework-commit", default="unknown")
    parser.add_argument("--pfp-commit", default="unknown")
    parser.add_argument("--preparation-report", type=Path, required=True)
    parser.add_argument("--embedding-report", type=Path)
    parser.add_argument("--embedding-post-report", type=Path)
    parser.add_argument("--result-root", type=Path)
    parser.add_argument("--evaluation-root", type=Path)
    parser.add_argument("--expected-metrics", type=Path)
    parser.add_argument("--reference-tolerance", type=float)
    parser.add_argument("--require-reference-match", action="store_true")
    parser.add_argument("--aspect", action="append", default=[])
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()

    aspects = selected_aspects(args.aspect)
    preparation = load_json(args.preparation_report)
    if preparation.get("status") != "passed":
        raise ValueError("Preparation report is not a passed report")
    embedding = load_json(args.embedding_report) if args.embedding_report else None
    if embedding and embedding.get("status") != "passed":
        raise ValueError("Embedding report is not a passed report")
    embedding_post = (
        load_json(args.embedding_post_report) if args.embedding_post_report else None
    )
    if embedding_post and embedding_post.get("status") != "passed":
        raise ValueError("Post-execution embedding report is not a passed report")

    results: Dict[str, Any] = {}
    if args.execution_mode == "train-eval":
        if not args.result_root or not args.evaluation_root:
            raise ValueError("--result-root and --evaluation-root are required for train-eval")
        for aspect in aspects:
            directory = expected_result_dir(args.result_root, aspect)
            checkpoint = directory / "best_model.pt"
            training_result_file = directory / "results.json"
            strict_result_file = args.evaluation_root / aspect / "results.json"
            if not checkpoint.is_file() or not training_result_file.is_file():
                raise FileNotFoundError(f"Incomplete PFP result for {aspect}: {directory}")
            if not strict_result_file.is_file():
                raise FileNotFoundError(
                    f"Strict post-training evaluation is missing for {aspect}: {strict_result_file}"
                )
            training_result = load_json(training_result_file)
            result = load_json(strict_result_file)
            validate_mandatory_metrics(result, f"Strict evaluation result for {aspect}")
            results[aspect] = {
                "metrics": result,
                "checkpoint_sha256": sha256_file(checkpoint),
                "results_sha256": sha256_file(strict_result_file),
                "training_results_sha256": sha256_file(training_result_file),
                "training_metrics": training_result,
            }
    elif args.execution_mode == "eval-only":
        if not args.evaluation_root:
            raise ValueError("--evaluation-root is required for eval-only")
        for aspect in aspects:
            result_file = args.evaluation_root / aspect / "results.json"
            result = load_json(result_file)
            validate_mandatory_metrics(result, f"Evaluation result for {aspect}")
            results[aspect] = {
                "metrics": result,
                "checkpoint_sha256": result["checkpoint_sha256"],
                "results_sha256": sha256_file(result_file),
            }

    reference: Optional[Dict[str, Any]] = None
    reference_passed = True
    if args.expected_metrics:
        expected = load_json(args.expected_metrics)
        tolerance = args.reference_tolerance
        if tolerance is None:
            tolerance = float(expected.get("comparison_policy", {}).get("eval_only_tolerance", 0.0005))
        comparisons: list[Dict[str, Any]] = []
        expected_by_aspect = expected.get("metrics", {})
        missing_reference_aspects = [aspect for aspect in aspects if aspect not in expected_by_aspect]
        if missing_reference_aspects:
            reference_passed = False
            comparisons.extend(
                {
                    "aspect": aspect,
                    "metric": "*",
                    "status": "missing_reference_aspect",
                }
                for aspect in missing_reference_aspects
            )
        if args.execution_mode == "prepare-only":
            reference_passed = False
            comparisons.append(
                {
                    "aspect": "*",
                    "metric": "*",
                    "status": "no_results_in_prepare_only_mode",
                }
            )
        for aspect, expected_metrics in expected.get("metrics", {}).items():
            if aspect not in results:
                continue
            observed = results[aspect]["metrics"]
            for metric, wanted in expected_metrics.items():
                if not math.isfinite(float(wanted)):
                    raise ValueError(f"Expected metric is non-finite: {aspect}/{metric}")
                if metric not in observed:
                    comparisons.append(
                        {"aspect": aspect, "metric": metric, "status": "missing", "expected": wanted}
                    )
                    reference_passed = False
                    continue
                observed_value = float(observed[metric])
                if not math.isfinite(observed_value):
                    raise ValueError(f"Observed metric is non-finite: {aspect}/{metric}")
                difference = abs(observed_value - float(wanted))
                passed = difference <= tolerance
                reference_passed = reference_passed and passed
                comparisons.append(
                    {
                        "aspect": aspect,
                        "metric": metric,
                        "expected": float(wanted),
                        "observed": float(observed[metric]),
                        "absolute_difference": difference,
                        "tolerance": tolerance,
                        "status": "passed" if passed else "outside_tolerance",
                    }
                )
        if not comparisons:
            reference_passed = False
            comparisons.append(
                {"aspect": "*", "metric": "*", "status": "no_reference_comparisons"}
            )
        reference = {
            "source": str(args.expected_metrics.resolve()),
            "source_sha256": sha256_file(args.expected_metrics),
            "passed": reference_passed,
            "required": args.require_reference_match,
            "comparisons": comparisons,
        }
    reference_gate_failed = bool(
        args.require_reference_match and (reference is None or not reference_passed)
    )
    distributions = (
        "torch",
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "cafaeval",
        "obonet",
        "networkx",
    )
    packages = {}
    for distribution in distributions:
        try:
            packages[distribution] = version(distribution)
        except PackageNotFoundError:
            packages[distribution] = "missing"

    report = {
        "schema_version": 1,
        "status": "failed_reference_match" if reference_gate_failed else "passed",
        "benchmark_id": args.benchmark_id,
        "benchmark_fingerprint": preparation["benchmark_fingerprint"],
        "execution_mode": args.execution_mode,
        "modality_mode": args.modality_mode,
        "seed": args.seed,
        "aspects": aspects,
        "framework_commit": args.framework_commit,
        "pfp_commit": args.pfp_commit,
        "environment": {"python": platform.python_version(), "packages": packages},
        "preparation_report_sha256": sha256_file(args.preparation_report),
        "embedding_report_sha256": sha256_file(args.embedding_report) if args.embedding_report else None,
        "embedding_post_report_sha256": (
            sha256_file(args.embedding_post_report) if args.embedding_post_report else None
        ),
        "results": results,
        "reference_comparison": reference,
        "metric_interpretation": {
            "primary_cross_benchmark": "cafa_fmax",
            "secondary_benchmark_specific": ["cafa_wfmax", "cafa_smin"],
            "reason": "IA-weighted metrics depend on each benchmark's training labels and ontology snapshot.",
        },
    }
    atomic_write_json(args.output_json, report)

    lines = [
        f"# PFP benchmark run: {args.benchmark_id}",
        "",
        f"- Status: **{'FAILED REFERENCE MATCH' if reference_gate_failed else 'PASS'}**",
        f"- Execution: `{args.execution_mode}`",
        f"- Modalities: `{args.modality_mode}`",
        f"- Seed: `{args.seed}`",
        f"- Benchmark fingerprint: `{preparation['benchmark_fingerprint']}`",
        "",
    ]
    if embedding:
        lines.extend(["## Embedding coverage", ""])
        for modality, value in embedding["modalities"].items():
            lines.append(
                f"- {modality}: {value['valid']}/{value['target_count']} "
                f"({100.0 * value['coverage']:.2f}%) valid"
            )
        lines.append("")
    if results:
        lines.extend(["## Results", "", "| Aspect | CAFA Fmax | Weighted Fmax | Smin |", "|---|---:|---:|---:|"])
        for aspect in aspects:
            metrics = results[aspect]["metrics"]
            lines.append(
                f"| {aspect} | {metrics['cafa_fmax']:.6f} | "
                f"{metrics['cafa_wfmax']:.6f} | {metrics['cafa_smin']:.6f} |"
            )
        lines.extend(
            [
                "",
                "Ordinary CAFA Fmax is the primary cross-benchmark comparison. Weighted Fmax and Smin are",
                "reported as benchmark-specific secondary metrics because ontology and information-accretion",
                "inputs can differ between benchmarks.",
                "",
            ]
        )
    if reference:
        lines.extend(
            [
                "## Reference comparison",
                "",
                f"- Result: **{'PASS' if reference_passed else 'OUTSIDE TOLERANCE'}**",
                f"- Required gate: `{args.require_reference_match}`",
                "",
            ]
        )
    atomic_write_text(args.output_md, "\n".join(lines).rstrip() + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 1 if reference_gate_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
