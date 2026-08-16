#!/usr/bin/env python3
"""Analyse a completed contemporary NK+LK PFP prediction publication."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np


ASPECT_TO_ONTOLOGY = {"BPO": "bp", "CCO": "cc", "MFO": "mf"}
GLOBAL_NK_BENCHMARK_ID = "contemporary-2025-01-to-2026-02-supervisor"
NK_LK_BENCHMARK_ID = "contemporary-2025-01-to-2026-02-supervisor-nk-lk"


def _load_helpers(framework_root: Path) -> tuple[Any, Any, Any, Any]:
    diagnostics = framework_root / "scripts" / "diagnostics"
    sys.path.insert(0, str(diagnostics))
    from evaluate_pfp_information_content import (  # type: ignore[import-not-found]
        bootstrap_fixed_threshold,
        read_information_accretion,
        threshold_metrics,
    )
    from pfp_sensitivity_common import (  # type: ignore[import-not-found]
        load_aspect_bundle,
        verify_artifact_manifest,
    )

    return (
        bootstrap_fixed_threshold,
        read_information_accretion,
        threshold_metrics,
        (verify_artifact_manifest, load_aspect_bundle),
    )


def _load_membership(path: Path) -> dict[tuple[str, str], str]:
    result: dict[tuple[str, str], str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t", strict=True):
            key = (row["ontology"], row["protein"])
            if key in result:
                raise ValueError(f"Duplicate cohort membership: {key}")
            cohort = row["cohort"]
            if cohort not in {"no-knowledge", "limited-knowledge"}:
                raise ValueError(f"Unsupported cohort {cohort!r} for {key}")
            result[key] = cohort
    if not result:
        raise ValueError("Cohort membership table is empty")
    return result


def _interval_index(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {str(row["metric_name"]): row for row in rows}


def _fmt(value: Any, digits: int = 3) -> str:
    if value is None:
        return "NA"
    return f"{float(value):.{digits}f}"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_completed_run_file(path: Path) -> Path:
    for root in path.parents:
        marker_path = root / "WORKFLOW_COMPLETE.json"
        manifest_path = root / "output_manifest.json"
        if not marker_path.is_file() or not manifest_path.is_file():
            continue
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
        if marker.get("complete") is not True:
            raise ValueError(f"Model run is not complete: {root}")
        if marker.get("manifest_sha256") != _sha256(manifest_path):
            raise ValueError(f"Model run manifest hash is invalid: {root}")
        relative = path.resolve().relative_to(root.resolve()).as_posix()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        matches = [item for item in manifest.get("files", []) if item.get("path") == relative]
        if len(matches) != 1 or matches[0].get("sha256") != _sha256(path):
            raise ValueError(f"Global-NK summary is not authenticated by its model run: {path}")
        return root
    raise ValueError(f"Global-NK summary is not inside a completed model run: {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--framework-root", type=Path, required=True)
    parser.add_argument("--prediction-manifest", type=Path, required=True)
    parser.add_argument("--membership-tsv", type=Path, required=True)
    parser.add_argument("--nk-evaluation-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260805)
    args = parser.parse_args()

    if args.output_dir.exists():
        raise ValueError(f"Output directory already exists: {args.output_dir}")
    if args.bootstrap_replicates < 1:
        raise ValueError("--bootstrap-replicates must be positive")

    (
        bootstrap_fixed_threshold,
        read_information_accretion,
        threshold_metrics,
        artifact_helpers,
    ) = _load_helpers(args.framework_root.resolve())
    verify_artifact_manifest, load_aspect_bundle = artifact_helpers
    manifest, artifact_root = verify_artifact_manifest(
        args.prediction_manifest.resolve()
    )
    if manifest.get("benchmark_id") != NK_LK_BENCHMARK_ID:
        raise ValueError("Unexpected benchmark ID")
    if manifest.get("mode") != "full" or manifest.get("evaluation_split") != "test":
        raise ValueError("NK+LK analysis requires a completed full-model test capture")
    if manifest.get("selected_aspects") != list(ASPECT_TO_ONTOLOGY):
        raise ValueError("NK+LK analysis requires exactly BPO, CCO, and MFO")

    membership_path = args.membership_tsv.resolve()
    build_manifest_path = membership_path.parent / "build_manifest.json"
    if not build_manifest_path.is_file():
        raise ValueError("Cohort membership must be accompanied by build_manifest.json")
    build_manifest = json.loads(build_manifest_path.read_text(encoding="utf-8"))
    if build_manifest.get("profile") != "supervisor-nk-lk":
        raise ValueError("Cohort membership is not from the NK+LK builder profile")
    membership = _load_membership(membership_path)

    nk_summary_path = args.nk_evaluation_summary.resolve()
    nk_run_root = _verify_completed_run_file(nk_summary_path)
    nk_run_report_path = nk_run_root / "reports" / "run_report.json"
    _verify_completed_run_file(nk_run_report_path)
    nk_run_report = json.loads(nk_run_report_path.read_text(encoding="utf-8"))
    if (
        nk_run_report.get("benchmark_id") != GLOBAL_NK_BENCHMARK_ID
        or nk_run_report.get("execution_mode") != "train-eval"
        or nk_run_report.get("modality_mode") != "full"
    ):
        raise ValueError("Global-NK comparison is not the accepted full training run")
    nk_summary = json.loads(nk_summary_path.read_text(encoding="utf-8"))
    if nk_summary.get("mode") != "full" or nk_summary.get("evaluation_split") != "test":
        raise ValueError("Global-NK comparison must be a completed full-model test evaluation")
    if set(nk_summary.get("aspects", {})) != set(ASPECT_TO_ONTOLOGY):
        raise ValueError("Global-NK comparison must contain BPO, CCO, and MFO")

    args.output_dir.mkdir(parents=True)
    metric_rows: list[dict[str, Any]] = []
    bootstrap_rows: list[dict[str, Any]] = []
    composition_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    report_aspects: dict[str, Any] = {}

    for aspect_number, aspect in enumerate(manifest["selected_aspects"]):
        bundle = load_aspect_bundle(manifest, artifact_root, aspect)
        ontology = ASPECT_TO_ONTOLOGY[aspect]
        cohorts = []
        for protein_id in bundle["protein_ids"]:
            key = (ontology, protein_id)
            if key not in membership:
                raise ValueError(f"Missing cohort membership for {key}")
            cohorts.append(membership[key])
        expected_keys = {
            (ontology, protein_id) for protein_id in bundle["protein_ids"]
        }
        observed_keys = {key for key in membership if key[0] == ontology}
        if expected_keys != observed_keys:
            raise ValueError(
                f"Prediction/membership population differs for {aspect}: "
                f"missing={len(expected_keys-observed_keys)}, "
                f"extra={len(observed_keys-expected_keys)}"
            )

        keep = np.ones(len(bundle["go_terms"]), dtype=bool)
        keep[bundle["root_index"]] = False
        truth = bundle["truth"][:, keep]
        scores = bundle["scores"][:, keep]
        ia, ia_contract = read_information_accretion(
            bundle["ia_path"], bundle["go_terms"]
        )
        weights = ia[keep]
        threshold = float(bundle["specification"]["canonical_cafa_metrics"]["threshold"])
        has_nonroot = truth.any(axis=1)
        cohort_array = np.asarray(cohorts)

        masks = {
            "all": np.ones(len(cohorts), dtype=bool),
            "no-knowledge": cohort_array == "no-knowledge",
            "limited-knowledge": cohort_array == "limited-knowledge",
        }
        aspect_metrics: dict[str, Any] = {}
        for cohort_number, (cohort, cohort_mask) in enumerate(masks.items()):
            total = int(cohort_mask.sum())
            if total == 0:
                raise ValueError(f"Empty {cohort} cohort for {aspect}")
            informative = int(np.logical_and(cohort_mask, has_nonroot).sum())
            root_only = total - informative
            composition = {
                "aspect": aspect,
                "cohort": cohort,
                "total_proteins": total,
                "root_only_proteins": root_only,
                "root_only_percent": 100.0 * root_only / total,
                "nonroot_proteins": informative,
                "nonroot_percent": 100.0 * informative / total,
            }
            composition_rows.append(composition)

            aspect_metrics[cohort] = {}
            for population_number, (population, mask) in enumerate(
                (
                    ("all-test", cohort_mask),
                    ("nonroot-only", np.logical_and(cohort_mask, has_nonroot)),
                )
            ):
                selected_truth = truth[mask]
                selected_scores = scores[mask]
                if selected_truth.shape[0] == 0:
                    raise ValueError(
                        f"Empty {population} population for {aspect} {cohort}"
                    )
                metrics = threshold_metrics(
                    selected_truth, selected_scores, weights, threshold
                )
                row = {
                    "aspect": aspect,
                    "cohort": cohort,
                    "population": population,
                    "threshold": threshold,
                    **metrics,
                }
                metric_rows.append(row)
                intervals = bootstrap_fixed_threshold(
                    selected_truth,
                    selected_scores,
                    weights,
                    threshold,
                    args.bootstrap_replicates,
                    args.bootstrap_seed
                    + aspect_number * 100
                    + cohort_number * 10
                    + population_number,
                )
                for interval in intervals:
                    bootstrap_rows.append(
                        {
                            "aspect": aspect,
                            "cohort": cohort,
                            "population": population,
                            "threshold": threshold,
                            **interval,
                        }
                    )
                aspect_metrics[cohort][population] = {
                    "metrics": metrics,
                    "bootstrap": _interval_index(intervals),
                }

        current = bundle["specification"]["canonical_cafa_metrics"]
        previous = nk_summary["aspects"][aspect]
        comparison = {
            "aspect": aspect,
            "nk_fmax": float(previous["cafa_fmax"]),
            "nk_lk_fmax": float(current["fmax"]),
            "fmax_delta": float(current["fmax"] - previous["cafa_fmax"]),
            "nk_wfmax": float(previous["cafa_wfmax"]),
            "nk_lk_wfmax": float(current["wfmax"]),
            "wfmax_delta": float(current["wfmax"] - previous["cafa_wfmax"]),
            "nk_smin": float(previous["cafa_smin"]),
            "nk_lk_smin": float(current["smin"]),
            "smin_delta": float(current["smin"] - previous["cafa_smin"]),
            "nk_wsmin": float(previous["cafa_wsmin"]),
            "nk_lk_wsmin": float(current["wsmin"]),
            "wsmin_delta": float(current["wsmin"] - previous["cafa_wsmin"]),
        }
        comparison_rows.append(comparison)
        report_aspects[aspect] = {
            "threshold": threshold,
            "go_terms": len(bundle["go_terms"]),
            "ia_sha256": ia_contract["sha256"],
            "composition": {
                row["cohort"]: row
                for row in composition_rows
                if row["aspect"] == aspect
            },
            "cohort_metrics": aspect_metrics,
            "canonical_nk_vs_nk_lk": comparison,
        }

    def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
        fields = list(rows[0])
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    write_tsv(args.output_dir / "cohort_composition.tsv", composition_rows)
    write_tsv(args.output_dir / "cohort_fixed_threshold_metrics.tsv", metric_rows)
    write_tsv(args.output_dir / "cohort_bootstrap_intervals.tsv", bootstrap_rows)
    write_tsv(args.output_dir / "nk_vs_nk_lk_canonical_metrics.tsv", comparison_rows)

    report = {
        "schema_version": 1,
        "status": "complete",
        "analysis_kind": "contemporary_nk_lk_completed_model_analysis",
        "benchmark_id": manifest["benchmark_id"],
        "mode": manifest["mode"],
        "threshold_policy": (
            "one descriptive full-test CAFA-optimal threshold per aspect, held fixed "
            "for all NK/LK and all-test/nonroot-only subgroup comparisons"
        ),
        "bootstrap": {
            "replicates": args.bootstrap_replicates,
            "seed_base": args.bootstrap_seed,
            "unit": "protein",
            "interval": "percentile 95%",
        },
        "interpretation_limits": [
            "NK+LK versus NK-only canonical deltas compare separately retrained models and changed split populations.",
            "Cohort metrics are flat fixed-threshold diagnostics, not independent canonical CAFA evaluations.",
            "Limited knowledge means qualifying t0 knowledge in another ontology, not same-aspect partial knowledge.",
            "The thresholds are selected on the full test set and are descriptive, not deployable validation-fixed thresholds.",
        ],
        "inputs": {
            "prediction_manifest": {
                "path": str(args.prediction_manifest.resolve()),
                "sha256": _sha256(args.prediction_manifest.resolve()),
            },
            "membership_tsv": {
                "path": str(membership_path),
                "sha256": _sha256(membership_path),
            },
            "benchmark_build_manifest": {
                "path": str(build_manifest_path.resolve()),
                "sha256": _sha256(build_manifest_path.resolve()),
            },
            "nk_evaluation_summary": {
                "path": str(nk_summary_path),
                "sha256": _sha256(nk_summary_path),
            },
            "nk_run_report": {
                "path": str(nk_run_report_path),
                "sha256": _sha256(nk_run_report_path),
            },
        },
        "aspects": report_aspects,
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    lines = [
        "# Contemporary NK+LK Completed-Model Analysis",
        "",
        "## Canonical Full-Test Comparison",
        "",
        "| Aspect | NK Fmax | NK+LK Fmax | Delta | NK wFmax | NK+LK wFmax | Delta |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in comparison_rows:
        lines.append(
            f"| {row['aspect']} | {_fmt(row['nk_fmax'])} | {_fmt(row['nk_lk_fmax'])} | "
            f"{float(row['fmax_delta']):+.3f} | {_fmt(row['nk_wfmax'])} | "
            f"{_fmt(row['nk_lk_wfmax'])} | {float(row['wfmax_delta']):+.3f} |"
        )
    lines.extend(
        [
            "",
            "These are canonical CAFA metrics from separately retrained models. The deltas therefore combine the effect of the expanded test population with the changed disjoint training split.",
            "",
            "## Test Composition",
            "",
            "| Aspect | Cohort | Total | Root-only | Root-only % | Non-root | Non-root % |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in composition_rows:
        lines.append(
            f"| {row['aspect']} | {row['cohort']} | {row['total_proteins']:,} | "
            f"{row['root_only_proteins']:,} | {_fmt(row['root_only_percent'], 2)}% | "
            f"{row['nonroot_proteins']:,} | {_fmt(row['nonroot_percent'], 2)}% |"
        )
    lines.extend(
        [
            "",
            "## Cohort Performance at the Shared Aspect Threshold",
            "",
            "The threshold is the single full-test CAFA-optimal threshold for that aspect and is not re-optimised for either cohort.",
            "",
            "| Aspect | Population | Cohort | N | F | 95% CI | wF | 95% CI |",
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for aspect in manifest["selected_aspects"]:
        for population in ("all-test", "nonroot-only"):
            for cohort in ("no-knowledge", "limited-knowledge"):
                value = report_aspects[aspect]["cohort_metrics"][cohort][population]
                metrics = value["metrics"]
                intervals = value["bootstrap"]
                f_ci = intervals["f"]
                wf_ci = intervals["weighted_f"]
                lines.append(
                    f"| {aspect} | {population} | {cohort} | "
                    f"{metrics['population_proteins']:,} | {_fmt(metrics['f'])} | "
                    f"[{_fmt(f_ci['ci_low'])}, {_fmt(f_ci['ci_high'])}] | "
                    f"{_fmt(metrics['weighted_f'])} | "
                    f"[{_fmt(wf_ci['ci_low'])}, {_fmt(wf_ci['ci_high'])}] |"
                )
    lines.extend(
        [
            "",
            "## Interpretation Guardrails",
            "",
            "- Limited knowledge means prior qualifying knowledge in another ontology only; it is not same-aspect annotation extension.",
            "- Root-only and non-root-only results answer different questions and must remain visible separately.",
            "- Cohort results use a shared threshold for a fair within-aspect comparison; they are diagnostics rather than new canonical Fmax claims.",
            "- Bootstrap intervals quantify protein-sampling uncertainty only, not training-seed uncertainty.",
            "",
        ]
    )
    (args.output_dir / "NK+LK Contemporary Benchmark Analysis.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    shutil.copy2(Path(__file__).resolve(), args.output_dir / "analysis_driver.py")
    core_outputs = sorted(
        path
        for path in args.output_dir.iterdir()
        if path.is_file()
    )
    output_manifest = {
        "schema_version": 1,
        "status": "complete",
        "analysis_kind": "contemporary_nk_lk_completed_model_analysis",
        "files": [
            {
                "path": path.name,
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in core_outputs
        ],
        "analysis_script": {
            "path": str(Path(__file__).resolve()),
            "sha256": _sha256(Path(__file__).resolve()),
        },
    }
    manifest_path = args.output_dir / "output_manifest.json"
    manifest_path.write_text(
        json.dumps(output_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (args.output_dir / "ANALYSIS_COMPLETE.json").write_text(
        json.dumps(
            {
                "complete": True,
                "output_manifest_sha256": _sha256(manifest_path),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
