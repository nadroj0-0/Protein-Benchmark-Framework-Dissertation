#!/usr/bin/env python3
"""Run opt-in root-only sensitivity analyses on captured PFP predictions."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any

from label_space_common import (
    ASPECTS,
    ASPECT_TO_ROOT,
    atomic_write_json,
    atomic_write_text,
    output_manifest,
    peak_rss_bytes,
    read_obo,
    sha256_file,
)
from pfp_sensitivity_common import (
    aspect_comparison_contract,
    cafaeval_contract,
    cohort_masks,
    flat_non_root_metrics,
    global_comparison_contract,
    load_aspect_bundle,
    run_cafa_evaluation,
    verify_artifact_manifest,
)


def selected_aspects(values: list[str], available: list[str]) -> list[str]:
    result = values or available
    unknown = sorted(set(result) - set(ASPECTS))
    unavailable = sorted(set(result) - set(available))
    if unknown:
        raise ValueError(f"Unknown PFP aspects: {unknown}")
    if unavailable:
        raise ValueError(f"Prediction artifact does not contain aspects: {unavailable}")
    if len(result) != len(set(result)):
        raise ValueError("Each PFP aspect may be selected only once")
    return result


def tsv_text(rows: list[dict[str, Any]]) -> str:
    fields = ["benchmark_id", "mode", "aspect", "protein_id", "cohort"]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def assert_canonical_reproduction(
    captured: dict[str, float],
    rerun: dict[str, Any],
    tolerance: float,
    aspect: str,
) -> dict[str, float]:
    deltas = {}
    for captured_key, rerun_key in (
        ("fmax", "fmax"),
        ("wfmax", "wfmax"),
        ("smin", "smin"),
        ("threshold", "threshold"),
    ):
        if captured_key not in captured:
            raise ValueError(
                f"Captured canonical metrics lack {captured_key} for {aspect}"
            )
        delta = abs(float(captured[captured_key]) - float(rerun[rerun_key]))
        deltas[captured_key] = delta
        if not math.isfinite(delta) or delta > tolerance:
            raise ValueError(
                f"Canonical {captured_key} did not reproduce for {aspect}: "
                f"delta={delta}, tolerance={tolerance}"
            )
    return deltas


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        f"# PFP Root-Only Sensitivity: {report['benchmark_id']} ({report['mode']})",
        "",
        "This is a separate diagnostic analysis. It does not replace or modify the canonical evaluation.",
        "",
        "| Aspect | Captured rows | CAFA-evaluable | All-zero | Root-only | Eligible after exclusion | Canonical Fmax | Excluded Fmax | At original threshold | Root baseline Fmax |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for aspect in report["selected_aspects"]:
        value = report["aspects"][aspect]
        excluded = value["root_only_excluded"]
        if excluded["status"] == "complete":
            excluded_f = f"{excluded['cafa']['fmax']:.4f}"
            fixed_f = f"{excluded['cafa']['fixed_at_canonical_threshold']['f']:.4f}"
        else:
            excluded_f = "n/a"
            fixed_f = "n/a"
        lines.append(
            "| {aspect} | {captured:,} | {evaluable:,} | {all_zero:,} | {root_only:,} | {eligible:,} | {canonical:.4f} | "
            "{excluded} | {fixed} | {baseline:.4f} |".format(
                aspect=aspect,
                captured=value["cohorts"]["captured_rows"],
                evaluable=value["cohorts"]["cafaeval_evaluable_targets"],
                all_zero=value["cohorts"]["all_zero"],
                root_only=value["cohorts"]["root_only"],
                eligible=value["cohorts"]["eligible_non_root"],
                canonical=value["canonical_recheck"]["fmax"],
                excluded=excluded_f,
                fixed=fixed_f,
                baseline=value["root_only_prediction_baseline"]["fmax"],
            )
        )
    lines.extend(
        [
            "",
            "- `canonical_recheck` reruns the captured predictions and must reproduce the original strict cafaeval metrics.",
            "- `root_only_excluded` removes test proteins with no positive non-root term; it retrains nothing.",
            "- All-zero rows remain reported as benchmark rows but have no truth record and are not cafaeval targets.",
            "- The original-threshold result uses this run's own canonical Fmax threshold, while the excluded Fmax re-optimizes the threshold and is labelled accordingly.",
            "- `flat_non_root_diagnostic` removes the root column without ontology propagation and is not a CAFA metric.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction-manifest", type=Path, required=True)
    parser.add_argument("--obo-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--aspect", action="append", default=[])
    parser.add_argument("--metric-tolerance", type=float, default=1e-8)
    args = parser.parse_args()

    manifest_path = args.prediction_manifest.resolve()
    obo_file = args.obo_file.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")
    if not obo_file.is_file():
        raise FileNotFoundError(f"GO OBO file is missing: {obo_file}")
    if not math.isfinite(args.metric_tolerance) or args.metric_tolerance < 0:
        raise ValueError("--metric-tolerance must be finite and non-negative")

    prediction_manifest_sha256 = sha256_file(manifest_path)
    manifest, artifact_root = verify_artifact_manifest(manifest_path)
    if manifest["obo"]["sha256"] != sha256_file(obo_file):
        raise ValueError("Sensitivity OBO differs from the prediction artifact OBO")
    read_obo(obo_file)
    aspects = selected_aspects(args.aspect, list(manifest["selected_aspects"]))

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.stage-", dir=str(output_dir.parent))
    )
    started = time.perf_counter()
    comparison_contract = global_comparison_contract(manifest)
    report = {
        "schema_version": 2,
        "status": "complete",
        "analysis_policy": (
            "canonical results retained; root-only exclusion is an opt-in sensitivity "
            "analysis with no retraining"
        ),
        "benchmark_id": manifest["benchmark_id"],
        "mode": manifest["mode"],
        "comparison_contract": comparison_contract,
        "code_provenance": {
            "framework_commit": manifest["provenance"]["framework_commit"],
            "pfp_commit": manifest["provenance"]["pfp_commit"],
        },
        "cafaeval_contract": cafaeval_contract(),
        "selected_aspects": aspects,
        "prediction_manifest": {
            "path": str(manifest_path),
            "sha256": prediction_manifest_sha256,
        },
        "obo": {"path": str(obo_file), "sha256": sha256_file(obo_file)},
        "metric_tolerance": args.metric_tolerance,
        "aspects": {},
    }
    cohort_rows: list[dict[str, Any]] = []
    try:
        for aspect in aspects:
            bundle = load_aspect_bundle(manifest, artifact_root, aspect)
            scores = bundle["scores"]
            truth = bundle["truth"]
            protein_ids = bundle["protein_ids"]
            go_terms = bundle["go_terms"]
            root_index = bundle["root_index"]
            masks = cohort_masks(truth, root_index)
            evaluable = ~masks["all_zero"]
            captured = bundle["specification"]["canonical_cafa_metrics"]
            canonical_threshold = float(captured["threshold"])

            canonical = run_cafa_evaluation(
                obo_file=obo_file,
                ia_file=bundle["ia_path"],
                destination=stage / aspect / "canonical_all_targets",
                protein_ids=protein_ids,
                go_terms=go_terms,
                truth=truth,
                fixed_threshold=canonical_threshold,
                scores=scores,
            )
            deltas = assert_canonical_reproduction(
                captured, canonical, args.metric_tolerance, aspect
            )

            eligible = masks["eligible_non_root"]
            if eligible.any():
                eligible_ids = [
                    protein_id
                    for protein_id, keep in zip(protein_ids, eligible.tolist())
                    if keep
                ]
                excluded_cafa = run_cafa_evaluation(
                    obo_file=obo_file,
                    ia_file=bundle["ia_path"],
                    destination=stage / aspect / "root_only_excluded",
                    protein_ids=eligible_ids,
                    go_terms=go_terms,
                    truth=truth[eligible],
                    fixed_threshold=canonical_threshold,
                    scores=scores[eligible],
                )
                excluded = {"status": "complete", "cafa": excluded_cafa}
            else:
                excluded = {"status": "not_evaluable_no_non_root_targets"}

            root_baseline = run_cafa_evaluation(
                obo_file=obo_file,
                ia_file=bundle["ia_path"],
                destination=stage / aspect / "root_only_prediction_baseline",
                protein_ids=[
                    protein_id
                    for protein_id, keep in zip(protein_ids, evaluable.tolist())
                    if keep
                ],
                go_terms=go_terms,
                truth=truth[evaluable],
                fixed_threshold=canonical_threshold,
                scores=None,
                root_index=root_index,
            )
            flat = flat_non_root_metrics(
                truth, scores, root_index, canonical_threshold
            )
            if (
                sha256_file(bundle["ia_path"])
                != bundle["specification"]["ia_file_sha256"]
            ):
                raise ValueError(f"IA file changed during sensitivity analysis for {aspect}")

            for row_index, protein_id in enumerate(protein_ids):
                if masks["root_only"][row_index]:
                    cohort = "root-only"
                elif masks["eligible_non_root"][row_index]:
                    cohort = "eligible-non-root"
                else:
                    cohort = "all-zero"
                cohort_rows.append(
                    {
                        "benchmark_id": manifest["benchmark_id"],
                        "mode": manifest["mode"],
                        "aspect": aspect,
                        "protein_id": protein_id,
                        "cohort": cohort,
                    }
                )

            report["aspects"][aspect] = {
                "root": ASPECT_TO_ROOT[aspect],
                "checkpoint_sha256": bundle["specification"]["checkpoint_sha256"],
                "ia_file_sha256": bundle["specification"]["ia_file_sha256"],
                "comparison_contract": aspect_comparison_contract(
                    bundle["specification"]
                ),
                "cohorts": {
                    "captured_rows": len(protein_ids),
                    "cafaeval_evaluable_targets": int(evaluable.sum()),
                    "eligible_non_root": int(masks["eligible_non_root"].sum()),
                    "eligible_non_root_fraction": float(masks["eligible_non_root"].mean()),
                    "root_only": int(masks["root_only"].sum()),
                    "root_only_fraction_of_captured": float(masks["root_only"].mean()),
                    "root_only_fraction_of_cafaeval_evaluable": (
                        float(masks["root_only"].sum() / evaluable.sum())
                        if evaluable.any()
                        else None
                    ),
                    "all_zero": int(masks["all_zero"].sum()),
                },
                "captured_canonical_metrics": captured,
                "canonical_recheck": canonical,
                "canonical_recheck_absolute_deltas": deltas,
                "root_only_excluded": excluded,
                "root_only_prediction_baseline": root_baseline,
                "flat_non_root_diagnostic": flat,
            }

        if sha256_file(obo_file) != manifest["obo"]["sha256"]:
            raise ValueError("Sensitivity OBO changed during analysis")
        if sha256_file(manifest_path) != prediction_manifest_sha256:
            raise ValueError("Prediction manifest changed during sensitivity analysis")
        report["resource_usage"] = {
            "wall_seconds": time.perf_counter() - started,
            "peak_rss_bytes": peak_rss_bytes(),
        }
        atomic_write_json(stage / "root_exclusion_sensitivity.json", report)
        atomic_write_text(stage / "root_exclusion_sensitivity.md", markdown_report(report))
        atomic_write_text(stage / "target_cohorts.tsv", tsv_text(cohort_rows))
        artifact_manifest = output_manifest(
            stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"}
        )
        atomic_write_json(stage / "output_manifest.json", artifact_manifest)
        atomic_write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
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
