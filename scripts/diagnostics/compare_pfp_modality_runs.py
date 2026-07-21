#!/usr/bin/env python3
"""Compare true PFP modality retraining runs against sequence-only."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Mapping

from label_space_common import (
    ASPECTS,
    atomic_write_json,
    atomic_write_text,
    file_snapshot,
    output_manifest,
    require_unchanged,
    sha256_file,
)
from pfp_sensitivity_common import (
    MODE_MODALITIES,
    embedding_content_contract,
    verify_artifact_manifest,
)


SUPPORTED_MODES = (
    "sequence-only",
    "sequence-text",
    "sequence-structure",
    "sequence-ppi",
    "full",
)
BASELINE_MODE = "sequence-only"
STRICT_EVALUATOR_POLICY = "strict-ia-norm-cafa-prop-max-no-fallback"


def _manifest_index(manifest: Mapping[str, Any], run_root: Path) -> dict[str, Any]:
    listed: dict[str, Any] = {}
    for item in manifest.get("files", []):
        relative = Path(str(item.get("path", "")))
        if not relative.parts or relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"PFP output manifest contains unsafe path: {relative}")
        key = relative.as_posix()
        if key in listed:
            raise ValueError(f"PFP output manifest repeats {key}: {run_root}")
        listed[key] = item
    return listed


def _load_bound_json(
    run_root: Path,
    listed: Mapping[str, Any],
    relative: str,
    expected_sha256: str | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    path = run_root / relative
    snapshot = file_snapshot(path)
    specification = listed.get(relative)
    if not isinstance(specification, dict):
        raise ValueError(f"PFP output manifest does not bind {relative}: {run_root}")
    if (
        specification.get("bytes") != snapshot["bytes"]
        or specification.get("sha256") != snapshot["sha256"]
    ):
        raise ValueError(f"PFP output manifest binds different bytes: {path}")
    if expected_sha256 is not None and snapshot["sha256"] != expected_sha256:
        raise ValueError(f"PFP run report binds different {relative} bytes")
    value = json.loads(path.read_text(encoding="utf-8"))
    require_unchanged(path, snapshot, relative)
    return value, snapshot


def _embedding_contract(report: Mapping[str, Any], mode: str) -> dict[str, Any]:
    if report.get("status") != "passed" or report.get("mode") != mode:
        raise ValueError(f"Embedding report mode/status differs for {mode}")
    expected = set(MODE_MODALITIES[mode])
    observed = set(report.get("modalities", {}))
    if observed != expected:
        raise ValueError(
            f"Embedding report modalities differ for {mode}: expected "
            f"{sorted(expected)}, found {sorted(observed)}"
        )
    contents = {
        modality: report["modalities"][modality]["valid_content_sha256"]
        for modality in MODE_MODALITIES[mode]
    }
    if any(not isinstance(value, str) or len(value) != 64 for value in contents.values()):
        raise ValueError(f"Embedding report lacks valid content hashes for {mode}")
    binding = report.get("embedding_evidence_binding", {})
    return {
        "active_modalities": list(MODE_MODALITIES[mode]),
        "valid_content_sha256": contents,
        "information_accretion": {
            aspect: value.get("sha256")
            for aspect, value in report.get("information_accretion", {}).items()
        },
        "evidence_contract_sha256": binding.get("contract_sha256"),
        "target_manifest_sha256": binding.get("target_manifest_sha256"),
        "pair_status_sha256": binding.get("pair_status_sha256"),
    }


def load_run_report(path: Path) -> dict[str, Any]:
    path = path.resolve()
    run_root = path.parent.parent
    completion_path = run_root / "WORKFLOW_COMPLETE.json"
    manifest_path = run_root / "output_manifest.json"
    completion_snapshot = file_snapshot(completion_path)
    manifest_snapshot = file_snapshot(manifest_path)
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    if not completion.get("complete"):
        raise ValueError(f"PFP workflow is not complete: {run_root}")
    if completion.get("manifest_sha256") != manifest_snapshot["sha256"]:
        raise ValueError(f"PFP completion marker does not bind its manifest: {run_root}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    listed = _manifest_index(manifest, run_root)
    report, report_snapshot = _load_bound_json(
        run_root, listed, path.relative_to(run_root).as_posix()
    )
    if report.get("schema_version") != 1 or report.get("status") != "passed":
        raise ValueError(f"Unsupported or incomplete PFP run report: {path}")
    if report.get("execution_mode") != "train-eval":
        raise ValueError(
            f"Canonical modality comparison requires retrained train-eval runs: {path}"
        )
    mode = report.get("modality_mode")
    if mode not in SUPPORTED_MODES:
        raise ValueError(f"Unsupported modality mode in PFP report: {mode!r}")
    preparation, _ = _load_bound_json(
        run_root,
        listed,
        "reports/preparation.json",
        report.get("preparation_report_sha256"),
    )
    embedding, _ = _load_bound_json(
        run_root,
        listed,
        "reports/embedding_cache.json",
        report.get("embedding_report_sha256"),
    )
    embedding_post, _ = _load_bound_json(
        run_root,
        listed,
        "reports/embedding_cache_post.json",
        report.get("embedding_post_report_sha256"),
    )
    config_path = run_root / "run_config.json"
    config_snapshot = file_snapshot(config_path)
    config_specification = listed.get("run_config.json")
    if (
        not isinstance(config_specification, dict)
        or config_specification.get("bytes") != config_snapshot["bytes"]
        or config_specification.get("sha256") != config_snapshot["sha256"]
    ):
        raise ValueError(f"PFP output manifest does not bind run_config.json: {run_root}")
    before = _embedding_contract(embedding, mode)
    after = _embedding_contract(embedding_post, mode)
    if before != after:
        raise ValueError(f"Embedding scientific contract changed during {mode} run")
    aspects = report.get("aspects")
    if not isinstance(aspects, list) or not aspects or len(aspects) != len(set(aspects)):
        raise ValueError(f"Invalid PFP aspect inventory for {mode}")
    for aspect in aspects:
        result = report.get("results", {}).get(aspect, {})
        metrics = result.get("metrics", {})
        if metrics.get("cafa_evaluator_policy") != STRICT_EVALUATOR_POLICY:
            raise ValueError(f"Canonical evaluator policy differs for {mode}/{aspect}")
        for metric in ("cafa_fmax", "cafa_wfmax", "cafa_smin", "cafa_threshold"):
            if metric not in metrics or not math.isfinite(float(metrics[metric])):
                raise ValueError(f"Canonical metric {metric} is invalid for {mode}/{aspect}")
    require_unchanged(completion_path, completion_snapshot, "PFP completion marker")
    require_unchanged(manifest_path, manifest_snapshot, "PFP output manifest")
    return {
        "path": path,
        "run_root": run_root,
        "report": report,
        "report_snapshot": report_snapshot,
        "preparation": preparation,
        "embedding_contract": before,
        "config_sha256": config_snapshot["sha256"],
    }


def _parse_prediction_specs(values: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--prediction-manifest must use MODE=PATH")
        mode, raw_path = value.split("=", 1)
        if mode not in SUPPORTED_MODES:
            raise ValueError(f"Unsupported prediction modality mode: {mode!r}")
        if mode in result:
            raise ValueError(f"Repeated prediction modality mode: {mode}")
        result[mode] = Path(raw_path).resolve()
    return result


def _same_number(left: Any, right: Any) -> bool:
    return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=1e-12)


def bind_prediction_manifest(
    bundle: Mapping[str, Any],
    manifest_path: Path,
    allow_framework_commit_drift: bool,
) -> dict[str, Any]:
    manifest, artifact_root = verify_artifact_manifest(manifest_path)
    report = bundle["report"]
    mode = report["modality_mode"]
    expected_fields = {
        "mode": mode,
        "benchmark_id": report["benchmark_id"],
        "seed": report["seed"],
    }
    mismatched = [field for field, value in expected_fields.items() if manifest.get(field) != value]
    provenance = manifest["provenance"]
    if provenance.get("benchmark_fingerprint") != report["benchmark_fingerprint"]:
        mismatched.append("benchmark_fingerprint")
    if provenance.get("pfp_commit") != report["pfp_commit"]:
        mismatched.append("pfp_commit")
    framework_commit_match = (
        provenance.get("framework_commit") == report.get("framework_commit")
    )
    if not framework_commit_match and not allow_framework_commit_drift:
        mismatched.append("framework_commit")
    if provenance.get("source_csv_sha256", {}) != bundle["preparation"].get(
        "source_csv_sha256", {}
    ):
        mismatched.append("source_csv_sha256")
    if manifest["config"]["sha256"] != bundle["config_sha256"]:
        mismatched.append("config_sha256")
    if manifest.get("selected_aspects") != report.get("aspects"):
        mismatched.append("selected_aspects")
    prediction_embedding = embedding_content_contract(manifest, artifact_root)
    if prediction_embedding != bundle["embedding_contract"]:
        mismatched.append("embedding_content_contract")
    if mismatched:
        raise ValueError(
            f"Prediction capture differs from canonical {mode} run: "
            + ", ".join(sorted(set(mismatched)))
        )
    for aspect in report["aspects"]:
        specification = manifest["aspects"][aspect]
        canonical = report["results"][aspect]
        if specification["checkpoint_sha256"] != canonical["checkpoint_sha256"]:
            raise ValueError(f"Prediction checkpoint differs from canonical {mode}/{aspect}")
        expected_ia = bundle["embedding_contract"]["information_accretion"].get(aspect)
        if expected_ia is not None and specification.get("ia_file_sha256") != expected_ia:
            raise ValueError(f"Prediction IA differs from canonical {mode}/{aspect}")
        captured = specification["canonical_cafa_metrics"]
        for short_name, report_name in (
            ("fmax", "cafa_fmax"),
            ("wfmax", "cafa_wfmax"),
            ("smin", "cafa_smin"),
            ("threshold", "cafa_threshold"),
        ):
            if short_name not in captured or not _same_number(
                captured[short_name], canonical["metrics"][report_name]
            ):
                raise ValueError(
                    f"Prediction canonical {short_name} differs for {mode}/{aspect}"
                )
    return {
        "path": str(manifest_path),
        "sha256": sha256_file(manifest_path),
        "run_framework_commit": provenance.get("framework_commit"),
        "canonical_framework_commit": report.get("framework_commit"),
        "framework_commit_match": framework_commit_match,
        "ia_binding": {
            aspect: (
                "exact_precomputed_sha256"
                if bundle["embedding_contract"]["information_accretion"].get(aspect)
                is not None
                else "computed_ia_bound_by_checkpoint_and_exact_canonical_metrics"
            )
            for aspect in report["aspects"]
        },
        "obo_sha256": manifest["obo"]["sha256"],
    }


def tsv_text(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return ""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PFP Modality-Contribution Comparison",
        "",
        "Canonical metrics come only from fresh retraining runs. Prediction artifacts are accepted only after exact checkpoint and scientific-contract binding.",
        "",
        "| Mode | Aspect | Fmax | wFmax | Smin |",
        "|---|---|---:|---:|---:|",
    ]
    for row in report["run_rows"]:
        lines.append(
            f"| {row['mode']} | {row['aspect']} | {row['cafa_fmax']:.6f} | "
            f"{row['cafa_wfmax']:.6f} | {row['cafa_smin']:.6f} |"
        )
    lines.extend(
        [
            "",
            "## Improvement Over Sequence-Only",
            "",
            "| Mode | Aspect | Fmax improvement | wFmax improvement | Smin improvement |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in report["delta_rows"]:
        lines.append(
            f"| {row['comparison_mode']} | {row['aspect']} | "
            f"{row['cafa_fmax_improvement']:+.6f} | "
            f"{row['cafa_wfmax_improvement']:+.6f} | "
            f"{row['cafa_smin_improvement']:+.6f} |"
        )
    lines.extend(
        [
            "",
            "Fmax is primary. wFmax and Smin remain benchmark-specific secondary metrics.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-report", action="append", type=Path, required=True)
    parser.add_argument("--prediction-manifest", action="append", default=[])
    parser.add_argument("--allow-framework-commit-drift", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if len(args.run_report) < 2:
        raise ValueError("At least two completed PFP run reports are required")
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")

    bundles = [load_run_report(path) for path in args.run_report]
    reports = [bundle["report"] for bundle in bundles]
    modes = [report["modality_mode"] for report in reports]
    if len(modes) != len(set(modes)):
        raise ValueError("PFP run reports repeat a modality mode")
    if BASELINE_MODE not in modes:
        raise ValueError("A sequence-only run is required as the comparison baseline")
    by_mode = {bundle["report"]["modality_mode"]: bundle for bundle in bundles}

    reference = bundles[0]
    reference_report = reference["report"]
    contract_fields = ("benchmark_id", "benchmark_fingerprint", "seed", "pfp_commit")
    for bundle in bundles[1:]:
        report = bundle["report"]
        mismatched = [
            field for field in contract_fields
            if report.get(field) != reference_report.get(field)
        ]
        if report.get("aspects") != reference_report.get("aspects"):
            mismatched.append("aspects")
        if report.get("environment") != reference_report.get("environment"):
            mismatched.append("environment")
        if bundle["config_sha256"] != reference["config_sha256"]:
            mismatched.append("config_sha256")
        if bundle["preparation"].get("source_csv_sha256", {}) != reference[
            "preparation"
        ].get("source_csv_sha256", {}):
            mismatched.append("source_csv_sha256")
        for shared in (
            "evidence_contract_sha256",
            "target_manifest_sha256",
            "pair_status_sha256",
        ):
            if bundle["embedding_contract"].get(shared) != reference[
                "embedding_contract"
            ].get(shared):
                mismatched.append(shared)
        if bundle["embedding_contract"]["valid_content_sha256"]["sequence"] != reference[
            "embedding_contract"
        ]["valid_content_sha256"]["sequence"]:
            mismatched.append("sequence_embedding_content")
        if bundle["embedding_contract"]["information_accretion"] != reference[
            "embedding_contract"
        ]["information_accretion"]:
            mismatched.append("information_accretion")
        if mismatched:
            raise ValueError(
                f"PFP run comparison contract differs for {report['modality_mode']}: "
                + ", ".join(sorted(set(mismatched)))
            )

    full = by_mode.get("full")
    if full is not None:
        full_hashes = full["embedding_contract"]["valid_content_sha256"]
        for mode, bundle in by_mode.items():
            for modality in MODE_MODALITIES[mode]:
                if bundle["embedding_contract"]["valid_content_sha256"][modality] != full_hashes[modality]:
                    raise ValueError(
                        f"Active {modality} embedding content differs for {mode} and full"
                    )

    framework_commits = {report["framework_commit"] for report in reports}
    if len(framework_commits) > 1 and not args.allow_framework_commit_drift:
        raise ValueError(
            "Framework commits differ; audit the drift and pass "
            "--allow-framework-commit-drift explicitly"
        )

    prediction_specs = _parse_prediction_specs(args.prediction_manifest)
    prediction_sources: list[dict[str, Any]] = []
    if prediction_specs:
        if set(prediction_specs) != set(modes):
            raise ValueError(
                "Prediction manifest modes must exactly match canonical run modes"
            )
        for mode in SUPPORTED_MODES:
            if mode not in prediction_specs:
                continue
            source = bind_prediction_manifest(
                by_mode[mode],
                prediction_specs[mode],
                args.allow_framework_commit_drift,
            )
            source["mode"] = mode
            prediction_sources.append(source)
        if len({source["obo_sha256"] for source in prediction_sources}) != 1:
            raise ValueError("Prediction manifests use different OBO files")

    run_rows: list[dict[str, Any]] = []
    indexed: dict[tuple[str, str], dict[str, Any]] = {}
    for report in reports:
        mode = report["modality_mode"]
        for aspect in report["aspects"]:
            metrics = report["results"][aspect]["metrics"]
            row = {
                "benchmark_id": report["benchmark_id"],
                "mode": mode,
                "aspect": aspect,
                "cafa_fmax": float(metrics["cafa_fmax"]),
                "cafa_wfmax": float(metrics["cafa_wfmax"]),
                "cafa_smin": float(metrics["cafa_smin"]),
                "checkpoint_sha256": report["results"][aspect]["checkpoint_sha256"],
                "framework_commit": report["framework_commit"],
                "pfp_commit": report["pfp_commit"],
            }
            run_rows.append(row)
            indexed[(mode, aspect)] = row

    delta_rows: list[dict[str, Any]] = []
    for mode in SUPPORTED_MODES:
        if mode == BASELINE_MODE or mode not in modes:
            continue
        for aspect in reference_report["aspects"]:
            candidate = indexed[(mode, aspect)]
            baseline = indexed[(BASELINE_MODE, aspect)]
            delta_rows.append(
                {
                    "benchmark_id": reference_report["benchmark_id"],
                    "comparison_mode": mode,
                    "baseline_mode": BASELINE_MODE,
                    "aspect": aspect,
                    "cafa_fmax_improvement": candidate["cafa_fmax"] - baseline["cafa_fmax"],
                    "cafa_wfmax_improvement": candidate["cafa_wfmax"] - baseline["cafa_wfmax"],
                    "cafa_smin_improvement": baseline["cafa_smin"] - candidate["cafa_smin"],
                    "framework_commit_match": candidate["framework_commit"] == baseline["framework_commit"],
                }
            )

    comparison = {
        "schema_version": 2,
        "status": "complete",
        "benchmark_id": reference_report["benchmark_id"],
        "baseline_mode": BASELINE_MODE,
        "comparison_contract": {
            **{field: reference_report.get(field) for field in contract_fields},
            "config_sha256": reference["config_sha256"],
            "source_csv_sha256": reference["preparation"].get("source_csv_sha256", {}),
            "sequence_embedding_content_sha256": reference["embedding_contract"]["valid_content_sha256"]["sequence"],
            "information_accretion": reference["embedding_contract"]["information_accretion"],
            "framework_commit_drift_allowed": args.allow_framework_commit_drift,
        },
        "source_reports": [
            {
                "mode": bundle["report"]["modality_mode"],
                "path": str(bundle["path"]),
                "sha256": bundle["report_snapshot"]["sha256"],
            }
            for bundle in bundles
        ],
        "prediction_sources": prediction_sources,
        "run_rows": sorted(
            run_rows,
            key=lambda row: (SUPPORTED_MODES.index(row["mode"]), ASPECTS.index(row["aspect"])),
        ),
        "delta_rows": delta_rows,
        "metric_policy": {
            "primary": "cafa_fmax",
            "secondary": ["cafa_wfmax", "cafa_smin"],
            "positive_improvement_is_better": True,
        },
    }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.stage-", dir=str(output_dir.parent)))
    try:
        atomic_write_json(stage / "modality_comparison.json", comparison)
        atomic_write_text(stage / "modality_comparison.md", markdown_report(comparison))
        atomic_write_text(stage / "canonical_metrics.tsv", tsv_text(comparison["run_rows"]))
        atomic_write_text(stage / "mode_vs_sequence_deltas.tsv", tsv_text(delta_rows))
        manifest = output_manifest(stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"})
        atomic_write_json(stage / "output_manifest.json", manifest)
        atomic_write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "output_manifest_sha256": sha256_file(stage / "output_manifest.json"),
            },
        )
        os.replace(stage, output_dir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps(comparison, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
