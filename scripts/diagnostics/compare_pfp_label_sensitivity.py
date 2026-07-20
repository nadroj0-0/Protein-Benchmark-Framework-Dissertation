#!/usr/bin/env python3
"""Compare root-only sensitivity reports across modes and benchmarks."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

from label_space_common import (
    ASPECTS,
    atomic_write_json,
    atomic_write_text,
    file_snapshot,
    output_manifest,
    require_unchanged,
    sha256_file,
)


def load_report(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    completion_path = path.parent / "RUN_COMPLETE.json"
    output_path = path.parent / "output_manifest.json"
    report_snapshot = file_snapshot(path)
    completion_snapshot = file_snapshot(completion_path)
    output_snapshot = file_snapshot(output_path)
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    if not completion.get("complete"):
        raise ValueError(f"Sensitivity completion marker is false: {path.parent}")
    if completion.get("output_manifest_sha256") != sha256_file(output_path):
        raise ValueError(f"Sensitivity completion marker does not bind its manifest: {path}")
    output = json.loads(output_path.read_text(encoding="utf-8"))
    listed = {item.get("path"): item for item in output.get("files", [])}
    if path.name not in listed:
        raise ValueError(f"Sensitivity output manifest does not bind its report: {path}")
    if (
        listed[path.name].get("bytes") != report_snapshot["bytes"]
        or listed[path.name].get("sha256") != report_snapshot["sha256"]
    ):
        raise ValueError(f"Sensitivity manifest binds different report bytes: {path}")
    for item in output.get("files", []):
        relative = Path(item["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Sensitivity manifest contains unsafe path: {relative}")
        artifact = path.parent / relative
        if not artifact.is_file():
            raise FileNotFoundError(f"Sensitivity artifact is missing: {artifact}")
        if (
            artifact.stat().st_size != item["bytes"]
            or sha256_file(artifact) != item["sha256"]
        ):
            raise ValueError(f"Sensitivity artifact changed after publication: {artifact}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("schema_version") != 2 or value.get("status") != "complete":
        raise ValueError(f"Unsupported or incomplete sensitivity report: {path}")
    require_unchanged(path, report_snapshot, "Sensitivity report")
    require_unchanged(completion_path, completion_snapshot, "Sensitivity completion marker")
    require_unchanged(output_path, output_snapshot, "Sensitivity output manifest")
    return value, report_snapshot


def tsv_text(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return ""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PFP Label-Sensitivity Comparison",
        "",
        "Canonical results remain primary. All exclusions and mode deltas below are sensitivity analyses.",
        "",
        "## Run Profiles",
        "",
        "| Benchmark | Mode | Aspect | Captured | CAFA-evaluable | All-zero | Root-only % of evaluable | Canonical threshold | Canonical Fmax | Excluded Fmax | At original threshold | Status |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in report["run_rows"]:
        excluded_fmax = (
            f"{row['root_excluded_fmax']:.4f}"
            if row["root_excluded_fmax"] is not None
            else "n/a"
        )
        original_threshold_f = (
            f"{row['root_excluded_at_mode_canonical_threshold_f']:.4f}"
            if row["root_excluded_at_mode_canonical_threshold_f"] is not None
            else "n/a"
        )
        root_fraction = row["root_only_fraction_of_cafaeval_evaluable"]
        lines.append(
            "| {benchmark_id} | {mode} | {aspect} | {captured_rows:,} | "
            "{cafaeval_evaluable_targets:,} | {all_zero:,} | {root_fraction} | "
            "{canonical_threshold:.4f} | {canonical_fmax:.4f} | {excluded_fmax} | "
            "{original_threshold_f} | {exclusion_status} |".format(
                **row,
                root_fraction=(
                    f"{root_fraction:.2%}" if root_fraction is not None else "n/a"
                ),
                excluded_fmax=excluded_fmax,
                original_threshold_f=original_threshold_f,
            )
        )
    lines.extend(
        [
            "",
            "## Full Minus Sequence-Only",
            "",
            "| Benchmark | Aspect | Full threshold | Sequence threshold | Canonical delta | Excluded optimized delta | Delta at each mode's original threshold | Status |",
            "|---|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    if report["mode_delta_rows"]:
        for row in report["mode_delta_rows"]:
            excluded_delta = (
                f"{row['root_excluded_fmax_delta']:.4f}"
                if row["root_excluded_fmax_delta"] is not None
                else "n/a"
            )
            original_delta = (
                f"{row['root_excluded_at_mode_canonical_threshold_delta']:.4f}"
                if row["root_excluded_at_mode_canonical_threshold_delta"] is not None
                else "n/a"
            )
            lines.append(
                "| {benchmark_id} | {aspect} | {full_canonical_threshold:.4f} | "
                "{sequence_canonical_threshold:.4f} | {canonical_fmax_delta:.4f} | "
                "{excluded_delta} | {original_delta} | {status} |".format(
                    **row, excluded_delta=excluded_delta, original_delta=original_delta
                )
            )
    else:
        lines.append("| No benchmark has both modes yet |  |  |  |  |  |  |  |")
    lines.extend(
        [
            "",
            "The original-threshold value avoids retuning after cohort exclusion within each mode. Full and sequence-only may have different original thresholds, which are shown explicitly; this is not a shared-threshold comparison.",
            "These results compare model modes only within the same benchmark and aspect; they are not cross-benchmark model rankings.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", action="append", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if len(args.report) < 2:
        raise ValueError("At least two sensitivity reports are required")
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")
    loaded = [
        (path.resolve(), *load_report(path.resolve())) for path in args.report
    ]
    reports = [report for _, report, _ in loaded]
    identities = [(value["benchmark_id"], value["mode"]) for value in reports]
    if len(identities) != len(set(identities)):
        raise ValueError("Sensitivity reports repeat a benchmark/mode identity")

    run_rows: list[dict[str, Any]] = []
    indexed: dict[tuple[str, str, str], dict[str, Any]] = {}
    for report in reports:
        for aspect in report["selected_aspects"]:
            value = report["aspects"][aspect]
            excluded = value["root_only_excluded"]
            exclusion_complete = excluded.get("status") == "complete"
            cohorts = value["cohorts"]
            row = {
                "benchmark_id": report["benchmark_id"],
                "mode": report["mode"],
                "aspect": aspect,
                "captured_rows": cohorts["captured_rows"],
                "cafaeval_evaluable_targets": cohorts["cafaeval_evaluable_targets"],
                "all_zero": cohorts["all_zero"],
                "eligible_non_root": cohorts["eligible_non_root"],
                "root_only": cohorts["root_only"],
                "root_only_fraction_of_cafaeval_evaluable": cohorts[
                    "root_only_fraction_of_cafaeval_evaluable"
                ],
                "canonical_fmax": value["canonical_recheck"]["fmax"],
                "canonical_threshold": value["canonical_recheck"]["threshold"],
                "exclusion_status": excluded["status"],
                "root_excluded_fmax": (
                    excluded["cafa"]["fmax"] if exclusion_complete else None
                ),
                "root_excluded_at_mode_canonical_threshold_f": (
                    excluded["cafa"]["fixed_at_canonical_threshold"]["f"]
                    if exclusion_complete
                    else None
                ),
                "root_baseline_fmax": value["root_only_prediction_baseline"]["fmax"],
                "global_comparison_contract_sha256": report["comparison_contract"][
                    "sha256"
                ],
                "aspect_comparison_contract_sha256": value["comparison_contract"][
                    "sha256"
                ],
                "framework_commit": report["code_provenance"]["framework_commit"],
                "pfp_commit": report["code_provenance"]["pfp_commit"],
            }
            row["canonical_minus_root_baseline"] = (
                row["canonical_fmax"] - row["root_baseline_fmax"]
            )
            run_rows.append(row)
            indexed[(row["benchmark_id"], row["mode"], aspect)] = row

    delta_rows = []
    for benchmark_id in sorted({row["benchmark_id"] for row in run_rows}):
        for aspect in ASPECTS:
            full = indexed.get((benchmark_id, "full", aspect))
            sequence = indexed.get((benchmark_id, "sequence-only", aspect))
            if full is None or sequence is None:
                continue
            if (
                full["global_comparison_contract_sha256"]
                != sequence["global_comparison_contract_sha256"]
                or full["aspect_comparison_contract_sha256"]
                != sequence["aspect_comparison_contract_sha256"]
            ):
                raise ValueError(
                    f"Full and sequence-only provenance differs for {benchmark_id} {aspect}"
                )
            if full["exclusion_status"] != sequence["exclusion_status"]:
                raise ValueError(
                    f"Full and sequence-only exclusion status differs for {benchmark_id} {aspect}"
                )
            complete = full["exclusion_status"] == "complete"
            delta_rows.append(
                {
                    "benchmark_id": benchmark_id,
                    "aspect": aspect,
                    "canonical_fmax_delta": full["canonical_fmax"]
                    - sequence["canonical_fmax"],
                    "full_canonical_threshold": full["canonical_threshold"],
                    "sequence_canonical_threshold": sequence["canonical_threshold"],
                    "root_excluded_fmax_delta": (
                        full["root_excluded_fmax"] - sequence["root_excluded_fmax"]
                        if complete
                        else None
                    ),
                    "root_excluded_at_mode_canonical_threshold_delta": (
                        full["root_excluded_at_mode_canonical_threshold_f"]
                        - sequence["root_excluded_at_mode_canonical_threshold_f"]
                        if complete
                        else None
                    ),
                    "status": full["exclusion_status"],
                    "framework_commit_match": (
                        full["framework_commit"] == sequence["framework_commit"]
                    ),
                }
            )

    comparison = {
        "schema_version": 2,
        "status": "complete",
        "source_reports": [
            {"path": str(path), "sha256": snapshot["sha256"]}
            for path, _, snapshot in loaded
        ],
        "run_rows": sorted(
            run_rows, key=lambda row: (row["benchmark_id"], row["mode"], row["aspect"])
        ),
        "mode_delta_rows": delta_rows,
        "interpretation_policy": (
            "compare modes only within the same benchmark/aspect; preserve canonical "
            "results; mode-specific original-threshold deltas are no-retuning controls, "
            "not shared-threshold comparisons; framework commits are reported but are "
            "not part of the scientific pairing key"
        ),
    }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.stage-", dir=str(output_dir.parent))
    )
    try:
        atomic_write_json(stage / "label_sensitivity_comparison.json", comparison)
        atomic_write_text(
            stage / "label_sensitivity_comparison.md", markdown_report(comparison)
        )
        atomic_write_text(stage / "run_profiles.tsv", tsv_text(comparison["run_rows"]))
        atomic_write_text(
            stage / "full_vs_sequence_deltas.tsv",
            tsv_text(comparison["mode_delta_rows"]),
        )
        manifest = output_manifest(
            stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"}
        )
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
