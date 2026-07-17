#!/usr/bin/env python3
"""Measure same-node embedding repeatability without choosing a tolerance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


COMPARISONS = (
    ("repeat_1_vs_repeat_2", "repeat_1", "repeat_2"),
    ("baseline_vs_repeat_1", "baseline", "repeat_1"),
    ("baseline_vs_repeat_2", "baseline", "repeat_2"),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_controls(path: Path, modality: str) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"protein_id", "modality"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(f"Invalid control table: {path}")
        result = []
        seen = set()
        for row in reader:
            if row["modality"] != modality:
                raise ValueError(
                    f"Control {row['protein_id']} has modality {row['modality']}, "
                    f"expected {modality}"
                )
            protein_id = row["protein_id"]
            if not protein_id or protein_id in seen:
                raise ValueError(f"Invalid or repeated control ID: {protein_id!r}")
            seen.add(protein_id)
            result.append(protein_id)
    if not result:
        raise ValueError("Control table is empty")
    return result


def load_vector(path: Path, dimension: int) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(path)
    value = np.squeeze(np.load(path, allow_pickle=False))
    if value.shape != (dimension,):
        raise ValueError(f"wrong shape {value.shape}, expected {(dimension,)}")
    if not np.isfinite(value).all():
        raise ValueError("array contains non-finite values")
    return np.asarray(value, dtype=np.float64)


def vector_metrics(left: np.ndarray, right: np.ndarray, rtol: float, atol: float) -> dict:
    difference = right - left
    absolute = np.abs(difference)
    left_norm = float(np.linalg.norm(left))
    right_norm = float(np.linalg.norm(right))
    denominator = left_norm * right_norm
    cosine = None
    if denominator > 0:
        cosine = float(np.dot(left, right) / denominator)
        cosine = max(-1.0, min(1.0, cosine))
    elif np.array_equal(left, right):
        cosine = 1.0
    relative = absolute / np.maximum(np.abs(left), np.finfo(np.float64).eps)
    return {
        "exact_equal": bool(np.array_equal(left, right)),
        "allclose": bool(np.allclose(left, right, rtol=rtol, atol=atol)),
        "max_abs_difference": float(np.max(absolute)),
        "mean_abs_difference": float(np.mean(absolute)),
        "rmse": float(math.sqrt(float(np.mean(np.square(difference))))),
        "l2_difference": float(np.linalg.norm(difference)),
        "max_relative_difference": float(np.max(relative)),
        "cosine_similarity": cosine,
        "left_l2_norm": left_norm,
        "right_l2_norm": right_norm,
    }


def aggregate(rows: list[dict]) -> dict:
    valid = [row for row in rows if row["status"] == "compared"]
    result = {
        "requested": len(rows),
        "compared": len(valid),
        "integrity_failures": len(rows) - len(valid),
    }
    if not valid:
        return result
    max_abs = [row["max_abs_difference"] for row in valid]
    mean_abs = [row["mean_abs_difference"] for row in valid]
    rmse = [row["rmse"] for row in valid]
    cosines = [row["cosine_similarity"] for row in valid if row["cosine_similarity"] is not None]
    result.update(
        {
            "exact_equal": sum(bool(row["exact_equal"]) for row in valid),
            "allclose": sum(bool(row["allclose"]) for row in valid),
            "max_abs_difference_max": max(max_abs),
            "max_abs_difference_median": float(np.median(max_abs)),
            "mean_abs_difference_max": max(mean_abs),
            "rmse_max": max(rmse),
            "cosine_similarity_min": min(cosines) if cosines else None,
        }
    )
    return result


def markdown_report(report: dict) -> str:
    lines = [
        f"# {report['modality'].title()} Embedding Reproducibility Diagnostic",
        "",
        "This diagnostic records numerical repeatability. It does not merge arrays or "
        "automatically choose an acceptance tolerance.",
        "",
        f"- Controls: `{report['control_count']}`",
        f"- Integrity passed: `{str(report['integrity_passed']).lower()}`",
        f"- Existing comparison tolerance: `rtol={report['rtol']}`, "
        f"`atol={report['atol']}`",
        "",
        "## Summary",
        "",
        "| Comparison | Compared | Exact | Existing allclose | Maximum absolute difference | Minimum cosine similarity |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, _, _ in COMPARISONS:
        value = report["summaries"][name]
        lines.append(
            "| {name} | {compared} | {exact} | {close} | {maximum:.12g} | {cosine} |".format(
                name=name.replace("_", " "),
                compared=value["compared"],
                exact=value.get("exact_equal", 0),
                close=value.get("allclose", 0),
                maximum=value.get("max_abs_difference_max", float("nan")),
                cosine=(
                    f"{value['cosine_similarity_min']:.12g}"
                    if value.get("cosine_similarity_min") is not None
                    else "n/a"
                ),
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "The repeat-1 versus repeat-2 distribution measures same-job, same-node "
            "numerical variation. Baseline comparisons additionally include any "
            "historical hardware or input differences. A scientific tolerance must be "
            "chosen after inspecting these distributions; this report does not silently "
            "relax the production gate.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--controls", type=Path, required=True)
    parser.add_argument("--modality", choices=("text", "structure"), required=True)
    parser.add_argument("--baseline-root", type=Path, required=True)
    parser.add_argument("--repeat-one-root", type=Path, required=True)
    parser.add_argument("--repeat-two-root", type=Path, required=True)
    parser.add_argument("--input-file", type=Path, action="append", default=[])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--rtol", type=float, default=1e-5)
    parser.add_argument("--atol", type=float, default=1e-6)
    parser.add_argument("--minimum-compared", type=int, default=20)
    args = parser.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    specification = contract["policy"]["modalities"][args.modality]
    directory = specification["cache_directory"]
    dimension = int(specification["dimension"])
    controls = load_controls(args.controls, args.modality)
    roots = {
        "baseline": args.baseline_root,
        "repeat_1": args.repeat_one_root,
        "repeat_2": args.repeat_two_root,
    }

    comparison_rows = []
    by_comparison = {name: [] for name, _, _ in COMPARISONS}
    for comparison, left_name, right_name in COMPARISONS:
        for protein_id in controls:
            row = {
                "comparison": comparison,
                "protein_id": protein_id,
                "status": "compared",
                "detail": "",
            }
            left_path = roots[left_name] / directory / f"{protein_id}.npy"
            right_path = roots[right_name] / directory / f"{protein_id}.npy"
            try:
                left = load_vector(left_path, dimension)
                right = load_vector(right_path, dimension)
                row.update(vector_metrics(left, right, args.rtol, args.atol))
            except Exception as error:
                row["status"] = "integrity_failure"
                row["detail"] = f"{type(error).__name__}:{error}"
            comparison_rows.append(row)
            by_comparison[comparison].append(row)

    input_rows = []
    seen_inputs = set()
    for path in args.input_file:
        resolved = path.resolve()
        if resolved in seen_inputs:
            raise SystemExit(f"Repeated input file: {resolved}")
        seen_inputs.add(resolved)
        if not resolved.is_file():
            raise SystemExit(f"Missing reproducibility input: {resolved}")
        input_rows.append(
            {
                "path": str(resolved),
                "size_bytes": resolved.stat().st_size,
                "sha256": sha256_file(resolved),
            }
        )

    summaries = {name: aggregate(rows) for name, rows in by_comparison.items()}
    integrity_passed = all(
        summary["integrity_failures"] == 0
        and summary["compared"] >= args.minimum_compared
        for summary in summaries.values()
    )
    report = {
        "schema_version": 1,
        "modality": args.modality,
        "control_count": len(controls),
        "dimension": dimension,
        "cache_directory": directory,
        "rtol": args.rtol,
        "atol": args.atol,
        "minimum_compared": args.minimum_compared,
        "integrity_passed": integrity_passed,
        "summaries": summaries,
        "inputs": input_rows,
        "rows": comparison_rows,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "embedding_reproducibility.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    columns = [
        "comparison",
        "protein_id",
        "status",
        "exact_equal",
        "allclose",
        "max_abs_difference",
        "mean_abs_difference",
        "rmse",
        "l2_difference",
        "max_relative_difference",
        "cosine_similarity",
        "left_l2_norm",
        "right_l2_norm",
        "detail",
    ]
    with (args.output_dir / "embedding_reproducibility.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in comparison_rows:
            writer.writerow({column: row.get(column, "") for column in columns})
    with (args.output_dir / "input_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("path", "size_bytes", "sha256"), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(input_rows)
    (args.output_dir / "embedding_reproducibility.md").write_text(
        markdown_report(report), encoding="utf-8"
    )
    print(json.dumps(report["summaries"], indent=2, sort_keys=True))
    return 0 if integrity_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
