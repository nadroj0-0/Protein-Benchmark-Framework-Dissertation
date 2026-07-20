#!/usr/bin/env python3
"""Shared validation and scoring helpers for PFP label sensitivity analyses."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from label_space_common import (
    ASPECT_TO_ROOT,
    file_snapshot,
    require_unchanged,
    sha256_file,
    sha256_json,
)


def sha256_lines(values: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def sha256_array(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(
        json.dumps(
            {"dtype": array.dtype.str, "shape": list(array.shape)},
            sort_keys=True,
            separators=(",", ":"),
        ).encode("ascii")
    )
    digest.update(b"\0")
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def verify_artifact_manifest(manifest_path: Path) -> tuple[dict[str, Any], Path]:
    root = manifest_path.resolve().parent
    completion_path = root / "RUN_COMPLETE.json"
    output_manifest_path = root / "output_manifest.json"
    manifest_snapshot = file_snapshot(manifest_path)
    completion_snapshot = file_snapshot(completion_path)
    output_snapshot = file_snapshot(output_manifest_path)
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    if not completion.get("complete"):
        raise ValueError(f"Prediction artifact completion marker is false: {root}")
    if completion.get("output_manifest_sha256") != sha256_file(output_manifest_path):
        raise ValueError(f"Prediction output manifest is not bound by completion marker: {root}")
    output = json.loads(output_manifest_path.read_text(encoding="utf-8"))
    listed = {item.get("path"): item for item in output.get("files", [])}
    if manifest_path.name not in listed:
        raise ValueError(f"Prediction output manifest does not bind its report: {manifest_path}")
    if (
        listed[manifest_path.name].get("bytes") != manifest_snapshot["bytes"]
        or listed[manifest_path.name].get("sha256") != manifest_snapshot["sha256"]
    ):
        raise ValueError(f"Prediction output manifest binds different report bytes: {manifest_path}")
    for item in output.get("files", []):
        relative = Path(item["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Prediction output manifest contains unsafe path: {relative}")
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"Prediction artifact is missing: {path}")
        if path.stat().st_size != item["bytes"] or sha256_file(path) != item["sha256"]:
            raise ValueError(f"Prediction artifact changed after publication: {path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 2 or manifest.get("status") != "complete":
        raise ValueError(f"Unsupported or incomplete prediction manifest: {manifest_path}")
    selected = manifest.get("selected_aspects")
    if (
        not isinstance(selected, list)
        or len(selected) != len(set(selected))
        or set(selected) != set(manifest.get("aspects", {}))
    ):
        raise ValueError(f"Prediction manifest aspect inventory is inconsistent: {manifest_path}")
    if manifest.get("mode") not in {"full", "sequence-only"}:
        raise ValueError(f"Prediction manifest has unsupported mode: {manifest_path}")
    provenance = manifest.get("provenance")
    if not isinstance(provenance, dict) or not provenance.get("benchmark_fingerprint"):
        raise ValueError(f"Prediction manifest lacks benchmark provenance: {manifest_path}")
    for field in ("preparation_report", "embedding_validation_report"):
        specification = provenance.get(field)
        if not isinstance(specification, dict):
            raise ValueError(f"Prediction manifest lacks {field}: {manifest_path}")
        relative = Path(str(specification.get("artifact_file", "")))
        if not relative.name or relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Prediction manifest has unsafe {field} path")
        artifact = root / relative
        if (
            not artifact.is_file()
            or artifact.stat().st_size != specification.get("bytes")
            or sha256_file(artifact) != specification.get("sha256")
        ):
            raise ValueError(f"Prediction {field} differs from provenance binding")
    require_unchanged(manifest_path, manifest_snapshot, "Prediction manifest")
    require_unchanged(completion_path, completion_snapshot, "Prediction completion marker")
    require_unchanged(output_manifest_path, output_snapshot, "Prediction output manifest")
    return manifest, root


def global_comparison_contract(manifest: Mapping[str, Any]) -> dict[str, Any]:
    provenance = manifest["provenance"]
    value = {
        "benchmark_id": manifest["benchmark_id"],
        "benchmark_fingerprint": provenance["benchmark_fingerprint"],
        "source_csv_sha256": provenance.get("source_csv_sha256", {}),
        "seed": manifest["seed"],
        "config_sha256": manifest["config"]["sha256"],
        "obo_sha256": manifest["obo"]["sha256"],
        "pfp_commit": provenance["pfp_commit"],
    }
    return {"value": value, "sha256": sha256_json(value)}


def aspect_comparison_contract(specification: Mapping[str, Any]) -> dict[str, Any]:
    value = {
        key: specification[key]
        for key in (
            "truth_content_sha256",
            "protein_ids_sha256",
            "go_terms_sha256",
            "ia_file_sha256",
        )
    }
    return {"value": value, "sha256": sha256_json(value)}


def cafaeval_contract() -> dict[str, Any]:
    try:
        version = importlib.metadata.version("cafaeval")
    except importlib.metadata.PackageNotFoundError:
        version = "unavailable-package-metadata"
    return {
        "package": "cafaeval",
        "version": version,
        "arguments": {
            "ia": "captured-exact-file",
            "no_orphans": False,
            "norm": "cafa",
            "prop": "max",
        },
    }


def load_aspect_bundle(
    manifest: Mapping[str, Any], root: Path, aspect: str
) -> dict[str, Any]:
    try:
        specification = manifest["aspects"][aspect]
    except KeyError as exc:
        raise ValueError(f"Prediction manifest has no {aspect} artifact") from exc
    array_relative = Path(specification["array_file"])
    ia_relative = Path(specification["ia_file"])
    for relative in (array_relative, ia_relative):
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Prediction manifest contains unsafe path: {relative}")
    array_path = root / array_relative
    array_snapshot = file_snapshot(array_path)
    if (
        array_snapshot["bytes"] != specification["array_file_bytes"]
        or array_snapshot["sha256"] != specification["array_file_sha256"]
    ):
        raise ValueError(f"Prediction array file hash differs for {aspect}")
    with np.load(array_path, allow_pickle=False) as archive:
        required = {"scores", "truth", "protein_ids", "go_terms"}
        if set(archive.files) != required:
            raise ValueError(
                f"Prediction artifact keys differ for {aspect}: {sorted(archive.files)}"
            )
        scores = np.asarray(archive["scores"])
        truth = np.asarray(archive["truth"])
        protein_ids = [str(value) for value in archive["protein_ids"].tolist()]
        go_terms = [str(value) for value in archive["go_terms"].tolist()]
    require_unchanged(array_path, array_snapshot, f"{aspect} prediction array")
    shape = tuple(specification["shape"])
    if scores.shape != truth.shape or scores.shape != shape:
        raise ValueError(f"Prediction/truth shape differs for {aspect}")
    if scores.ndim != 2 or scores.shape != (len(protein_ids), len(go_terms)):
        raise ValueError(f"Prediction metadata dimensions differ for {aspect}")
    if not np.isfinite(scores).all() or np.any(scores < 0) or np.any(scores > 1):
        raise ValueError(f"Prediction scores are invalid for {aspect}")
    if not np.isin(truth, (0, 1)).all():
        raise ValueError(f"Prediction truth is non-binary for {aspect}")
    if len(set(protein_ids)) != len(protein_ids) or len(set(go_terms)) != len(go_terms):
        raise ValueError(f"Prediction IDs or GO terms are duplicated for {aspect}")
    expected_hashes = {
        "scores_content_sha256": sha256_array(scores),
        "truth_content_sha256": sha256_array(truth.astype(np.uint8, copy=False)),
        "protein_ids_sha256": sha256_lines(protein_ids),
        "go_terms_sha256": sha256_lines(go_terms),
    }
    for field, observed in expected_hashes.items():
        if observed != specification[field]:
            raise ValueError(f"Prediction {field} differs for {aspect}")
    ia_path = root / ia_relative
    ia_snapshot = file_snapshot(ia_path)
    if ia_snapshot["sha256"] != specification["ia_file_sha256"]:
        raise ValueError(f"Prediction IA file differs for {aspect}")
    require_unchanged(ia_path, ia_snapshot, f"{aspect} IA file")
    required_root = ASPECT_TO_ROOT[aspect]
    if required_root not in go_terms:
        raise ValueError(f"Prediction GO terms lack {aspect} root {required_root}")
    return {
        "specification": dict(specification),
        "scores": scores,
        "truth": truth.astype(np.uint8, copy=False),
        "protein_ids": protein_ids,
        "go_terms": go_terms,
        "ia_path": ia_path,
        "root_index": go_terms.index(required_root),
    }


def cohort_masks(truth: np.ndarray, root_index: int) -> dict[str, np.ndarray]:
    non_root = np.ones(truth.shape[1], dtype=bool)
    non_root[root_index] = False
    has_non_root = truth[:, non_root].any(axis=1)
    has_any = truth.any(axis=1)
    root_only = (truth[:, root_index] == 1) & ~has_non_root
    return {
        "eligible_non_root": has_non_root,
        "root_only": root_only,
        "all_zero": ~has_any,
        "non_root_columns": non_root,
    }


def _threshold_metrics(truth: np.ndarray, scores: np.ndarray, threshold: float) -> dict[str, float]:
    predicted = scores >= threshold
    true = truth.astype(bool, copy=False)
    true_positive = np.logical_and(predicted, true).sum(axis=1)
    predicted_count = predicted.sum(axis=1)
    true_count = true.sum(axis=1)
    covered = predicted_count > 0
    precision = (
        float(np.mean(true_positive[covered] / predicted_count[covered]))
        if covered.any()
        else 0.0
    )
    recall = float(np.mean(true_positive / true_count)) if len(true_count) else 0.0
    macro_f = (
        2.0 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    tp = int(true_positive.sum())
    fp = int(np.logical_and(predicted, ~true).sum())
    fn = int(np.logical_and(~predicted, true).sum())
    micro_precision = tp / (tp + fp) if tp + fp else 0.0
    micro_recall = tp / (tp + fn) if tp + fn else 0.0
    micro_f = (
        2.0 * micro_precision * micro_recall / (micro_precision + micro_recall)
        if micro_precision + micro_recall
        else 0.0
    )
    return {
        "threshold": float(threshold),
        "macro_precision": precision,
        "macro_recall": recall,
        "macro_f": macro_f,
        "micro_precision": micro_precision,
        "micro_recall": micro_recall,
        "micro_f": micro_f,
        "coverage": float(np.mean(covered)) if len(covered) else 0.0,
    }


def flat_non_root_metrics(
    truth: np.ndarray,
    scores: np.ndarray,
    root_index: int,
    fixed_threshold: float,
) -> dict[str, Any]:
    masks = cohort_masks(truth, root_index)
    eligible = masks["eligible_non_root"]
    columns = masks["non_root_columns"]
    cohort_truth = truth[eligible][:, columns]
    cohort_scores = scores[eligible][:, columns]
    if cohort_truth.shape[0] == 0:
        return {"status": "not_evaluable_no_non_root_targets"}
    thresholds = sorted(
        set(np.linspace(0.01, 0.99, 100).tolist() + [float(fixed_threshold)])
    )
    rows = [_threshold_metrics(cohort_truth, cohort_scores, value) for value in thresholds]
    best = max(rows, key=lambda row: (row["macro_f"], -row["threshold"]))
    fixed = min(rows, key=lambda row: abs(row["threshold"] - fixed_threshold))
    return {
        "status": "complete",
        "policy": (
            "protein-centric flat metrics after removing the ontology root; "
            "no GO propagation; diagnostic only, not a CAFA metric"
        ),
        "best": best,
        "fixed_at_canonical_threshold": fixed,
    }


def write_prediction_file(
    path: Path,
    protein_ids: Sequence[str],
    go_terms: Sequence[str],
    scores: np.ndarray | None,
    root_index: int | None = None,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        if scores is None:
            if root_index is None:
                raise ValueError("Root-only predictions require a root index")
            root = go_terms[root_index]
            for protein_id in protein_ids:
                handle.write(f"{protein_id}\t{root}\t1.000000\n")
            return
        lines: list[str] = []
        for row_index, protein_id in enumerate(protein_ids):
            for term_index in np.flatnonzero(scores[row_index] > 0):
                lines.append(
                    f"{protein_id}\t{go_terms[int(term_index)]}\t"
                    f"{float(scores[row_index, term_index]):.6f}\n"
                )
                if len(lines) >= 16_384:
                    handle.write("".join(lines))
                    lines.clear()
        if lines:
            handle.write("".join(lines))


def write_truth_file(
    path: Path,
    protein_ids: Sequence[str],
    go_terms: Sequence[str],
    truth: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        lines: list[str] = []
        for row_index, protein_id in enumerate(protein_ids):
            for term_index in np.flatnonzero(truth[row_index] == 1):
                lines.append(f"{protein_id}\t{go_terms[int(term_index)]}\n")
                if len(lines) >= 16_384:
                    handle.write("".join(lines))
                    lines.clear()
        if lines:
            handle.write("".join(lines))


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_cafa_results(result_dir: Path, fixed_threshold: float) -> dict[str, Any]:
    all_path = result_dir / "evaluation_all.tsv"
    best_s_path = result_dir / "evaluation_best_s.tsv"
    best_w_path = result_dir / "evaluation_best_f_w.tsv"
    for path in (all_path, best_s_path, best_w_path):
        if not path.is_file():
            raise FileNotFoundError(f"cafaeval did not produce required result: {path}")
    rows = _read_tsv(all_path)
    if not rows:
        raise ValueError(f"cafaeval produced no threshold rows: {all_path}")
    best = max(rows, key=lambda row: float(row["f"]))
    fixed = min(rows, key=lambda row: abs(float(row["tau"]) - fixed_threshold))
    if not math.isclose(float(fixed["tau"]), fixed_threshold, abs_tol=1e-12):
        raise ValueError(
            f"Canonical threshold {fixed_threshold} is absent from cafaeval output"
        )
    best_s = _read_tsv(best_s_path)[0]
    best_w = _read_tsv(best_w_path)[0]

    def selected(row: Mapping[str, str], fields: Sequence[str]) -> dict[str, float]:
        return {
            field: float(row[field])
            for field in fields
            if field in row and row[field] not in ("", None)
        }

    return {
        "fmax": float(best["f"]),
        "threshold": float(best["tau"]),
        "precision": float(best["pr"]),
        "recall": float(best["rc"]),
        "coverage": float(best["cov"]) if "cov" in best else None,
        "wfmax": float(best_w["f_w"]),
        "wthreshold": float(best_w["tau"]),
        "wprecision": float(best_w["pr_w"]),
        "wrecall": float(best_w["rc_w"]),
        "smin": float(best_s["s"]),
        "fixed_at_canonical_threshold": selected(
            fixed, ("tau", "f", "pr", "rc", "cov", "f_w", "pr_w", "rc_w", "s")
        ),
    }


def run_cafa_evaluation(
    *,
    obo_file: Path,
    ia_file: Path,
    destination: Path,
    protein_ids: Sequence[str],
    go_terms: Sequence[str],
    truth: np.ndarray,
    fixed_threshold: float,
    scores: np.ndarray | None,
    root_index: int | None = None,
) -> dict[str, Any]:
    from cafaeval.evaluation import cafa_eval, write_results

    started = time.perf_counter()
    destination.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(prefix="pfp-label-sensitivity-") as name:
        temporary = Path(name)
        prediction_dir = temporary / "predictions"
        prediction_dir.mkdir()
        write_prediction_file(
            prediction_dir / "model.tsv",
            protein_ids,
            go_terms,
            scores,
            root_index=root_index,
        )
        truth_file = temporary / "truth.tsv"
        write_truth_file(truth_file, protein_ids, go_terms, truth)
        results = cafa_eval(
            str(obo_file),
            str(prediction_dir),
            str(truth_file),
            ia=str(ia_file),
            no_orphans=False,
            norm="cafa",
            prop="max",
        )
        write_results(*results, out_dir=str(destination))
    parsed = parse_cafa_results(destination, fixed_threshold)
    parsed["wall_seconds"] = time.perf_counter() - started
    return parsed
