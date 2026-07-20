#!/usr/bin/env python3
"""Opt-in capture and publication of reusable PFP evaluation arrays."""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from common import atomic_write_json, sha256_file, sha256_lines


SCHEMA_VERSION = 1


def sha256_array(value: np.ndarray) -> str:
    """Hash an array's shape, dtype and canonical C-order bytes."""
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


def atomic_savez_compressed(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def directory_manifest(root: Path, excluded: Sequence[str] = ()) -> dict[str, Any]:
    omitted = set(excluded)
    files = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        if relative in omitted:
            continue
        files.append(
            {
                "path": relative,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "payload_file_count": len(files),
        "payload_bytes": sum(item["bytes"] for item in files),
        "files": files,
    }


class EvaluationArrayCapture:
    """Temporarily observe PFP's normal CAFA writers without changing their output."""

    def __init__(self, evaluation_module: Any, aspect: str, stage: Path) -> None:
        self.module = evaluation_module
        self.aspect = aspect
        self.stage = stage
        self.predictions: np.ndarray | None = None
        self.labels: np.ndarray | None = None
        self.protein_ids: list[str] | None = None
        self.go_terms: list[str] | None = None
        self.computed_ia_path: Path | None = None
        self._originals: dict[str, Any] = {}

    def __enter__(self) -> "EvaluationArrayCapture":
        names = (
            "save_predictions_cafa_format",
            "save_ground_truth_cafa_format",
            "save_ia_file",
        )
        missing = [name for name in names if not callable(getattr(self.module, name, None))]
        if missing:
            raise RuntimeError(
                "Prediction capture requires the standard PFP evaluation writers; "
                f"missing: {missing}"
            )
        self._originals = {name: getattr(self.module, name) for name in names}

        def save_predictions(
            predictions: np.ndarray,
            protein_ids: Sequence[str],
            go_terms: Sequence[str],
            output_file: str | Path,
        ) -> None:
            self.predictions = np.asarray(predictions)
            self.protein_ids = [str(value) for value in protein_ids]
            self.go_terms = [str(value) for value in go_terms]
            self._originals["save_predictions_cafa_format"](
                predictions, protein_ids, go_terms, output_file
            )

        def save_truth(
            labels: np.ndarray,
            protein_ids: Sequence[str],
            go_terms: Sequence[str],
            output_file: str | Path,
        ) -> None:
            observed_ids = [str(value) for value in protein_ids]
            observed_terms = [str(value) for value in go_terms]
            if self.protein_ids is not None and observed_ids != self.protein_ids:
                raise ValueError("PFP prediction and truth protein order differ")
            if self.go_terms is not None and observed_terms != self.go_terms:
                raise ValueError("PFP prediction and truth GO-term order differ")
            self.labels = np.asarray(labels)
            self.protein_ids = observed_ids
            self.go_terms = observed_terms
            self._originals["save_ground_truth_cafa_format"](
                labels, protein_ids, go_terms, output_file
            )

        def save_ia(values: Mapping[str, float], output_file: str | Path) -> None:
            self._originals["save_ia_file"](values, output_file)
            source = Path(output_file)
            destination = self.stage / f"{self.aspect}_ia.txt"
            shutil.copyfile(source, destination)
            self.computed_ia_path = destination

        self.module.save_predictions_cafa_format = save_predictions
        self.module.save_ground_truth_cafa_format = save_truth
        self.module.save_ia_file = save_ia
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        for name, value in self._originals.items():
            setattr(self.module, name, value)

    def persist(
        self,
        *,
        expected_protein_ids: Sequence[str],
        expected_go_terms: Sequence[str],
        checkpoint: Path,
        expected_checkpoint_sha256: str,
        cafa_metrics: Mapping[str, Any],
        ia_file: Path | None,
        expected_ia_sha256: str | None,
    ) -> dict[str, Any]:
        if self.predictions is None or self.labels is None:
            raise RuntimeError(f"PFP did not expose prediction and truth arrays for {self.aspect}")
        protein_ids = [str(value) for value in expected_protein_ids]
        go_terms = [str(value) for value in expected_go_terms]
        if self.protein_ids != protein_ids or self.go_terms != go_terms:
            raise ValueError(f"Captured PFP order differs from prepared data for {self.aspect}")

        scores = np.asarray(self.predictions)
        labels = np.asarray(self.labels)
        if scores.shape != labels.shape or scores.shape != (len(protein_ids), len(go_terms)):
            raise ValueError(
                f"Captured array shape mismatch for {self.aspect}: "
                f"scores={scores.shape}, labels={labels.shape}"
            )
        if scores.dtype.kind != "f":
            raise ValueError(f"Captured prediction scores are not floating point for {self.aspect}")
        if not np.isfinite(scores).all() or np.any(scores < 0) or np.any(scores > 1):
            raise ValueError(f"Captured prediction scores are invalid for {self.aspect}")
        if not np.isin(labels, (0, 1)).all():
            raise ValueError(f"Captured truth is non-binary for {self.aspect}")
        truth = labels.astype(np.uint8, copy=False)
        if len(set(protein_ids)) != len(protein_ids):
            raise ValueError(f"Captured protein IDs are not unique for {self.aspect}")
        if len(set(go_terms)) != len(go_terms):
            raise ValueError(f"Captured GO terms are not unique for {self.aspect}")

        threshold = cafa_metrics.get("threshold")
        if not isinstance(threshold, (int, float, np.integer, np.floating)) or not math.isfinite(
            float(threshold)
        ):
            raise ValueError(
                f"Prediction capture requires PFP's canonical CAFA threshold for {self.aspect}"
            )

        array_path = self.stage / f"{self.aspect}_evaluation_arrays.npz"
        atomic_savez_compressed(
            array_path,
            scores=scores,
            truth=truth,
            protein_ids=np.asarray(protein_ids, dtype=str),
            go_terms=np.asarray(go_terms, dtype=str),
        )

        ia_source = "computed-from-training-labels"
        effective_ia = self.computed_ia_path
        if ia_file is not None:
            if expected_ia_sha256 is None or sha256_file(ia_file) != expected_ia_sha256:
                raise ValueError(f"Precomputed IA changed during evaluation for {self.aspect}")
            effective_ia = self.stage / f"{self.aspect}_ia.txt"
            shutil.copyfile(ia_file, effective_ia)
            ia_source = "precomputed"
        if effective_ia is None or not effective_ia.is_file():
            raise ValueError(
                f"Prediction capture requires the exact IA file used for {self.aspect}"
            )
        if ia_file is not None and sha256_file(effective_ia) != expected_ia_sha256:
            raise ValueError(f"Copied IA differs from its source for {self.aspect}")
        if sha256_file(checkpoint) != expected_checkpoint_sha256:
            raise ValueError(f"Checkpoint changed during evaluation for {self.aspect}")

        numeric_metrics = {
            key: float(value)
            for key, value in cafa_metrics.items()
            if isinstance(value, (int, float, np.integer, np.floating))
            and math.isfinite(float(value))
        }
        return {
            "aspect": self.aspect,
            "array_file": array_path.name,
            "array_file_bytes": array_path.stat().st_size,
            "array_file_sha256": sha256_file(array_path),
            "shape": list(scores.shape),
            "scores_dtype": scores.dtype.str,
            "truth_dtype": truth.dtype.str,
            "scores_content_sha256": sha256_array(scores),
            "truth_content_sha256": sha256_array(truth),
            "protein_ids_sha256": sha256_lines(protein_ids),
            "go_terms_sha256": sha256_lines(go_terms),
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": expected_checkpoint_sha256,
            "ia_file": effective_ia.name,
            "ia_file_sha256": sha256_file(effective_ia),
            "ia_source": ia_source,
            "canonical_cafa_metrics": numeric_metrics,
        }


def publish_prediction_artifacts(
    stage: Path,
    destination: Path,
    manifest: dict[str, Any],
) -> None:
    if destination.exists():
        raise ValueError(f"Prediction artifact directory already exists: {destination}")
    atomic_write_json(stage / "prediction_artifact_manifest.json", manifest)
    output = directory_manifest(
        stage, excluded=("output_manifest.json", "RUN_COMPLETE.json")
    )
    atomic_write_json(stage / "output_manifest.json", output)
    atomic_write_json(
        stage / "RUN_COMPLETE.json",
        {
            "schema_version": SCHEMA_VERSION,
            "complete": True,
            "benchmark_id": manifest["benchmark_id"],
            "mode": manifest["mode"],
            "output_manifest_sha256": sha256_file(stage / "output_manifest.json"),
        },
    )
    os.replace(stage, destination)
