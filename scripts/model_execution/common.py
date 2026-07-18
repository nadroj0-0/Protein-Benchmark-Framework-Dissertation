#!/usr/bin/env python3
"""Shared helpers for benchmark-agnostic PFP execution wrappers."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping


FRAMEWORK_ROOT = Path(__file__).resolve().parents[2]
ASPECTS = ("BPO", "CCO", "MFO")
ASPECT_TO_CSV = {"BPO": "bp", "CCO": "cc", "MFO": "mf"}
ASPECT_TO_NAMESPACE = {
    "BPO": "biological_process",
    "CCO": "cellular_component",
    "MFO": "molecular_function",
}
CSV_SPLITS = ("training", "validation", "test")
PFP_SPLITS = {"training": "train", "validation": "valid", "test": "test"}
MANDATORY_CAFA_METRICS = ("cafa_fmax", "cafa_wfmax", "cafa_smin")


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_run_config(path: Path) -> Dict[str, Any]:
    value = load_json(path)
    if not isinstance(value, dict) or value.get("schema_version") != 1:
        raise ValueError(f"Unsupported model-execution config: {path}")
    for key in ("benchmark_contract", "modalities", "training"):
        if not isinstance(value.get(key), dict):
            raise ValueError(f"Config is missing object {key!r}: {path}")
    expected_dimensions = {"sequence": 1024, "text": 768, "structure": 512, "ppi": 512}
    if set(value["modalities"]) != set(expected_dimensions):
        raise ValueError(f"Config must define exactly four PFP modalities: {path}")
    for modality, dimension in expected_dimensions.items():
        if int(value["modalities"][modality].get("dimension", -1)) != dimension:
            raise ValueError(f"Unexpected {modality} dimension in {path}")
    expected_training = {
        "seq_model": "prott5",
        "fusion_type": "gated_bilinear",
        "use_late_fusion": True,
        "late_output_mode": "hybrid",
        "hidden_dim": 512,
        "dropout": 0.4,
        "modality_dropout": 0.1,
        "aux_loss_weight": 0.8,
    }
    for key, expected in expected_training.items():
        if value["training"].get(key) != expected:
            raise ValueError(
                f"Run config changes frozen PFP behavior for {key}: "
                f"expected {expected!r}, found {value['training'].get(key)!r}"
            )
    evaluation = value.get("evaluation", {})
    if not isinstance(evaluation, dict) or not isinstance(
        evaluation.get("require_precomputed_ia", False), bool
    ):
        raise ValueError(f"Config evaluation.require_precomputed_ia must be boolean: {path}")
    reference = value.get("reference_preparation")
    if reference is not None:
        required_reference = {
            "name",
            "size_bytes",
            "checksum_algorithm",
            "checksum",
        }
        if not isinstance(reference, dict) or not required_reference.issubset(reference):
            raise ValueError(
                f"Config reference_preparation is incomplete: {path}"
            )
        if reference["checksum_algorithm"] not in hashlib.algorithms_available:
            raise ValueError(
                f"Unsupported reference checksum algorithm in {path}: "
                f"{reference['checksum_algorithm']}"
            )
        if not isinstance(reference["size_bytes"], int) or reference["size_bytes"] <= 0:
            raise ValueError(f"Invalid reference archive size in {path}")
        checksum = str(reference["checksum"])
        expected_length = hashlib.new(reference["checksum_algorithm"]).digest_size * 2
        if len(checksum) != expected_length or not re.fullmatch(r"[0-9a-f]+", checksum):
            raise ValueError(f"Invalid reference archive checksum in {path}")
    return value


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def atomic_write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_name, path)
    except BaseException:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def require_empty_output(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"Output directory must be absent or empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def selected_aspects(values: Iterable[str]) -> list[str]:
    result = list(values)
    if not result:
        result = list(ASPECTS)
    unknown = sorted(set(result) - set(ASPECTS))
    if unknown:
        raise ValueError(f"Unknown PFP aspects: {', '.join(unknown)}")
    if len(result) != len(set(result)):
        raise ValueError("Each PFP aspect may be selected only once")
    return result


def modality_paths(cache_root: Path, config: Mapping[str, Any]) -> Dict[str, Path]:
    paths: Dict[str, Path] = {}
    for name, spec in config["modalities"].items():
        directory = spec.get("directory")
        if not isinstance(directory, str) or not directory:
            raise ValueError(f"Modality {name!r} has no directory")
        path = Path(directory)
        paths[name] = path if path.is_absolute() else cache_root / path
    return paths


def expected_result_dir(output_base: Path, aspect: str) -> Path:
    return output_base / "fusion_comparison" / "prott5" / aspect / "gated_bilinear"


def validate_mandatory_metrics(metrics: Mapping[str, Any], context: str) -> None:
    missing = [key for key in MANDATORY_CAFA_METRICS if key not in metrics]
    if missing:
        raise ValueError(f"{context} lacks mandatory CAFA metrics: {missing}")
    values = {key: float(metrics[key]) for key in MANDATORY_CAFA_METRICS}
    non_finite = [key for key, value in values.items() if not math.isfinite(value)]
    if non_finite:
        raise ValueError(f"{context} has non-finite CAFA metrics: {non_finite}")
    bounded = [key for key in ("cafa_fmax", "cafa_wfmax") if not 0.0 <= values[key] <= 1.0]
    if bounded:
        raise ValueError(f"{context} has CAFA Fmax outside [0, 1]: {bounded}")
    if values["cafa_smin"] < 0.0:
        raise ValueError(f"{context} has negative cafa_smin")
