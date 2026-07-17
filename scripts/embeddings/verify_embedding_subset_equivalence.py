#!/usr/bin/env python3
"""Verify that subset generation reproduces accepted per-protein arrays."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


DEFAULT_RTOL = 1e-5
DEFAULT_ATOL = 1e-6
REGENERATED_EMBEDDING_ATOL = 1e-4


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-root", type=Path, required=True)
    parser.add_argument("--reference-cache-root", type=Path)
    parser.add_argument("--generated-cache-root", type=Path, required=True)
    parser.add_argument("--control-pairs", type=Path, required=True)
    parser.add_argument("--modality", required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--rtol", type=float)
    parser.add_argument("--atol", type=float)
    parser.add_argument("--minimum-compared", type=int, default=5)
    args = parser.parse_args()

    contract = json.loads((args.state_root / "contract.json").read_text(encoding="utf-8"))
    specification = contract["policy"]["modalities"].get(args.modality)
    if specification is None:
        raise SystemExit(f"Unknown modality: {args.modality}")
    rtol = DEFAULT_RTOL if args.rtol is None else args.rtol
    default_atol = (
        REGENERATED_EMBEDDING_ATOL
        if args.modality in {"text", "structure", "ppi"}
        else DEFAULT_ATOL
    )
    atol = default_atol if args.atol is None else args.atol
    directory = specification["cache_directory"]
    dimension = int(specification["dimension"])
    reference_root = args.reference_cache_root or (args.state_root / "cache")
    rows = []
    with args.control_pairs.open(encoding="utf-8", newline="") as handle:
        controls = list(csv.DictReader(handle, delimiter="\t"))
    if not controls:
        raise SystemExit("No accepted controls are available for subset-equivalence testing")

    failed = 0
    compared = 0
    unavailable = 0
    equivalent = 0
    for control in controls:
        protein_id = control["protein_id"]
        reference_path = reference_root / directory / f"{protein_id}.npy"
        generated_path = args.generated_cache_root / directory / f"{protein_id}.npy"
        status = "equivalent"
        detail = ""
        max_abs = None
        if not reference_path.is_file():
            status = "missing_reference"
        elif not generated_path.is_file():
            status = "unavailable_in_control_run"
            unavailable += 1
        else:
            try:
                reference = np.squeeze(np.load(reference_path, allow_pickle=False))
                generated = np.squeeze(np.load(generated_path, allow_pickle=False))
                if reference.shape != (dimension,) or generated.shape != (dimension,):
                    status = "wrong_shape"
                    detail = f"reference={reference.shape};generated={generated.shape}"
                elif not np.isfinite(reference).all() or not np.isfinite(generated).all():
                    status = "non_finite"
                else:
                    compared += 1
                    max_abs = float(np.max(np.abs(reference - generated)))
                    if not np.allclose(reference, generated, rtol=rtol, atol=atol):
                        status = "different"
                    else:
                        equivalent += 1
            except Exception as error:
                status = "unreadable"
                detail = f"{type(error).__name__}:{error}"
        if status not in {"equivalent", "unavailable_in_control_run"}:
            failed += 1
        rows.append(
            {
                "protein_id": protein_id,
                "status": status,
                "max_abs_difference": max_abs,
                "detail": detail,
            }
        )

    report = {
        "schema_version": 1,
        "modality": args.modality,
        "control_count": len(rows),
        "equivalent": equivalent,
        "compared": compared,
        "unavailable": unavailable,
        "failed": failed,
        "minimum_compared": args.minimum_compared,
        "rtol": rtol,
        "atol": atol,
        "tolerance_source": {
            "rtol": "cli" if args.rtol is not None else "default",
            "atol": (
                "cli"
                if args.atol is not None
                else (
                    "regenerated_embedding_compatibility"
                    if args.modality in {"text", "structure", "ppi"}
                    else "default"
                )
            ),
        },
        "rows": rows,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if failed == 0 and compared >= args.minimum_compared else 1


if __name__ == "__main__":
    raise SystemExit(main())
