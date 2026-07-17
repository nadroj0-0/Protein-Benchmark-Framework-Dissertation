#!/usr/bin/env python3
"""Build a minimal PFP embedding workspace for contemporary retry pairs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from prepare_regeneration_workspace import validate_target_benchmark, write_workspace


def load_pair_ids(path: Path, modality: str) -> set[str]:
    if not path.is_file():
        raise ValueError(f"Missing pair table: {path}")
    result = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"protein_id", "modality"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(f"Invalid pair-table schema: {path}")
        for row in reader:
            if row["modality"] != modality:
                raise ValueError(
                    f"Pair table {path} contains {row['modality']}, expected {modality}"
                )
            if row["protein_id"] in result:
                raise ValueError(f"Repeated pair ID in {path}: {row['protein_id']}")
            result.add(row["protein_id"])
    return result


def load_plan_rows(plan_dir: Path) -> dict[str, dict[str, str]]:
    result = {}
    for action in ("reuse", "regenerate"):
        path = plan_dir / f"{action}_proteins.tsv"
        if not path.is_file():
            raise ValueError(f"Missing planner action table: {path}")
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "protein_id",
                "sequence",
                "sequence_sha256",
                "target_memberships",
                "action",
            }
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(f"Invalid planner action table: {path}")
            for row in reader:
                protein_id = row["protein_id"]
                if row["action"] != action:
                    raise ValueError(f"Action mismatch for {protein_id}: {row['action']}")
                if protein_id in result:
                    raise ValueError(f"Protein occurs in both planner actions: {protein_id}")
                result[protein_id] = row
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--target-benchmark-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--requested-pairs", type=Path, required=True)
    parser.add_argument("--control-pairs", type=Path, required=True)
    parser.add_argument(
        "--modality", choices=("sequence", "text", "structure", "ppi"), required=True
    )
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    target_hashes = validate_target_benchmark(
        args.plan_dir, args.target_benchmark_dir
    )
    requested = load_pair_ids(args.requested_pairs, args.modality)
    controls = load_pair_ids(args.control_pairs, args.modality)
    if requested & controls:
        raise SystemExit("Requested and subset-equivalence control IDs overlap")
    if not requested:
        raise SystemExit("Contemporary retry selection is empty")

    plan_rows = load_plan_rows(args.plan_dir)
    selected = requested | controls
    missing = sorted(selected - set(plan_rows))
    if missing:
        raise SystemExit(
            f"Retry IDs are absent from the bound planner tables: {missing[:10]} "
            f"(total={len(missing)})"
        )
    ordered = [plan_rows[protein_id] for protein_id in sorted(controls)]
    ordered.extend(plan_rows[protein_id] for protein_id in sorted(requested))
    workspace = write_workspace(ordered, args.data_dir)
    workspace.update(
        {
            "modality": args.modality,
            "requested_count": len(requested),
            "control_count": len(controls),
            "target_csv_sha256": target_hashes,
            "plan_dir": str(args.plan_dir.resolve()),
            "target_benchmark_dir": str(args.target_benchmark_dir.resolve()),
        }
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(workspace, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(workspace, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
