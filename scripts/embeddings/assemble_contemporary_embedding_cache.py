#!/usr/bin/env python3
"""Assemble and validate a PFP cache from reuse and regeneration sources."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
from collections import Counter
from pathlib import Path

import numpy as np


MODALITIES = {
    "prott5": {"directory": "prott5", "dimension": 1024},
    "text": {"directory": "exp_text_embeddings_temporal", "dimension": 768},
    "structure": {"directory": "IF1", "dimension": 512},
    "ppi": {"directory": "ppi", "dimension": 512},
}


def load_actions(plan_dir: Path) -> tuple[dict[str, dict], dict[str, dict]]:
    actions = {}
    by_action = {"reuse": {}, "regenerate": {}}
    for action in by_action:
        path = plan_dir / f"{action}_proteins.tsv"
        if not path.is_file():
            raise ValueError(f"Missing planner action table: {path}")
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                protein_id = row["protein_id"]
                if (
                    not protein_id
                    or protein_id in {".", ".."}
                    or Path(protein_id).name != protein_id
                    or "/" in protein_id
                    or "\\" in protein_id
                    or any(character.isspace() for character in protein_id)
                ):
                    raise ValueError(f"Unsafe protein ID for cache filename: {protein_id}")
                if row["action"] != action:
                    raise ValueError(f"Action mismatch for {protein_id}: {row['action']} != {action}")
                if protein_id in actions:
                    raise ValueError(f"Protein appears in both action tables: {protein_id}")
                actions[protein_id] = row
                by_action[action][protein_id] = row

    summary = json.loads((plan_dir / "summary.json").read_text(encoding="utf-8"))
    expected = int(summary["counts"]["target_proteins"])
    if len(actions) != expected:
        raise ValueError(f"Action partition has {len(actions)} proteins, expected {expected}")
    return by_action["reuse"], by_action["regenerate"]


def validate_array(path: Path, dimension: int) -> tuple[str, tuple[int, ...]]:
    try:
        array = np.load(path, allow_pickle=False)
    except Exception as exc:
        raise ValueError(f"Unreadable embedding {path}: {exc}") from exc
    if array.shape != (dimension,):
        raise ValueError(f"Unexpected embedding shape {path}: {array.shape} != ({dimension},)")
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError(f"Non-numeric embedding dtype for {path}: {array.dtype}")
    if not np.isfinite(array).all():
        raise ValueError(f"Non-finite values in embedding: {path}")
    return str(array.dtype), tuple(array.shape)


def link_or_copy(source: Path, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, destination)
        return "hardlink"
    except OSError:
        shutil.copy2(source, destination)
        return "copy"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--published-cache", type=Path, required=True)
    parser.add_argument("--generated-cache", type=Path, required=True)
    parser.add_argument("--output-cache", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    return parser.parse_args()


def validate_generated_scope(generated_cache: Path, regenerate: set[str]) -> None:
    for modality, specification in MODALITIES.items():
        directory = generated_cache / specification["directory"]
        if not directory.is_dir():
            continue
        unexpected = sorted(
            path.stem
            for path in directory.glob("*.npy")
            if path.stem not in regenerate
        )
        if unexpected:
            sample = ", ".join(unexpected[:5])
            raise ValueError(
                f"Generated {modality} cache contains {len(unexpected)} proteins "
                f"outside the regenerate partition; sample: {sample}"
            )


def main() -> int:
    args = parse_args()
    if args.output_cache.exists():
        raise SystemExit(f"Refusing to overwrite output cache: {args.output_cache}")
    args.output_cache.mkdir(parents=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    reuse, regenerate = load_actions(args.plan_dir)
    validate_generated_scope(args.generated_cache, set(regenerate))
    action_sets = {"reuse": reuse, "regenerate": regenerate}
    source_roots = {"reuse": args.published_cache, "regenerate": args.generated_cache}
    counts = Counter()
    details_path = args.report_dir / "embedding_assembly.tsv.gz"
    missing_handles = {}

    try:
        for action in action_sets:
            for modality in MODALITIES:
                path = args.report_dir / f"{action}_missing_{modality}.txt"
                missing_handles[(action, modality)] = path.open("w", encoding="utf-8")

        with gzip.open(details_path, "wt", encoding="utf-8", newline="") as detail_handle:
            fields = [
                "protein_id",
                "action",
                "modality",
                "status",
                "source",
                "destination",
                "dtype",
                "dimension",
                "transfer",
            ]
            writer = csv.DictWriter(detail_handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()

            for action, proteins in action_sets.items():
                for protein_id in sorted(proteins):
                    for modality, specification in MODALITIES.items():
                        directory = specification["directory"]
                        source = source_roots[action] / directory / f"{protein_id}.npy"
                        destination = args.output_cache / directory / f"{protein_id}.npy"
                        record = {
                            "protein_id": protein_id,
                            "action": action,
                            "modality": modality,
                            "source": str(source),
                            "destination": str(destination),
                            "dtype": "",
                            "dimension": specification["dimension"],
                            "transfer": "",
                        }
                        if not source.is_file():
                            record["status"] = "missing"
                            missing_handles[(action, modality)].write(protein_id + "\n")
                            counts[(action, modality, "missing")] += 1
                        else:
                            dtype, shape = validate_array(source, specification["dimension"])
                            record["status"] = "available"
                            record["dtype"] = dtype
                            record["dimension"] = shape[0]
                            record["transfer"] = link_or_copy(source, destination)
                            counts[(action, modality, "available")] += 1
                            counts[("transfer", record["transfer"])] += 1
                        writer.writerow(record)
    finally:
        for handle in missing_handles.values():
            handle.close()

    regenerate_count = len(regenerate)
    if counts[("regenerate", "prott5", "available")] != regenerate_count:
        raise SystemExit(
            "ProtT5 generation is incomplete: "
            f"{counts[('regenerate', 'prott5', 'available')]} / {regenerate_count}"
        )
    for modality in ("text", "structure", "ppi"):
        if regenerate_count and counts[("regenerate", modality, "available")] == 0:
            raise SystemExit(f"No regenerated {modality} embeddings were produced")

    summary = {
        "schema_version": 1,
        "target_proteins": len(reuse) + len(regenerate),
        "reuse_proteins": len(reuse),
        "regenerate_proteins": len(regenerate),
        "modalities": {},
        "transfer": {
            "hardlinks": counts[("transfer", "hardlink")],
            "copies": counts[("transfer", "copy")],
        },
        "policy": {
            "reuse": "Use an authenticated published array only for an exact planner-approved ID and sequence.",
            "regenerate": "Use only a newly generated array; never fall back to a published array.",
            "missing": "Leave absent so PFP applies its existing zero-vector and mask behavior.",
        },
    }
    for modality in MODALITIES:
        summary["modalities"][modality] = {}
        for action, proteins in action_sets.items():
            available = counts[(action, modality, "available")]
            missing = counts[(action, modality, "missing")]
            summary["modalities"][modality][action] = {
                "proteins": len(proteins),
                "available": available,
                "missing": missing,
                "coverage": available / len(proteins) if proteins else 0.0,
            }
        total_available = sum(
            counts[(action, modality, "available")] for action in action_sets
        )
        summary["modalities"][modality]["combined"] = {
            "available": total_available,
            "missing": summary["target_proteins"] - total_available,
            "coverage": total_available / summary["target_proteins"],
        }

    summary_path = args.report_dir / "assembly_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
