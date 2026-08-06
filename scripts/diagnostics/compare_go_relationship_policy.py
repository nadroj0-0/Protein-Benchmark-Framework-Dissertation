#!/usr/bin/env python3
"""Compare final PFP CSVs built under two GO relationship policies."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd


ASPECTS = ("bp", "cc", "mf")
SPLITS = ("training", "validation", "test")
ROOTS = {"bp": "GO:0008150", "cc": "GO:0005575", "mf": "GO:0003674"}
SCHEMA_NAME = "go-relationship-policy-final-csv-comparison"
SCHEMA_VERSION = 1


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _read_csv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype={"proteins": str, "sequences": str})
    if list(frame.columns[:2]) != ["proteins", "sequences"]:
        raise ValueError(f"Unexpected leading columns in {path}")
    if frame["proteins"].duplicated().any():
        raise ValueError(f"Duplicate proteins in {path}")
    frame = frame.set_index("proteins", verify_integrity=True)
    label_columns = list(frame.columns[1:])
    if not label_columns or any(not re.fullmatch(r"GO:\d{7}", item) for item in label_columns):
        raise ValueError(f"Invalid GO label columns in {path}")
    values = frame[label_columns].to_numpy(copy=False)
    if not np.isin(values, (0, 1)).all():
        raise ValueError(f"Non-binary label values in {path}")
    return frame


def _summarise_frame(frame: pd.DataFrame, aspect: str) -> dict[str, object]:
    labels = frame.iloc[:, 1:].to_numpy(dtype=np.uint8, copy=False)
    positive_counts = labels.sum(axis=1)
    root = ROOTS[aspect]
    root_index = list(frame.columns[1:]).index(root) if root in frame.columns else None
    root_only = 0
    if root_index is not None:
        root_only = int(np.count_nonzero((positive_counts == 1) & (labels[:, root_index] == 1)))
    return {
        "proteins": int(len(frame)),
        "terms": int(labels.shape[1]),
        "positive_labels": int(positive_counts.sum()),
        "root_only_proteins": root_only,
    }


def _compare_file(broad_path: Path, narrow_path: Path, aspect: str) -> dict[str, object]:
    broad = _read_csv(broad_path)
    narrow = _read_csv(narrow_path)
    broad_terms = list(broad.columns[1:])
    narrow_terms = list(narrow.columns[1:])
    common_terms = sorted(set(broad_terms) & set(narrow_terms))
    broad_only_terms = sorted(set(broad_terms) - set(narrow_terms))
    narrow_only_terms = sorted(set(narrow_terms) - set(broad_terms))
    common_proteins = broad.index.intersection(narrow.index, sort=False)
    broad_only_proteins = broad.index.difference(narrow.index, sort=False)
    narrow_only_proteins = narrow.index.difference(broad.index, sort=False)

    broad_common = broad.loc[common_proteins]
    narrow_common = narrow.loc[common_proteins]
    sequence_disagreements = int(np.count_nonzero(
        broad_common["sequences"].to_numpy() != narrow_common["sequences"].to_numpy()
    ))
    if sequence_disagreements:
        raise ValueError(
            f"Sequence content differs for {sequence_disagreements} common proteins in {aspect}"
        )

    broad_matrix = broad_common[common_terms].to_numpy(dtype=np.uint8, copy=False)
    narrow_matrix = narrow_common[common_terms].to_numpy(dtype=np.uint8, copy=False)
    changed = np.any(broad_matrix != narrow_matrix, axis=1)
    broad_only_positive = int(np.count_nonzero((broad_matrix == 1) & (narrow_matrix == 0)))
    narrow_only_positive = int(np.count_nonzero((narrow_matrix == 1) & (broad_matrix == 0)))
    if broad_only_terms:
        extra = broad_common[broad_only_terms].to_numpy(dtype=np.uint8, copy=False)
        broad_only_positive += int(extra.sum())
        changed |= np.any(extra == 1, axis=1)
    if narrow_only_terms:
        extra = narrow_common[narrow_only_terms].to_numpy(dtype=np.uint8, copy=False)
        narrow_only_positive += int(extra.sum())
        changed |= np.any(extra == 1, axis=1)

    result = {
        "broad": _summarise_frame(broad, aspect),
        "narrow": _summarise_frame(narrow, aspect),
        "common_proteins": int(len(common_proteins)),
        "broad_only_proteins": int(len(broad_only_proteins)),
        "narrow_only_proteins": int(len(narrow_only_proteins)),
        "common_proteins_with_changed_labels": int(np.count_nonzero(changed)),
        "broad_only_positive_labels_on_common_proteins": broad_only_positive,
        "narrow_only_positive_labels_on_common_proteins": narrow_only_positive,
        "common_terms": len(common_terms),
        "broad_only_terms": broad_only_terms,
        "narrow_only_terms": narrow_only_terms,
        "files": {
            "broad": {"path": str(broad_path), "sha256": _sha256(broad_path)},
            "narrow": {"path": str(narrow_path), "sha256": _sha256(narrow_path)},
        },
    }
    del broad, narrow, broad_common, narrow_common, broad_matrix, narrow_matrix
    return result


def _obo_relationship_counts(path: Path) -> dict[str, object]:
    namespaces: dict[str, str] = {}
    relationships: list[tuple[str, str, str]] = []
    current = ""
    namespace = ""
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line == "[Term]":
                if current:
                    namespaces[current] = namespace
                current, namespace = "", ""
            elif line == "[Typedef]":
                if current:
                    namespaces[current] = namespace
                current, namespace = "", ""
            elif line.startswith("id: GO:"):
                current = line.split("id: ", 1)[1]
            elif line.startswith("namespace: "):
                namespace = line.split(": ", 1)[1]
            elif current and line.startswith("relationship: "):
                fields = line.split()
                if len(fields) >= 3:
                    relationships.append((current, fields[1], fields[2]))
    if current:
        namespaces[current] = namespace
    counts: dict[str, dict[str, int]] = {}
    for child, relation, parent in relationships:
        item = counts.setdefault(relation, {"edges": 0, "cross_namespace_edges": 0})
        item["edges"] += 1
        item["cross_namespace_edges"] += int(
            bool(namespaces.get(child))
            and bool(namespaces.get(parent))
            and namespaces[child] != namespaces[parent]
        )
    return {
        "path": str(path),
        "sha256": _sha256(path),
        "relationship_types": dict(sorted(counts.items())),
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    broad_dir = args.broad_dir.resolve()
    narrow_dir = args.narrow_dir.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(output)
    results: dict[str, object] = {}
    for aspect in ASPECTS:
        results[aspect] = {}
        for split_name in SPLITS:
            name = f"{aspect}-{split_name}.csv"
            results[aspect][split_name] = _compare_file(
                broad_dir / name, narrow_dir / name, aspect
            )
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": _utc_now(),
        "broad": {"label": args.broad_label, "path": str(broad_dir)},
        "narrow": {"label": args.narrow_label, "path": str(narrow_dir)},
        "policy_boundary": {
            "broad": "is_a plus every relationship line",
            "narrow": "is_a plus part_of only",
            "interpretation": "Exact final model-input difference; no model was retrained.",
        },
        "ontology": _obo_relationship_counts(args.obo_file.resolve()),
        "comparisons": results,
    }
    output.mkdir(parents=True)
    audit_path = output / "audit.json"
    audit_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (output / "split_summary.tsv").open("w", encoding="utf-8") as handle:
        handle.write(
            "aspect\tsplit\tbroad_proteins\tnarrow_proteins\tcommon_proteins\t"
            "broad_only_proteins\tnarrow_only_proteins\tchanged_common_proteins\t"
            "broad_positive_labels\tnarrow_positive_labels\tbroad_only_common_labels\t"
            "narrow_only_common_labels\tbroad_root_only\tnarrow_root_only\n"
        )
        for aspect in ASPECTS:
            for split_name in SPLITS:
                row = results[aspect][split_name]
                handle.write(
                    f"{aspect}\t{split_name}\t{row['broad']['proteins']}\t"
                    f"{row['narrow']['proteins']}\t{row['common_proteins']}\t"
                    f"{row['broad_only_proteins']}\t{row['narrow_only_proteins']}\t"
                    f"{row['common_proteins_with_changed_labels']}\t"
                    f"{row['broad']['positive_labels']}\t{row['narrow']['positive_labels']}\t"
                    f"{row['broad_only_positive_labels_on_common_proteins']}\t"
                    f"{row['narrow_only_positive_labels_on_common_proteins']}\t"
                    f"{row['broad']['root_only_proteins']}\t{row['narrow']['root_only_proteins']}\n"
                )
    lines = [
        "# GO relationship-policy audit", "", "## Scope", "",
        "This read-only audit compares the accepted all-relationship benchmark with an exact diagnostic rebuild using only `is_a + part_of`. It measures final CSV changes and does not retrain a model or replace an accepted result.",
        "", "| Aspect | Split | Broad rows | Narrow rows | Changed common rows | Broad-only positive labels | Broad root-only | Narrow root-only |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for aspect in ASPECTS:
        for split_name in SPLITS:
            row = results[aspect][split_name]
            lines.append(
                f"| {aspect.upper()} | {split_name} | {row['broad']['proteins']:,} | "
                f"{row['narrow']['proteins']:,} | "
                f"{row['common_proteins_with_changed_labels']:,} | "
                f"{row['broad_only_positive_labels_on_common_proteins']:,} | "
                f"{row['broad']['root_only_proteins']:,} | {row['narrow']['root_only_proteins']:,} |"
            )
    lines.extend(["", "## Interpretation boundary", "", "- These are exact construction differences, not performance differences.", "- Any material BP change would motivate a separately labelled sensitivity retrain; it does not invalidate or silently replace the accepted compatibility baseline."])
    summary_path = output / "summary.md"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest_files = (audit_path, output / "split_summary.tsv", summary_path)
    manifest = {path.name: {"sha256": _sha256(path), "size_bytes": path.stat().st_size} for path in manifest_files}
    (output / "output_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "RUN_COMPLETE.json").write_text(json.dumps({"complete": True, "schema_name": SCHEMA_NAME, "schema_version": SCHEMA_VERSION, "audit_sha256": _sha256(audit_path)}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--broad-label", required=True)
    parser.add_argument("--broad-dir", type=Path, required=True)
    parser.add_argument("--narrow-label", required=True)
    parser.add_argument("--narrow-dir", type=Path, required=True)
    parser.add_argument("--obo-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
