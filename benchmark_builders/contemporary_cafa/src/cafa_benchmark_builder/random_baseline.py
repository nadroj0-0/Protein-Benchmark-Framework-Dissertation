from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import uuid
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from .builder import propagate_annotations
from .config import PREFIX_TO_NAMESPACE, SUPERVISOR_EXP_CODES
from .goa import load_normalized_annotation_map
from .ontology import Ontology
from .parsers import iter_uniprot


SOURCE_PICKLES = {
    "training": "train_data_train.pkl",
    "validation": "train_data_valid.pkl",
    "test": "test_data.pkl",
}
ASPECTS = ("bp", "cc", "mf")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_source(
    source_dir: Path,
) -> tuple[pd.DataFrame, dict[str, int], list[str], dict[str, dict[str, object]]]:
    frames = []
    inputs = {}
    for split, filename in SOURCE_PICKLES.items():
        path = source_dir / filename
        if not path.is_file():
            raise FileNotFoundError(f"Missing source benchmark pickle: {path}")
        frame = pd.read_pickle(path)
        required = {"proteins", "sequences", "annotations"}
        if not required.issubset(frame.columns):
            raise ValueError(f"{path} lacks required columns: {sorted(required - set(frame.columns))}")
        frame = frame.loc[:, ["proteins", "sequences", "annotations"]].copy()
        frame["annotations"] = frame["annotations"].map(lambda value: tuple(sorted(set(value))))
        frames.append(frame)
        inputs[filename] = {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }

    combined = pd.concat(frames, ignore_index=True)
    duplicates = combined[combined.proteins.duplicated(keep=False)]
    if not duplicates.empty:
        raise ValueError(
            "Source benchmark is not globally protein-disjoint: "
            + ", ".join(sorted(duplicates.proteins.astype(str).unique())[:10])
        )
    if combined.sequences.eq("").any():
        raise ValueError("Source benchmark contains an empty protein sequence")

    split_members = {}
    for split in SOURCE_PICKLES:
        members = set()
        for prefix in ASPECTS:
            path = source_dir / f"{prefix}-{split}.csv"
            if not path.is_file():
                raise FileNotFoundError(f"Missing source benchmark CSV: {path}")
            members.update(pd.read_csv(path, usecols=["proteins"]).proteins.astype(str))
        split_members[split] = members
    for left, right in (("training", "validation"), ("training", "test"), ("validation", "test")):
        overlap = split_members[left] & split_members[right]
        if overlap:
            raise ValueError(
                f"Source benchmark is not globally split-disjoint ({left}/{right}): "
                + ", ".join(sorted(overlap)[:10])
            )
    visible = set().union(*split_members.values())
    missing = sorted(visible - set(combined.proteins.astype(str)))
    if missing:
        raise ValueError("Source CSV proteins are missing from source pickles: " + ", ".join(missing[:10]))
    combined = combined[combined.proteins.astype(str).isin(visible)].copy()

    terms_path = source_dir / "terms.pkl"
    if not terms_path.is_file():
        raise FileNotFoundError(f"Missing source term universe: {terms_path}")
    terms = pd.read_pickle(terms_path).terms.astype(str).tolist()
    if not terms or len(terms) != len(set(terms)):
        raise ValueError("Source term universe is empty or contains duplicate GO terms")
    inputs["terms.pkl"] = {
        "path": str(terms_path.resolve()),
        "size_bytes": terms_path.stat().st_size,
        "sha256": _sha256(terms_path),
    }
    counts = {split: len(members) for split, members in split_members.items()}
    return (
        combined.sort_values("proteins", kind="stable").reset_index(drop=True),
        counts,
        terms,
        inputs,
    )


def _relabel_at_snapshot(
    source: pd.DataFrame,
    uniprot_paths: tuple[Path, ...],
    goa_path: Path,
    source_go: Ontology,
    benchmark_go: Ontology,
    evidence_codes: frozenset[str],
) -> tuple[pd.DataFrame, dict[str, object]]:
    source_ids = set(source.proteins.astype(str))
    source_sequences = dict(zip(source.proteins.astype(str), source.sequences.astype(str)))
    records = {}
    source_candidates: dict[str, set[str]] = defaultdict(set)
    alias_candidates: dict[str, set[str]] = defaultdict(set)

    for path in uniprot_paths:
        for record in iter_uniprot(path):
            identifiers = set(record.accessions) | {record.protein_id}
            matched = identifiers & source_ids
            if not matched:
                continue
            existing = records.get(record.protein_id)
            if existing is not None and existing.sequence != record.sequence:
                raise ValueError(f"Conflicting t1 sequences for {record.protein_id}")
            records[record.protein_id] = record
            for source_id in matched:
                source_candidates[source_id].add(record.protein_id)
            for alias in identifiers:
                alias_candidates[alias].add(record.protein_id)

    missing_sequences = sorted(source_ids - set(source_candidates))

    # Keep the paired source population fixed when an accession has left UniProt.
    source_to_primary = {source_id: source_id for source_id in missing_sequences}
    sequence_disambiguated = 0
    for source_id, candidates in source_candidates.items():
        if len(candidates) == 1:
            source_to_primary[source_id] = next(iter(candidates))
            continue
        exact = {
            candidate
            for candidate in candidates
            if records[candidate].sequence == source_sequences[source_id]
        }
        if len(exact) != 1:
            raise ValueError(
                f"Ambiguous t1 mapping for source protein {source_id}: "
                + ", ".join(sorted(candidates))
            )
        source_to_primary[source_id] = next(iter(exact))
        sequence_disambiguated += 1

    primary_source_counts = Counter(source_to_primary.values())
    alias_to_primary = {
        alias: next(iter(primaries))
        for alias, primaries in alias_candidates.items()
        if len(primaries) == 1
    }
    alias_to_primary.update({source_id: source_id for source_id in missing_sequences})
    annotations = load_normalized_annotation_map(
        goa_path,
        alias_to_primary=alias_to_primary,
        source_ontology=source_go,
        benchmark_ontology=benchmark_go,
        evidence_codes=evidence_codes,
        snapshot="random-t1",
    )

    rows = []
    missing_labels = []
    for source_id in sorted(source_ids):
        primary = source_to_primary[source_id]
        direct = annotations.annotations.get(primary)
        propagated = propagate_annotations(benchmark_go, direct or set())
        if not propagated:
            missing_labels.append(source_id)
            continue
        rows.append({
            "proteins": source_id,
            "sequences": source_sequences[source_id],
            "annotations": tuple(sorted(propagated)),
        })
    return pd.DataFrame(rows), {
        "source_population_proteins": len(source_to_primary),
        "mapped_source_proteins": len(source_to_primary) - len(missing_sequences),
        "sequence_disambiguated_source_proteins": sequence_disambiguated,
        "source_sequence_fallback_proteins": len(missing_sequences),
        "merged_t1_primary_groups": sum(
            count > 1 for count in primary_source_counts.values()
        ),
        "source_proteins_in_merged_t1_primary_groups": sum(
            count for count in primary_source_counts.values() if count > 1
        ),
        "qualifying_t1_proteins": len(rows),
        "nonqualifying_t1_proteins_excluded": len(missing_labels),
        "nonqualifying_t1_protein_sample": missing_labels[:10],
        "goa_counters": dict(sorted(annotations.counters.items())),
        "goa_evidence_counts": dict(sorted(annotations.evidence_counts.items())),
    }


def _take_groups(groups: list[list[int]], target: int) -> tuple[list[list[int]], list[list[int]]]:
    selected = []
    remaining = []
    selected_rows = 0
    for group in groups:
        if selected_rows < target and len(group) <= target - selected_rows:
            selected.append(group)
            selected_rows += len(group)
        else:
            remaining.append(group)
    if selected_rows != target:
        raise ValueError(
            f"Exact-sequence groups cannot satisfy requested split size {target}; reached {selected_rows}"
        )
    return selected, remaining


def _random_split(
    source: pd.DataFrame,
    source_counts: dict[str, int],
    seed: int,
) -> dict[str, pd.DataFrame]:
    sequence_groups: dict[str, list[int]] = defaultdict(list)
    for index, sequence in enumerate(source.sequences):
        sequence_groups[str(sequence)].append(index)
    groups = list(sequence_groups.values())
    np.random.RandomState(seed).shuffle(groups)
    test_groups, groups = _take_groups(groups, source_counts["test"])
    validation_groups, training_groups = _take_groups(groups, source_counts["validation"])
    grouped = {
        "training": training_groups,
        "validation": validation_groups,
        "test": test_groups,
    }
    frames = {}
    for split, split_groups in grouped.items():
        indices = [index for group in split_groups for index in group]
        frames[split] = source.iloc[indices].sort_values("proteins", kind="stable").reset_index(drop=True)
    return frames


def _write_outputs(
    output_dir: Path,
    frames: dict[str, pd.DataFrame],
    terms: list[str],
    go: Ontology,
) -> tuple[dict[str, dict[str, int]], dict[str, int]]:
    development = pd.concat([frames["training"], frames["validation"]], ignore_index=True)
    development.sort_values("proteins", kind="stable").to_pickle(output_dir / "train_data.pkl")
    frames["training"].to_pickle(output_dir / "train_data_train.pkl")
    frames["validation"].to_pickle(output_dir / "train_data_valid.pkl")
    frames["test"].to_pickle(output_dir / "test_data.pkl")
    pd.DataFrame({"terms": terms}).to_pickle(output_dir / "terms.pkl")

    stats = {}
    visible_by_split: dict[str, set[str]] = defaultdict(set)
    for prefix, namespace in PREFIX_TO_NAMESPACE.items():
        aspect_terms = [term for term in terms if go.get_namespace(term) == namespace]
        for split, frame in frames.items():
            path = output_dir / f"{prefix}-{split}.csv"
            rows = 0
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(["proteins", "sequences", *aspect_terms])
                for row in frame.itertuples(index=False):
                    labels = set(row.annotations)
                    if not labels.intersection(aspect_terms):
                        continue
                    writer.writerow([
                        row.proteins,
                        row.sequences,
                        *(int(term in labels) for term in aspect_terms),
                    ])
                    rows += 1
                    visible_by_split[split].add(str(row.proteins))
            if rows == 0:
                raise ValueError(f"Random baseline produced an empty split: {path.name}")
            stats[path.name] = {"rows": rows, "terms": len(aspect_terms)}
    return stats, {split: len(visible_by_split[split]) for split in frames}


def _output_manifest(directory: Path) -> dict[str, object]:
    files = []
    for path in sorted(directory.iterdir()):
        if not path.is_file() or path.name in {"output_manifest.json", "RUN_COMPLETE.json"}:
            continue
        files.append({
            "path": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        })
    return {"schema_version": 1, "payload_file_count": len(files), "files": files}


def build_random_baseline(
    *,
    mode: str,
    benchmark_id: str,
    source_dir: Path,
    go_obo: Path,
    output_dir: Path,
    seed: int = 0,
    evidence_codes: frozenset[str] = SUPERVISOR_EXP_CODES,
    uniprot_paths: tuple[Path, ...] = (),
    goa_path: Path | None = None,
    source_go_obo: Path | None = None,
    input_manifest: Path | None = None,
) -> Path:
    if mode not in {"prepared", "t1-snapshot"}:
        raise ValueError(f"Unsupported random baseline mode: {mode}")
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")

    source, source_counts, terms, source_inputs = _load_source(source_dir)
    go = Ontology(go_obo, with_rels=True)
    snapshot_stats = None
    if mode == "t1-snapshot":
        if not uniprot_paths or goa_path is None:
            raise ValueError("t1-snapshot mode requires UniProt inputs and a GOA file")
        source_go = Ontology(source_go_obo or go_obo, with_rels=True)
        source, snapshot_stats = _relabel_at_snapshot(
            source, uniprot_paths, goa_path, source_go, go, evidence_codes
        )

    requested_counts = dict(source_counts)
    if mode == "t1-snapshot":
        requested_counts["training"] = (
            len(source) - source_counts["validation"] - source_counts["test"]
        )
        if requested_counts["training"] <= 0:
            raise ValueError(
                "The qualifying t1 population is too small to preserve the paired "
                "validation and test counts"
            )

    frames = _random_split(source, requested_counts, seed)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = output_dir.parent / f".{output_dir.name}.part-{uuid.uuid4().hex}"
    stage.mkdir()
    try:
        csv_stats, csv_union_counts = _write_outputs(stage, frames, terms, go)
        actual_counts = {split: len(frame) for split, frame in frames.items()}
        checks = {
            "selected_population_preserved": sum(actual_counts.values()) == len(source),
            "requested_split_sizes_met": actual_counts == requested_counts,
            "paired_validation_size_preserved": (
                actual_counts["validation"] == source_counts["validation"]
            ),
            "paired_test_size_preserved": actual_counts["test"] == source_counts["test"],
            "selected_population_remains_model_visible": csv_union_counts == actual_counts,
            "global_protein_disjoint": all(
                set(frames[left].proteins).isdisjoint(frames[right].proteins)
                for left, right in (("training", "validation"), ("training", "test"), ("validation", "test"))
            ),
            "global_exact_sequence_disjoint": all(
                set(frames[left].sequences).isdisjoint(frames[right].sequences)
                for left, right in (("training", "validation"), ("training", "test"), ("validation", "test"))
            ),
        }
        if not all(checks.values()):
            raise ValueError(f"Random baseline validation failed: {checks}")

        manifest = {
            "schema_version": 1,
            "benchmark_id": benchmark_id,
            "benchmark_design": "random-exact-sequence-grouped",
            "label_source": "inherited-prepared-labels" if mode == "prepared" else "t1-snapshot-membership",
            "sequence_source_policy": "frozen-from-source-paired-control",
            "mode": mode,
            "seed": seed,
            "evidence_codes": sorted(evidence_codes),
            "relationship_policy": "all-relationships",
            "term_universe_policy": "frozen-from-source-paired-control",
            "source_benchmark_dir": str(source_dir.resolve()),
            "source_split_counts": source_counts,
            "requested_output_split_counts": requested_counts,
            "output_split_counts": actual_counts,
            "output_csv_union_counts": csv_union_counts,
            "population_policy": (
                "t1-qualifying-only; preserve source validation and test counts; "
                "exclude nonqualifying proteins from training capacity"
                if mode == "t1-snapshot"
                else "preserve accepted source population and split counts"
            ),
            "source_files": source_inputs,
            "go_obo": {"path": str(go_obo.resolve()), "sha256": _sha256(go_obo)},
            "snapshot": snapshot_stats,
            "snapshot_inputs": (
                {
                    "uniprot": [str(path.resolve()) for path in uniprot_paths],
                    "goa": str(goa_path.resolve()) if goa_path else None,
                    "source_go_obo": str((source_go_obo or go_obo).resolve()),
                }
                if mode == "t1-snapshot" else None
            ),
            "frozen_input_manifest": (
                {"path": str(input_manifest.resolve()), "sha256": _sha256(input_manifest)}
                if input_manifest else None
            ),
            "csv_outputs": csv_stats,
        }
        (stage / "build_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        validation = {"valid": True, "checks": checks}
        (stage / "validation_report.json").write_text(
            json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        scientific_fingerprint = hashlib.sha256(
            json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        output_manifest = _output_manifest(stage)
        (stage / "output_manifest.json").write_text(
            json.dumps(output_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        completion = {
            "complete": True,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "manifest": "output_manifest.json",
            "manifest_sha256": _sha256(stage / "output_manifest.json"),
            "scientific_fingerprint": scientific_fingerprint,
        }
        (stage / "RUN_COMPLETE.json").write_text(
            json.dumps(completion, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(stage, output_dir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a random split paired control from an accepted PFP benchmark."
    )
    parser.add_argument("--mode", choices=("prepared", "t1-snapshot"), required=True)
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument("--source-benchmark-dir", type=Path, required=True)
    parser.add_argument("--go-obo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--evidence-code", action="append", default=[])
    parser.add_argument("--uniprot", type=Path, action="append", default=[])
    parser.add_argument("--goa", type=Path)
    parser.add_argument("--source-go-obo", type=Path)
    parser.add_argument("--input-manifest", type=Path)
    args = parser.parse_args()
    build_random_baseline(
        mode=args.mode,
        benchmark_id=args.benchmark_id,
        source_dir=args.source_benchmark_dir,
        go_obo=args.go_obo,
        output_dir=args.output_dir,
        seed=args.seed,
        evidence_codes=(
            frozenset(args.evidence_code) if args.evidence_code else SUPERVISOR_EXP_CODES
        ),
        uniprot_paths=tuple(args.uniprot),
        goa_path=args.goa,
        source_go_obo=args.source_go_obo,
        input_manifest=args.input_manifest,
    )


if __name__ == "__main__":
    main()
