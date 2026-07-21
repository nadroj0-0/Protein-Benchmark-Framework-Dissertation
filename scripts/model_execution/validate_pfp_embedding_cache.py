#!/usr/bin/env python3
"""Exhaustively validate the exact embedding cache PFP will consume."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict

import numpy as np

from common import (
    ASPECTS,
    MODALITY_MODES,
    PFP_SPLITS,
    active_modalities,
    atomic_write_json,
    file_snapshot,
    load_run_config,
    modality_paths,
    require_unchanged,
    selected_aspects,
    sha256_file,
)


def load_memberships(
    data_dir: Path, aspects: list[str]
) -> tuple[set[str], Dict[tuple[str, str], list[str]]]:
    union: set[str] = set()
    memberships: Dict[tuple[str, str], list[str]] = {}
    for aspect in aspects:
        for split in PFP_SPLITS.values():
            path = data_dir / f"{aspect}_{split}_names.npy"
            if not path.is_file():
                raise FileNotFoundError(f"Prepared split IDs are missing: {path}")
            values = [str(value) for value in np.load(path, allow_pickle=True)]
            memberships[(aspect, split)] = values
            union.update(values)
    if not union:
        raise ValueError("Prepared benchmark contains no protein IDs")
    return union, memberships


def inspect_array(path: Path, dimension: int) -> tuple[bool, str, str]:
    try:
        value = np.load(path, allow_pickle=False)
    except Exception as exc:
        return False, "unreadable", f"{type(exc).__name__}: {exc}"
    if value.shape != (dimension,):
        return False, "wrong_shape", str(tuple(value.shape))
    if value.dtype.kind != "f":
        return False, "unsupported_dtype", str(value.dtype)
    if not np.isfinite(value).all():
        return False, "non_finite", str(value.dtype)
    with np.errstate(over="ignore", invalid="ignore"):
        converted = value.astype(np.float32)
    if not np.isfinite(converted).all():
        return False, "non_finite_after_float32", str(value.dtype)
    return True, "valid", str(value.dtype)


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def prepared_target_hashes(data_dir: Path) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    for aspect in ASPECTS:
        for split in PFP_SPLITS.values():
            path = data_dir / f"{aspect}_{split}_sequences.json"
            payload = json.loads(path.read_text(encoding="utf-8"))
            for protein_id, sequence in payload.items():
                previous = sequences.get(str(protein_id))
                if previous is not None and previous != str(sequence):
                    raise ValueError(f"Conflicting prepared sequences for {protein_id}")
                sequences[str(protein_id)] = str(sequence)
    return {
        protein_id: hashlib.sha256(sequence.encode("utf-8")).hexdigest()
        for protein_id, sequence in sorted(sequences.items())
    }


def read_target_manifest(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["protein_id", "sequence_sha256"]:
            raise ValueError("targets.tsv has an unexpected schema")
        for row in reader:
            protein_id = row["protein_id"]
            if protein_id in result:
                raise ValueError(f"targets.tsv repeats protein {protein_id}")
            result[protein_id] = row["sequence_sha256"]
    return result


def read_pair_status(
    path: Path,
    targets: Dict[str, str],
    modalities: set[str],
) -> tuple[Dict[str, set[str]], Dict[str, Dict[str, str]]]:
    accepted: Dict[str, set[str]] = {modality: set() for modality in modalities}
    accepted_hashes: Dict[str, Dict[str, str]] = {
        modality: {} for modality in modalities
    }
    seen: set[tuple[str, str]] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "protein_id",
            "modality",
            "state",
            "sequence_sha256",
            "embedding_sha256",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("pair_status.tsv has an unexpected schema")
        for row in reader:
            protein_id = row["protein_id"]
            modality = row["modality"]
            key = (protein_id, modality)
            if protein_id not in targets or modality not in modalities:
                raise ValueError(f"pair_status.tsv contains unexpected pair {key}")
            if key in seen:
                raise ValueError(f"pair_status.tsv repeats pair {key}")
            seen.add(key)
            if row["sequence_sha256"] != targets[protein_id]:
                raise ValueError(f"pair_status.tsv sequence hash differs for {protein_id}")
            if row["state"] not in {"accepted", "needs_retry"}:
                raise ValueError(f"pair_status.tsv has invalid status for {key}")
            if row["state"] == "accepted":
                digest = row["embedding_sha256"]
                if len(digest) != 64 or any(
                    character not in "0123456789abcdef" for character in digest
                ):
                    raise ValueError(f"pair_status.tsv lacks a valid embedding hash for {key}")
                accepted[modality].add(protein_id)
                accepted_hashes[modality][protein_id] = digest
            elif row["embedding_sha256"]:
                raise ValueError(f"needs_retry pair unexpectedly has an embedding hash: {key}")
    expected_pairs = {(protein_id, modality) for protein_id in targets for modality in modalities}
    if seen != expected_pairs:
        raise ValueError("pair_status.tsv does not contain every target/modality pair exactly once")
    return accepted, accepted_hashes


def validate_ia_inputs(
    data_dir: Path,
    ia_dir: Path | None,
    aspects: list[str],
    required: bool,
) -> tuple[Dict[str, Any], list[str]]:
    result: Dict[str, Any] = {}
    failures: list[str] = []
    if ia_dir is None:
        if required:
            failures.append("This run config requires a precomputed IA directory")
        return result, failures
    for aspect in aspects:
        path = ia_dir / f"{aspect}_ia.txt"
        if not path.is_file():
            failures.append(f"Precomputed IA file is missing for {aspect}: {path}")
            continue
        values: Dict[str, float] = {}
        try:
            with path.open(encoding="utf-8", newline="") as handle:
                for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
                    if len(row) != 2 or not row[0]:
                        raise ValueError(f"invalid row {line_number}")
                    if row[0] in values:
                        raise ValueError(f"duplicate term {row[0]}")
                    value = float(row[1])
                    if not math.isfinite(value) or value < 0:
                        raise ValueError(f"invalid IA value for {row[0]}")
                    values[row[0]] = value
        except (OSError, ValueError) as exc:
            failures.append(f"Cannot validate IA file for {aspect}: {exc}")
            continue
        terms = json.loads((data_dir / f"{aspect}_go_terms.json").read_text(encoding="utf-8"))
        if set(values) != set(terms):
            failures.append(
                f"Precomputed IA term set differs from prepared GO terms for {aspect}"
            )
            continue
        result[aspect] = {
            "path": str(path.resolve()),
            "sha256": sha256_file(path),
            "terms": len(values),
        }
    return result, failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--mode", choices=MODALITY_MODES, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--issues-tsv", type=Path, required=True)
    parser.add_argument("--embedding-evidence", type=Path, action="append", default=[])
    parser.add_argument("--require-embedding-evidence", action="store_true")
    parser.add_argument("--preparation-report", type=Path, required=True)
    parser.add_argument("--ia-file-dir", type=Path)
    parser.add_argument("--aspect", action="append", default=[])
    args = parser.parse_args()

    config = load_run_config(args.config)
    cache_root = args.cache_root.resolve()
    directories = modality_paths(cache_root, config)
    aspects = selected_aspects(args.aspect)
    targets, memberships = load_memberships(args.data_dir, aspects)
    active = list(active_modalities(args.mode))
    summaries: Dict[str, Any] = {}
    valid_ids: Dict[str, set[str]] = defaultdict(set)
    valid_hashes: Dict[str, Dict[str, str]] = defaultdict(dict)
    issues: list[tuple[str, str, str, str, str]] = []

    for modality in active:
        specification = config["modalities"][modality]
        dimension = int(specification["dimension"])
        directory = directories[modality]
        if not directory.is_dir():
            raise FileNotFoundError(f"Embedding directory is missing for {modality}: {directory}")
        reasons: Dict[str, int] = defaultdict(int)
        dtypes: Dict[str, int] = defaultdict(int)
        missing = 0
        content_digest = hashlib.sha256()
        for protein_id in sorted(targets):
            path = directory / f"{protein_id}.npy"
            if not path.is_file():
                missing += 1
                reasons["missing"] += 1
                issues.append((protein_id, modality, "missing", "", str(path)))
                continue
            valid, reason, detail = inspect_array(path, dimension)
            reasons[reason] += 1
            if valid:
                valid_ids[modality].add(protein_id)
                dtypes[detail] += 1
                file_digest = sha256_file(path)
                valid_hashes[modality][protein_id] = file_digest
                content_digest.update(protein_id.encode("utf-8"))
                content_digest.update(b"\t")
                content_digest.update(file_digest.encode("ascii"))
                content_digest.update(b"\n")
            else:
                issues.append((protein_id, modality, reason, detail, str(path)))
        valid_count = len(valid_ids[modality])
        summaries[modality] = {
            "directory": str(directory),
            "dimension": dimension,
            "target_count": len(targets),
            "valid": valid_count,
            "missing": missing,
            "invalid_present": sum(
                count for reason, count in reasons.items() if reason not in {"valid", "missing"}
            ),
            "coverage": valid_count / len(targets),
            "reasons": dict(sorted(reasons.items())),
            "dtypes": dict(sorted(dtypes.items())),
            "valid_content_sha256": content_digest.hexdigest(),
        }

    failures: list[str] = []
    preparation_path = args.preparation_report.resolve()
    preparation_snapshot = file_snapshot(preparation_path)
    preparation = json.loads(preparation_path.read_text(encoding="utf-8"))
    require_unchanged(preparation_path, preparation_snapshot, "Preparation report")
    if preparation.get("status") != "passed":
        failures.append("Preparation report does not declare status=passed")
    evidence = []
    for path in args.embedding_evidence:
        resolved = path.resolve()
        if not resolved.is_file():
            failures.append(f"Embedding evidence is missing: {resolved}")
            continue
        evidence.append(
            {"path": str(resolved), "name": resolved.name, "sha256": sha256_file(resolved)}
        )
    evidence_by_name: Dict[str, Dict[str, str]] = {}
    for item in evidence:
        if item["name"] in evidence_by_name:
            failures.append(f"Repeated embedding evidence filename: {item['name']}")
        else:
            evidence_by_name[item["name"]] = item
    evidence_bound = False
    binding: Dict[str, Any] = {}
    required_state_files = {
        "coverage.json", "contract.json", "targets.tsv", "pair_status.tsv"
    }
    if required_state_files.issubset(evidence_by_name):
        try:
            coverage_payload = json.loads(
                Path(evidence_by_name["coverage.json"]["path"]).read_text(encoding="utf-8")
            )
            contract_payload = json.loads(
                Path(evidence_by_name["contract.json"]["path"]).read_text(encoding="utf-8")
            )
            target_path = Path(evidence_by_name["targets.tsv"]["path"])
            pair_status_path = Path(evidence_by_name["pair_status.tsv"]["path"])
            target_payload = read_target_manifest(target_path)
            prepared_targets = prepared_target_hashes(args.data_dir)
            if target_payload != prepared_targets:
                raise ValueError("targets.tsv does not match prepared protein IDs and sequences")
            contract_copy = dict(contract_payload)
            recorded_contract_hash = contract_copy.pop("contract_sha256", None)
            if recorded_contract_hash != canonical_sha256(contract_copy):
                raise ValueError("contract.json self-hash is invalid")
            if coverage_payload.get("contract_sha256") != recorded_contract_hash:
                raise ValueError("coverage.json is not bound to contract.json")
            contract_targets = contract_payload.get("targets", {})
            if int(contract_targets.get("count", -1)) != len(prepared_targets):
                raise ValueError("contract.json target count differs from prepared benchmark")
            if contract_targets.get("manifest_sha256") != sha256_file(target_path):
                raise ValueError("contract.json is not bound to targets.tsv")
            source_hashes = preparation.get("source_csv_sha256", {})
            source_sizes = preparation.get("source_csv_bytes", {})
            contract_csvs = contract_payload.get("benchmark_csvs", [])
            observed_hashes = {
                item.get("name"): item.get("sha256")
                for item in contract_csvs
                if isinstance(item, dict)
            }
            if observed_hashes != source_hashes:
                raise ValueError("contract.json is not bound to the selected nine CSVs")
            observed_sizes = {
                item.get("name"): item.get("size_bytes")
                for item in contract_csvs
                if isinstance(item, dict)
            }
            if observed_sizes != source_sizes:
                raise ValueError("contract.json CSV sizes differ from the selected nine CSVs")
            if contract_payload.get("pfp_commit") != "1e04fd6d6d3c40458fd41ec1a881ed6e24de768e":
                raise ValueError("contract.json records an unexpected PFP commit")
            state_modalities = contract_payload.get("policy", {}).get("modalities", {})
            for modality, specification in config["modalities"].items():
                state_specification = state_modalities.get(modality, {})
                if int(state_specification.get("dimension", -1)) != int(specification["dimension"]):
                    raise ValueError(f"Embedding contract dimension differs for {modality}")
                expected_directory = Path(str(specification["directory"])).name
                if state_specification.get("cache_directory") != expected_directory:
                    raise ValueError(f"Embedding contract cache directory differs for {modality}")
            if int(coverage_payload.get("target_count", -1)) != len(prepared_targets):
                raise ValueError("coverage.json target count differs from prepared benchmark")
            coverage_values = coverage_payload.get("coverage", {})
            accepted_ids, accepted_hashes = read_pair_status(
                pair_status_path, prepared_targets, set(config["modalities"])
            )
            count_mismatches = [
                modality
                for modality in config["modalities"]
                if int(coverage_values.get(modality, {}).get("accepted", -1))
                != len(accepted_ids[modality])
            ]
            if count_mismatches:
                raise ValueError(
                    "coverage.json differs from pair_status.tsv for: "
                    + ", ".join(count_mismatches)
                )
            id_mismatches = [
                modality
                for modality in active
                if (accepted_ids[modality] & targets) != valid_ids[modality]
            ]
            if id_mismatches:
                raise ValueError(
                    "pair_status.tsv accepted IDs differ from validated cache for: "
                    + ", ".join(id_mismatches)
                )
            hash_mismatches = [
                modality
                for modality in active
                if {
                    protein_id: digest
                    for protein_id, digest in accepted_hashes[modality].items()
                    if protein_id in targets
                }
                != valid_hashes[modality]
            ]
            if hash_mismatches:
                raise ValueError(
                    "pair_status.tsv embedding hashes differ from validated cache for: "
                    + ", ".join(hash_mismatches)
                )
            evidence_bound = True
            binding = {
                "passed": True,
                "contract_sha256": recorded_contract_hash,
                "target_manifest_sha256": sha256_file(target_path),
                "target_count": len(prepared_targets),
                "pair_status_sha256": sha256_file(pair_status_path),
            }
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
            failures.append(f"Embedding state evidence binding failed: {exc}")
    if args.require_embedding_evidence and not evidence_bound:
        failures.append(
            "Required embedding evidence must include coverage.json, contract.json, "
            "targets.tsv and pair_status.tsv bound to this benchmark and cache"
        )
    sequence_summary = summaries["sequence"]
    if sequence_summary["valid"] != len(targets):
        failures.append(
            "Sequence embeddings must be valid for every benchmark protein: "
            f"{sequence_summary['valid']}/{len(targets)}"
        )
    for modality in active:
        if modality != "sequence":
            summary = summaries[modality]
            if summary["valid"] == 0:
                failures.append(
                    f"{args.mode} mode has zero valid {modality} embeddings"
                )
            if summary["invalid_present"]:
                failures.append(
                    f"{modality} contains {summary['invalid_present']} malformed present arrays"
                )

    split_coverage: Dict[str, Any] = {}
    for (aspect, split), protein_ids in memberships.items():
        key = f"{aspect}_{split}"
        split_coverage[key] = {}
        protein_set = set(protein_ids)
        for modality in active:
            valid_rows = len(protein_set & valid_ids[modality])
            split_coverage[key][modality] = {
                "valid": valid_rows,
                "total": len(protein_ids),
                "coverage": valid_rows / len(protein_ids),
            }
    normalization_samples: Dict[str, Any] = {}
    if len(active) > 1:
        for aspect in aspects:
            sample_ids = memberships[(aspect, "train")][:1000]
            sample_set = set(sample_ids)
            normalization_samples[aspect] = {}
            for modality in active:
                valid_rows = len(sample_set & valid_ids[modality])
                normalization_samples[aspect][modality] = {
                    "valid": valid_rows,
                    "sampled_training_rows": len(sample_ids),
                }
                if modality != "sequence" and valid_rows == 0:
                    failures.append(
                        f"{aspect} has zero valid {modality} embeddings in PFP's first "
                        "1,000 training rows used for normalization"
                    )

    ia_inputs, ia_failures = validate_ia_inputs(
        args.data_dir,
        args.ia_file_dir.resolve() if args.ia_file_dir else None,
        aspects,
        bool(config.get("evaluation", {}).get("require_precomputed_ia", False)),
    )
    failures.extend(ia_failures)

    args.issues_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.issues_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("protein_id", "modality", "reason", "detail", "path"))
        writer.writerows(issues)

    report = {
        "schema_version": 1,
        "status": "failed" if failures else "passed",
        "mode": args.mode,
        "cache_root": str(cache_root),
        "target_count": len(targets),
        "modalities": summaries,
        "split_coverage": split_coverage,
        "normalization_samples": normalization_samples,
        "failures": failures,
        "embedding_evidence": evidence,
        "embedding_evidence_binding": binding,
        "provenance_evidence_supplied": bool(evidence),
        "provenance_evidence_bound": evidence_bound,
        "information_accretion": ia_inputs,
        "preparation_report": {
            **preparation_snapshot,
            "benchmark_fingerprint": preparation.get("benchmark_fingerprint"),
        },
        "policy": {
            "sequence": "100% valid sequence coverage is mandatory",
            "non_sequence_full_mode": (
                "Missing arrays are retained as PFP masks; every present array must be valid, "
                "and each enabled modality must occur in each selected aspect's PFP "
                "training normalization sample"
            ),
            "sequence_only": "Text, structure and PPI are intentionally disabled and not read",
        },
    }
    atomic_write_json(args.report, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    if failures:
        raise SystemExit("Embedding cache validation failed: " + "; ".join(failures))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
