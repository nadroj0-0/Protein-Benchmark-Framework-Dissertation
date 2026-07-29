#!/usr/bin/env python3
"""Reconstruct policy-bound temporal knowledge cohorts for a benchmark snapshot."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence


FRAMEWORK_ROOT = Path(__file__).resolve().parents[2]
BUILDER_SRC = FRAMEWORK_ROOT / "benchmark_builders" / "contemporary_cafa" / "src"
if str(BUILDER_SRC) not in sys.path:
    sys.path.insert(0, str(BUILDER_SRC))

from cafa_benchmark_builder.goa import load_normalized_annotation_map  # noqa: E402
from cafa_benchmark_builder.ontology import Ontology  # noqa: E402
from cafa_benchmark_builder.parsers import load_protein_catalog  # noqa: E402
from cafa_benchmark_builder.snapshot import (  # noqa: E402
    _build_identity_crosswalk,
    _drop_protein_binding_only,
)
from label_space_common import (  # noqa: E402
    ASPECT_TO_NAMESPACE,
    atomic_write_json,
    atomic_write_text,
    output_manifest,
    peak_rss_bytes,
    sha256_file,
)


ASPECTS = ("BPO", "CCO", "MFO")
PREFIX_TO_ASPECT = {"bp": "BPO", "cc": "CCO", "mf": "MFO"}
NAMESPACE_TO_ASPECT = {
    namespace: aspect for aspect, namespace in ASPECT_TO_NAMESPACE.items()
}
EXPOSURE_FIELDS = (
    "protein_id",
    "train_id_member",
    "valid_id_member",
    "train_sequence_member",
    "valid_sequence_member",
    "train_homology_cluster_member",
    "modality_availability",
    "feature_temporal_policy",
)


def _require_file(path: Path, label: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise ValueError(f"{label} is missing or empty: {resolved}")
    return resolved


def _validate_manifest(manifest: Mapping[str, Any]) -> None:
    expected = {
        "profile": "supervisor",
        "test_eligibility_policy": "global-no-knowledge",
        "target_universe_policy": "reconstructed-all-qualifying",
        "training_reviewed_only": False,
        "target_reviewed_only": False,
        "exclude_t1_backfill": False,
        "t1_endpoint_policy": "snapshot-membership",
        "require_t0_presence": True,
        "sequence_change_policy": "exclude",
        "protein_binding_policy": "drop-mf-protein-binding-only",
        "include_relationships": True,
    }
    mismatches = {
        key: {"expected": value, "observed": manifest.get(key)}
        for key, value in expected.items()
        if manifest.get(key) != value
    }
    if mismatches:
        raise ValueError(
            "Accepted benchmark manifest does not match the locked supervisor "
            f"cohort policy: {mismatches}"
        )
    target_taxa = tuple(str(value) for value in manifest.get("target_taxa", []))
    training_taxa = tuple(str(value) for value in manifest.get("training_taxa", []))
    if not target_taxa or target_taxa != training_taxa:
        raise ValueError("Supervisor benchmark must use the same non-empty taxa in both snapshots")
    evidence_codes = tuple(str(value) for value in manifest.get("evidence_codes", []))
    if not evidence_codes:
        raise ValueError("Accepted benchmark manifest has no evidence-code policy")


def _source_contract(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path),
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _split_annotations(
    annotations: Mapping[str, set[str]], ontology: Ontology
) -> dict[str, dict[str, frozenset[str]]]:
    output: dict[str, dict[str, frozenset[str]]] = {}
    for protein_id, terms in annotations.items():
        by_aspect = {aspect: set() for aspect in ASPECTS}
        for term in terms:
            namespace = ontology.get_namespace(term)
            aspect = NAMESPACE_TO_ASPECT.get(namespace)
            if aspect is None:
                raise ValueError(f"GO term has an unsupported namespace: {term}/{namespace}")
            by_aspect[aspect].add(term)
        output[protein_id] = {
            aspect: frozenset(by_aspect[aspect]) for aspect in ASPECTS
        }
    return output


def _propagate_annotations(
    direct: Mapping[str, Mapping[str, frozenset[str]]], ontology: Ontology
) -> dict[str, dict[str, frozenset[str]]]:
    ancestor_cache: dict[str, frozenset[str]] = {}
    output: dict[str, dict[str, frozenset[str]]] = {}
    for protein_id, values in direct.items():
        propagated: dict[str, frozenset[str]] = {}
        for aspect in ASPECTS:
            terms: set[str] = set()
            for term in values.get(aspect, frozenset()):
                if term not in ancestor_cache:
                    ancestor_cache[term] = frozenset(ontology.get_ancestors(term))
                terms.update(ancestor_cache[term])
            propagated[aspect] = frozenset(terms)
        output[protein_id] = propagated
    return output


def _write_scope(path: Path, protein_ids: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id"])
        writer.writerows((protein_id,) for protein_id in protein_ids)
    os.replace(temporary, path)


def _write_annotations(
    path: Path,
    protein_ids: Sequence[str],
    annotations: Mapping[str, Mapping[str, frozenset[str]]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "aspect", "go_term"])
        for protein_id in protein_ids:
            for aspect in ASPECTS:
                for term in sorted(annotations.get(protein_id, {}).get(aspect, ())):
                    writer.writerow([protein_id, aspect, term])
    os.replace(temporary, path)


def _read_benchmark_split(
    path: Path,
) -> tuple[dict[str, str], dict[str, frozenset[str]], tuple[str, ...]]:
    proteins: dict[str, str] = {}
    positives: dict[str, frozenset[str]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Benchmark CSV is empty: {path}") from exc
        if len(header) < 3 or header[:2] != ["proteins", "sequences"]:
            raise ValueError(f"Unexpected benchmark CSV schema: {path}")
        terms = tuple(header[2:])
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Malformed benchmark row at {path}:{line_number}")
            protein_id, sequence = row[:2]
            if not protein_id or protein_id in proteins:
                raise ValueError(f"Empty or duplicate protein at {path}:{line_number}")
            proteins[protein_id] = sequence
            positives[protein_id] = frozenset(
                term for term, value in zip(terms, row[2:]) if value == "1"
            )
            invalid = sorted(set(row[2:]) - {"0", "1"})
            if invalid:
                raise ValueError(f"Non-binary label at {path}:{line_number}: {invalid}")
    return proteins, positives, terms


def _load_exposure_sets(
    benchmark_dir: Path,
) -> tuple[set[str], set[str], set[str], set[str], dict[str, Any]]:
    train_ids: set[str] = set()
    valid_ids: set[str] = set()
    train_sequences: set[str] = set()
    valid_sequences: set[str] = set()
    files: dict[str, dict[str, Any]] = {}
    for prefix in PREFIX_TO_ASPECT:
        for split, ids, sequences in (
            ("training", train_ids, train_sequences),
            ("validation", valid_ids, valid_sequences),
        ):
            path = _require_file(benchmark_dir / f"{prefix}-{split}.csv", f"{prefix} {split} CSV")
            proteins = _read_benchmark_membership(path)
            ids.update(proteins)
            sequences.update(proteins.values())
            files[path.name] = {
                **_source_contract(path),
                "proteins": len(proteins),
            }
    overlap = sorted(train_ids & valid_ids)
    if overlap:
        raise ValueError(f"Accepted training and validation IDs overlap: {overlap[:5]}")
    return train_ids, valid_ids, train_sequences, valid_sequences, files


def _read_benchmark_membership(path: Path) -> dict[str, str]:
    proteins: dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Benchmark CSV is empty: {path}") from exc
        if len(header) < 3 or header[:2] != ["proteins", "sequences"]:
            raise ValueError(f"Unexpected benchmark CSV schema: {path}")
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Malformed benchmark row at {path}:{line_number}")
            protein_id, sequence = row[:2]
            if not protein_id or protein_id in proteins:
                raise ValueError(f"Empty or duplicate protein at {path}:{line_number}")
            proteins[protein_id] = sequence
    return proteins


def _write_exposure(
    path: Path,
    protein_ids: Sequence[str],
    sequences: Mapping[str, str],
    train_ids: set[str],
    valid_ids: set[str],
    train_sequences: set[str],
    valid_sequences: set[str],
    feature_policy: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=EXPOSURE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for protein_id in protein_ids:
            sequence = sequences[protein_id]
            writer.writerow(
                {
                    "protein_id": protein_id,
                    "train_id_member": int(protein_id in train_ids),
                    "valid_id_member": int(protein_id in valid_ids),
                    "train_sequence_member": int(sequence in train_sequences),
                    "valid_sequence_member": int(sequence in valid_sequences),
                    "train_homology_cluster_member": "unknown",
                    "modality_availability": "not_assessed_census_only",
                    "feature_temporal_policy": feature_policy,
                }
            )
    os.replace(temporary, path)


def _verify_current_test_truth(
    benchmark_dir: Path,
    direct_t0: Mapping[str, Mapping[str, frozenset[str]]],
    closure_t1: Mapping[str, Mapping[str, frozenset[str]]],
) -> dict[str, Any]:
    report: dict[str, Any] = {}
    for prefix, aspect in PREFIX_TO_ASPECT.items():
        path = _require_file(benchmark_dir / f"{prefix}-test.csv", f"{prefix} test CSV")
        proteins, positives, terms = _read_benchmark_split(path)
        term_set = set(terms)
        mismatches: list[str] = []
        t0_known: list[str] = []
        for protein_id in proteins:
            if any(direct_t0.get(protein_id, {}).get(value, ()) for value in ASPECTS):
                t0_known.append(protein_id)
            observed = set(closure_t1.get(protein_id, {}).get(aspect, ())) & term_set
            if observed != set(positives[protein_id]):
                mismatches.append(protein_id)
        if t0_known:
            raise ValueError(
                f"Accepted {aspect} test contains qualifying t0 knowledge: {t0_known[:5]}"
            )
        if mismatches:
            raise ValueError(
                f"Reconstructed {aspect} t1 truth differs from accepted CSV for "
                f"{len(mismatches)} proteins: {mismatches[:5]}"
            )
        report[aspect] = {
            "test_proteins": len(proteins),
            "term_columns": len(terms),
            "t0_global_no_knowledge_verified": True,
            "t1_truth_exact_match_verified": True,
            "source": _source_contract(path),
        }
    return report


def _relationship_policy(path: Path) -> tuple[str, ...]:
    relationships = {"is_a"}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("relationship: "):
                fields = raw.split()
                if len(fields) >= 2:
                    relationships.add(fields[1])
    return tuple(sorted(relationships))


def _census_tsv(summary: Mapping[str, Any]) -> str:
    stream = io.StringIO(newline="")
    fields = ("aspect", "knowledge_state", "all_t1_annotated", "with_new_terms")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for aspect in ASPECTS:
        for state, count in summary["cohort_counts"][aspect].items():
            writer.writerow(
                {
                    "aspect": aspect,
                    "knowledge_state": state,
                    "all_t1_annotated": count,
                    "with_new_terms": summary["gainer_cohort_counts"][aspect][state],
                }
            )
    return stream.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted-benchmark-dir", type=Path, required=True)
    parser.add_argument("--t0-sprot", type=Path, required=True)
    parser.add_argument("--t0-trembl", type=Path, required=True)
    parser.add_argument("--t1-sprot", type=Path, required=True)
    parser.add_argument("--t1-trembl", type=Path, required=True)
    parser.add_argument("--goa-t0", type=Path, required=True)
    parser.add_argument("--goa-t1", type=Path, required=True)
    parser.add_argument("--benchmark-obo", type=Path, required=True)
    parser.add_argument("--t0-source-obo", type=Path, required=True)
    parser.add_argument("--t1-source-obo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--benchmark-id",
        default="contemporary/2025_01_to_2026_02_supervisor",
    )
    parser.add_argument(
        "--feature-temporal-policy",
        default="text-cutoff-2025-03-08__ppi-paper-faithful",
    )
    args = parser.parse_args()

    started = time.perf_counter()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")
    benchmark_dir = args.accepted_benchmark_dir.resolve()
    manifest_path = _require_file(benchmark_dir / "build_manifest.json", "accepted build manifest")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    _validate_manifest(manifest)

    sources = {
        "t0_sprot": _require_file(args.t0_sprot, "t0 Swiss-Prot"),
        "t0_trembl": _require_file(args.t0_trembl, "t0 TrEMBL subset"),
        "t1_sprot": _require_file(args.t1_sprot, "t1 Swiss-Prot"),
        "t1_trembl": _require_file(args.t1_trembl, "t1 TrEMBL subset"),
        "goa_t0": _require_file(args.goa_t0, "t0 GOA"),
        "goa_t1": _require_file(args.goa_t1, "t1 GOA"),
        "benchmark_obo": _require_file(args.benchmark_obo, "frozen benchmark ontology"),
        "t0_source_obo": _require_file(args.t0_source_obo, "t0 source ontology"),
        "t1_source_obo": _require_file(args.t1_source_obo, "t1 source ontology"),
    }
    relationships = _relationship_policy(sources["benchmark_obo"])
    if relationships == ("is_a",):
        raise ValueError("Accepted all-relationships policy found no OBO relationship edges")

    benchmark_go = Ontology(sources["benchmark_obo"], with_rels=True)
    t0_go = Ontology(sources["t0_source_obo"], with_rels=True)
    t1_go = Ontology(sources["t1_source_obo"], with_rels=True)
    taxa = frozenset(str(value) for value in manifest["target_taxa"])
    evidence_codes = frozenset(str(value) for value in manifest["evidence_codes"])

    t0_catalog = load_protein_catalog(
        (sources["t0_sprot"], sources["t0_trembl"]), taxa, False
    )
    t1_catalog = load_protein_catalog(
        (sources["t1_sprot"], sources["t1_trembl"]), taxa, False
    )
    matches, t1_to_t0 = _build_identity_crosswalk(t0_catalog, t1_catalog, "exclude")
    match_counts = Counter((match.status, match.reason) for match in matches)

    t0_result = load_normalized_annotation_map(
        sources["goa_t0"],
        alias_to_primary=t0_catalog.alias_to_primary,
        source_ontology=t0_go,
        benchmark_ontology=benchmark_go,
        other_ontology=t1_go,
        snapshot="t0",
        allow_frozen_source_fallback=bool(manifest["allow_frozen_source_fallback"]),
        evidence_codes=evidence_codes,
    )
    t1_result = load_normalized_annotation_map(
        sources["goa_t1"],
        alias_to_primary=t1_catalog.alias_to_primary,
        source_ontology=t1_go,
        benchmark_ontology=benchmark_go,
        other_ontology=t0_go,
        snapshot="t1",
        allow_frozen_source_fallback=bool(manifest["allow_frozen_source_fallback"]),
        evidence_codes=evidence_codes,
        target_taxa=taxa,
    )
    if t0_result.unmapped_terms or t1_result.unmapped_terms:
        raise ValueError("Strict cohort reconstruction found unresolved source GO IDs")

    translated_t1 = {
        t1_to_t0[t1_id]: _drop_protein_binding_only(
            set(terms), benchmark_go, str(manifest["protein_binding_policy"])
        )
        for t1_id, terms in t1_result.annotations.items()
        if t1_id in t1_to_t0
    }
    scope = tuple(sorted(protein_id for protein_id, terms in translated_t1.items() if terms))
    if not scope:
        raise ValueError("No matched t1-annotated proteins remain in cohort scope")
    scope_set = set(scope)
    t0_terms = {
        protein_id: set(t0_result.annotations.get(protein_id, set()))
        for protein_id in scope
        if t0_result.annotations.get(protein_id)
    }
    t1_terms = {protein_id: set(translated_t1[protein_id]) for protein_id in scope}
    t0_direct = _split_annotations(t0_terms, benchmark_go)
    t1_direct = _split_annotations(t1_terms, benchmark_go)
    t0_closure = _propagate_annotations(t0_direct, benchmark_go)
    t1_closure = _propagate_annotations(t1_direct, benchmark_go)

    train_ids, valid_ids, train_sequences, valid_sequences, exposure_sources = (
        _load_exposure_sets(benchmark_dir)
    )
    current_test_alignment = _verify_current_test_truth(
        benchmark_dir, t0_direct, t1_closure
    )
    accepted_test_ids: set[str] = set()
    for prefix in PREFIX_TO_ASPECT:
        proteins, _, _ = _read_benchmark_split(benchmark_dir / f"{prefix}-test.csv")
        accepted_test_ids.update(proteins)
    if not accepted_test_ids <= scope_set:
        missing = sorted(accepted_test_ids - scope_set)
        raise ValueError(f"Accepted test proteins are absent from cohort scope: {missing[:5]}")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.stage-", dir=output_dir.parent))
    input_dir = Path(tempfile.mkdtemp(prefix="cohort-inputs-", dir=stage.parent))
    ledger_dir = stage / "ledger"
    try:
        _write_scope(input_dir / "protein_scope.tsv", scope)
        _write_scope(input_dir / "t0_presence.tsv", scope)
        _write_scope(input_dir / "t1_presence.tsv", scope)
        _write_annotations(input_dir / "t0_direct.tsv", scope, t0_direct)
        _write_annotations(input_dir / "t1_direct.tsv", scope, t1_direct)
        _write_annotations(input_dir / "t0_closure.tsv", scope, t0_closure)
        _write_annotations(input_dir / "t1_closure.tsv", scope, t1_closure)
        _write_exposure(
            input_dir / "exposure.tsv",
            scope,
            t0_catalog.sequences,
            train_ids,
            valid_ids,
            train_sequences,
            valid_sequences,
            args.feature_temporal_policy,
        )
        command = [
            sys.executable,
            str(FRAMEWORK_ROOT / "scripts" / "diagnostics" / "build_temporal_annotation_ledger.py"),
            "--t0-direct-annotations", str(input_dir / "t0_direct.tsv"),
            "--t1-direct-annotations", str(input_dir / "t1_direct.tsv"),
            "--t0-closure-annotations", str(input_dir / "t0_closure.tsv"),
            "--t1-closure-annotations", str(input_dir / "t1_closure.tsv"),
            "--t0-protein-presence", str(input_dir / "t0_presence.tsv"),
            "--t1-protein-presence", str(input_dir / "t1_presence.tsv"),
            "--exposure-table", str(input_dir / "exposure.tsv"),
            "--protein-scope", str(input_dir / "protein_scope.tsv"),
            "--output-dir", str(ledger_dir),
            "--t0-snapshot", "UniProt-2025_01/GOA-225",
            "--t1-snapshot", "UniProt-2026_02/GOA-234",
            "--evidence-policy-id", "supervisor-explicit-17-code-policy",
            "--graph-policy-id", "accepted-frozen-t0-all-relationships",
            "--benchmark-id", args.benchmark_id,
        ]
        for relationship in relationships:
            command.extend(("--relationship", relationship))
        subprocess.run(command, check=True)
        for child in ledger_dir.iterdir():
            shutil.move(str(child), stage / child.name)
        ledger_dir.rmdir()

        summary_path = stage / "temporal_annotation_ledger.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        preparation = {
            "schema_version": 1,
            "status": "complete",
            "analysis_kind": "contemporary_knowledge_cohort_census_preparation",
            "scope_policy": "matched unchanged t0/t1 proteins with at least one policy-eligible direct t1 term",
            "accepted_benchmark_manifest": {
                "path": str(manifest_path),
                "sha256": sha256_file(manifest_path),
            },
            "source_files": {key: _source_contract(value) for key, value in sources.items()},
            "source_policy": {
                "profile": manifest["profile"],
                "evidence_codes": sorted(evidence_codes),
                "test_eligibility_policy": manifest["test_eligibility_policy"],
                "protein_binding_policy": manifest["protein_binding_policy"],
                "graph_relationships": list(relationships),
                "graph_semantics": "exactly mirrors the accepted benchmark builder, including every relationship line",
            },
            "catalog_counts": {
                "t0_proteins": len(t0_catalog.records),
                "t1_proteins": len(t1_catalog.records),
                "matched_unchanged": len(t1_to_t0),
                "cohort_scope": len(scope),
            },
            "identity_match_counts": {
                f"{status}:{reason}": count
                for (status, reason), count in sorted(match_counts.items())
            },
            "annotation_load_counts": {
                "t0": dict(sorted(t0_result.counters.items())),
                "t1": dict(sorted(t1_result.counters.items())),
            },
            "exposure_sources": exposure_sources,
            "accepted_test_alignment": current_test_alignment,
            "resource_usage": {
                "wall_seconds": time.perf_counter() - started,
                "peak_rss_bytes": peak_rss_bytes(),
            },
        }
        summary["preparation"] = preparation
        atomic_write_json(summary_path, summary)
        atomic_write_json(stage / "cohort_census.json", preparation)
        atomic_write_text(stage / "cohort_census.tsv", _census_tsv(summary))
        (stage / "output_manifest.json").unlink(missing_ok=True)
        (stage / "RUN_COMPLETE.json").unlink(missing_ok=True)
        artifacts = output_manifest(stage, exclude={"output_manifest.json", "RUN_COMPLETE.json"})
        atomic_write_json(stage / "output_manifest.json", artifacts)
        atomic_write_json(
            stage / "RUN_COMPLETE.json",
            {
                "schema_version": 1,
                "complete": True,
                "benchmark_id": args.benchmark_id,
                "analysis_kind": "contemporary_knowledge_cohort_census",
                "output_manifest_sha256": sha256_file(stage / "output_manifest.json"),
            },
        )
        os.replace(stage, output_dir)
    finally:
        if input_dir.exists():
            shutil.rmtree(input_dir)
        if stage.exists():
            shutil.rmtree(stage)

    print(f"Published contemporary cohort census: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
