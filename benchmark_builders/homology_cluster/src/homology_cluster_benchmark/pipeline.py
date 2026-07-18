from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
import csv
from dataclasses import replace
import gzip
import io
import json
import logging
import math
import hashlib
import os
from pathlib import Path
import re
import shutil
import statistics
import time
from typing import Iterable
import uuid

from .attrition import (
    METRIC_DEFINITIONS,
    evaluate_attrition,
    load_attrition_policy,
    observation,
)
from .clustering import (
    connect_proteins_to_clusters,
    mapping_counters,
    mapping_counters_by_source,
    retained_cluster_info,
)
from .common_cache import (
    CACHE_MARKER,
    common_cache_root,
    inspect_common_preprocessing_cache,
    load_common_preprocessing_cache,
)
from .config import PREFIX_TO_NAMESPACE, SPLITS, SUPPORTED_IDENTITIES, BuildConfig
from .export import export_all
from .frozen_inputs import (
    FrozenInputManifest,
    bind_frozen_inputs,
    load_frozen_input_manifest,
    verify_frozen_manifest_unchanged,
    write_synthetic_fixture_manifest,
)
from .goa import iter_annotation_records, iter_excluded_annotations, load_goa
from .idmapping import load_uniref90_mappings
from .inputs import resolve_input, sha256_file
from .labels import build_labels
from .mapping import canonicalize_goa_accessions, load_requested_proteins_from_sources
from .models import BuildResult, MappingDecision, ResolvedInput
from .mmseqs import (
    ClusterIndex,
    build_mmseqs_commands,
    execute_commands,
    resolve_mmseqs_runtime,
    validate_exact_mmseqs_version,
    validate_recorded_exact_mmseqs_version,
    validate_mmseqs_version,
    verify_mmseqs_executable_unchanged,
    write_command_manifest,
)
from .ontology import Ontology
from .provenance import (
    PUBLICATION_MARKER_KEYS,
    git_state,
    publish,
    runtime_provenance,
    staging_output,
    verify_output_manifest,
    write_completion_marker,
    write_output_manifest,
    write_publication_metadata,
    utc_now,
)
from .splitting import assign_development_test, assign_training_validation
from .uniref import UniRefIndex
from .validation import (
    split_balance_metrics,
    validate_outputs,
    validate_term_universe_artifacts,
    write_validation_report,
)


LOGGER = logging.getLogger(__name__)

ROOT_POLICY = "retain GO:0008150, GO:0005575, GO:0003674 for PFP/TEMPROT compatibility"
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def _json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


@contextmanager
def _output_text(path: Path):
    if path.suffix != ".gz":
        with path.open("w", encoding="utf-8", newline="") as handle:
            yield handle
        return
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                yield handle


def _tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with _output_text(path) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _parameters(config: BuildConfig) -> dict[str, object]:
    return {
        "identity_fraction": config.identity,
        "identity_percent": int(config.identity * 100),
        "coverage": config.coverage,
        "coverage_interpretation": "alignment residues / max(query length, target length) >= 0.8",
        "cov_mode": config.cov_mode,
        "cluster_mode": config.cluster_mode,
        "alignment_mode": config.alignment_mode,
        "seq_id_mode": 0,
        "createdb_shuffle": 0,
        "cluster_reassign": config.cluster_reassign,
        "sensitivity": config.sensitivity,
        "evalue": config.evalue,
        "threads": config.threads,
        "requested_slots": config.requested_slots,
        "allocated_slots": config.allocated_slots,
        "run_id": config.run_id,
        "uniprot_source_scope": config.uniprot_source_scope,
        "framework_revision": config.framework_revision,
        "benchmark_scope": config.benchmark_scope,
        "split_policy": config.split_policy,
        "development_fraction": config.development_fraction,
        "training_fraction_within_development": config.training_fraction_within_development,
        "seed": config.seed,
        "training_population": config.training_population,
        "min_count": config.min_count,
        "evidence_codes": sorted(config.evidence_codes),
        "include_relationships": config.include_relationships,
        "root_policy": ROOT_POLICY,
        "uniprot_release": config.release_uniprot,
        "goa_release": config.release_goa,
        "ontology_release": config.release_ontology,
        "fixture_mode": config.fixture_mode,
        "precomputed_cluster_assignments": str(config.cluster_assignments.resolve()) if config.cluster_assignments else None,
        "scratch_safety_multiplier": config.scratch_safety_multiplier,
        "minimum_free_disk_bytes": config.minimum_free_disk_bytes,
        "mmseqs_work_multiplier": config.mmseqs_work_multiplier,
        "publication_safety_multiplier": config.publication_safety_multiplier,
        "excluded_sample_per_reason": config.excluded_sample_per_reason,
        "expected_mmseqs_version": config.expected_mmseqs_version,
    }


def _scientific_fingerprint_payload(
    config: BuildConfig,
    frozen_manifest: FrozenInputManifest,
    mmseqs_runtime,
    repository_commit: str | None,
    attrition_policy_sha256: str,
) -> dict[str, object]:
    """Return the identity-independent contract shared by all six threshold jobs."""
    return {
        "frozen_input_source_fingerprint": frozen_manifest.source_fingerprint,
        "uniprot_source_scope": config.uniprot_source_scope,
        "attrition_policy_sha256": attrition_policy_sha256,
        "uniprot_release": config.release_uniprot,
        "goa_release": config.release_goa,
        "ontology_release": config.release_ontology,
        "split_policy": config.split_policy,
        "development_fraction": config.development_fraction,
        "training_fraction_within_development": (
            config.training_fraction_within_development
        ),
        "training_population": config.training_population,
        "seed": config.seed,
        "min_count": config.min_count,
        "evidence_codes": sorted(config.evidence_codes),
        "include_relationships": config.include_relationships,
        "root_policy": ROOT_POLICY,
        "coverage": config.coverage,
        "cov_mode": config.cov_mode,
        "cluster_mode": config.cluster_mode,
        "alignment_mode": config.alignment_mode,
        "seq_id_mode": 0,
        "createdb_shuffle": 0,
        "cluster_reassign": config.cluster_reassign,
        "sensitivity": config.sensitivity,
        "evalue": config.evalue,
        "expected_mmseqs_version": config.expected_mmseqs_version,
        "observed_mmseqs_version": mmseqs_runtime.version_token,
        "mmseqs_executable_sha256": mmseqs_runtime.executable_sha256,
        "repository_commit": repository_commit,
        "framework_revision": config.framework_revision or repository_commit,
        "requested_slots": config.requested_slots,
        "allocated_slots": config.allocated_slots,
        "mmseqs_threads": config.threads,
    }


def _scientific_fingerprint(payload: dict[str, object]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _validate_frozen_metadata(config: BuildConfig, ontology: Ontology, headers: dict[str, str]) -> None:
    if config.release_ontology and ontology.data_version != config.release_ontology:
        raise ValueError(
            "Frozen GO ontology version mismatch: "
            f"expected {config.release_ontology!r}, observed {ontology.data_version!r}"
        )
    if config.release_goa == "234":
        date_generated = headers.get("date_generated", "")
        go_version = headers.get("go_version", "")
        if not date_generated.startswith("2026-06-17"):
            raise ValueError(
                "GOA 234 convention requires !date-generated: 2026-06-17; "
                f"observed {date_generated!r}"
            )
        if "2026-06-15" not in go_version:
            raise ValueError(
                "GOA 234 convention requires a 2026-06-15 GO version; "
                f"observed {go_version!r}"
            )


def _resolved_inputs(
    config: BuildConfig,
    work: Path,
    frozen_manifest: FrozenInputManifest | None,
) -> dict[str, ResolvedInput]:
    specs = _input_specs(config)
    if config.common_preprocessing_cache is not None:
        if frozen_manifest is None:
            raise ValueError(
                "A common preprocessing cache requires a frozen-input manifest"
            )
        expected_hashes = {
            name: str(entry["sha256"])
            for name, entry in frozen_manifest.entries.items()
        }
        cache_root = common_cache_root(config.common_preprocessing_cache)
        inspect_common_preprocessing_cache(
            cache_root,
            expected_source_scope=config.uniprot_source_scope,
            expected_input_sha256=expected_hashes,
            verify_file_hashes=False,
        )
        resolved = {
            name: resolve_input(specs[name], work / "downloads", config.allow_downloads)
            for name in ("uniref90_fasta", "go_obo")
        }
        marker = cache_root / CACHE_MARKER
        for name, spec in specs.items():
            if name in resolved:
                continue
            entry = frozen_manifest.entries[name]
            resolved[name] = ResolvedInput(
                name=name,
                resolved_path=marker,
                source_url=str(entry["url"]),
                release=str(entry["release"]),
                size_bytes=int(entry["size_bytes"]),
                sha256=str(entry["sha256"]),
                expected_sha256=spec.expected_sha256,
                acquisition="authenticated-common-preprocessing-cache",
                source_population=str(entry["source_population"]),
            )
        return resolved
    return {
        name: resolve_input(spec, work / "downloads", config.allow_downloads)
        for name, spec in specs.items()
    }


def _input_specs(config: BuildConfig):
    specs = {
        "uniref90_fasta": config.uniref90_fasta,
        "idmapping": config.idmapping,
        "goa": config.goa,
        "go_obo": config.go_obo,
    }
    for name in config.selected_uniprot_input_names:
        spec = getattr(config, name)
        assert spec is not None
        specs[name] = spec
    return specs


def _validate_manifest_embedded_metadata(
    manifest: FrozenInputManifest, ontology: Ontology, goa_headers: dict[str, str]
) -> None:
    declared_goa = manifest.entries["goa"]["embedded_metadata"]
    for key, expected in declared_goa.items():
        observed = goa_headers.get(key)
        if expected and str(expected) not in str(observed):
            raise ValueError(
                f"Frozen-input manifest embedded GOA metadata mismatch for {key}: "
                f"expected {expected!r}, observed {observed!r}"
            )
    declared_ontology = manifest.entries["go_obo"]["embedded_metadata"]
    expected_data_version = declared_ontology.get("data_version")
    if expected_data_version and ontology.data_version != expected_data_version:
        raise ValueError(
            "Frozen-input manifest ontology data_version mismatch: "
            f"expected {expected_data_version!r}, observed {ontology.data_version!r}"
        )


def _disk_preflight(
    config: BuildConfig, work: Path, inputs: dict[str, ResolvedInput]
) -> dict[str, object]:
    scratch_usage = shutil.disk_usage(work)
    publication_root = config.output_dir.expanduser().resolve()
    publication_root.mkdir(parents=True, exist_ok=True)
    publication_usage = shutil.disk_usage(publication_root)
    logical_input_bytes = sum(item.size_bytes for item in inputs.values())
    uniref_bytes = inputs["uniref90_fasta"].size_bytes
    common_cache_bytes = 0
    if config.common_preprocessing_cache is not None:
        cache_root = common_cache_root(config.common_preprocessing_cache)
        common_cache_bytes = sum(
            path.stat().st_size for path in cache_root.rglob("*") if path.is_file()
        )
        staged_input_bytes = (
            uniref_bytes + inputs["go_obo"].size_bytes + common_cache_bytes
        )
        goa_bytes = 0
    else:
        staged_input_bytes = logical_input_bytes
        goa_bytes = inputs["goa"].size_bytes
    mmseqs_work_estimate = int(uniref_bytes * config.mmseqs_work_multiplier)
    parser_index_estimate = (
        0 if config.common_preprocessing_cache is not None else int(
            (goa_bytes + inputs["idmapping"].size_bytes)
            * config.scratch_safety_multiplier
        )
    )
    estimated_scratch = (
        staged_input_bytes
        + mmseqs_work_estimate
        + parser_index_estimate
        + config.minimum_free_disk_bytes
    )
    if config.common_preprocessing_cache is not None:
        # The cache is an input, not a publication. Cluster/member reports dominate the
        # threshold output, so use the UniRef scaffold as the conservative output basis.
        publication_basis = uniref_bytes
    else:
        publication_basis = (
            uniref_bytes
            + goa_bytes
            + sum(inputs[name].size_bytes for name in config.selected_uniprot_input_names)
        )
    estimated_publication = int(
        publication_basis * config.publication_safety_multiplier
        + config.minimum_free_disk_bytes
    )
    same_filesystem = os.stat(work).st_dev == os.stat(publication_root).st_dev
    required_scratch_filesystem = (
        estimated_scratch + estimated_publication if same_filesystem else estimated_scratch
    )
    estimate_enforced = not config.diagnostic_pilot
    estimate_exceeds_available_space = scratch_usage.free < required_scratch_filesystem
    if estimate_enforced and estimate_exceeds_available_space:
        raise OSError(
            "Scratch preflight failed: "
            f"free={scratch_usage.free} bytes, estimate={required_scratch_filesystem} bytes. "
            "Increase scratch, reduce only the explicit safety multiplier with measured evidence, "
            "or set a site-specific minimum."
        )
    if estimate_enforced and not same_filesystem and publication_usage.free < estimated_publication:
        raise OSError(
            "Publication-filesystem preflight failed: "
            f"free={publication_usage.free} bytes, estimate={estimated_publication} bytes."
        )
    persistent_path = (
        config.persistent_results_root.expanduser().resolve()
        if config.persistent_results_root else publication_root
    )
    persistent_path.mkdir(parents=True, exist_ok=True)
    persistent_usage = shutil.disk_usage(persistent_path)
    if estimate_enforced and persistent_usage.free < estimated_publication:
        raise OSError(
            "Persistent-results preflight failed before benchmark execution: "
            f"free={persistent_usage.free} bytes, estimate={estimated_publication} bytes."
        )
    return {
        "scratch_free_bytes_at_preflight": scratch_usage.free,
        "publication_free_bytes_at_preflight": publication_usage.free,
        "same_filesystem": same_filesystem,
        "logical_input_bytes": logical_input_bytes,
        "staged_input_bytes": staged_input_bytes,
        "common_preprocessing_cache_bytes": common_cache_bytes,
        "common_preprocessing_cache_used": config.common_preprocessing_cache is not None,
        "uniref90_bytes": uniref_bytes,
        "goa_bytes": goa_bytes,
        "persistent_results_root": str(persistent_path),
        "persistent_free_bytes_at_preflight": persistent_usage.free,
        "scratch_safety_multiplier": config.scratch_safety_multiplier,
        "mmseqs_work_multiplier": config.mmseqs_work_multiplier,
        "publication_safety_multiplier": config.publication_safety_multiplier,
        "minimum_free_disk_bytes": config.minimum_free_disk_bytes,
        "estimated_mmseqs_work_bytes": mmseqs_work_estimate,
        "estimated_parser_index_bytes": parser_index_estimate,
        "estimated_scratch_bytes": estimated_scratch,
        "estimated_publication_bytes": estimated_publication,
        "required_on_scratch_filesystem_bytes": required_scratch_filesystem,
        "estimate_is_exact": False,
        "estimate_enforced": estimate_enforced,
        "estimate_exceeds_available_space": estimate_exceeds_available_space,
        "measurement_note": (
            "Diagnostic pilots record but do not enforce speculative multiplier estimates; "
            "use the runtime disk-usage report to calibrate production requests."
            if config.diagnostic_pilot else
            "Production runs enforce the configured conservative estimate."
        ),
    }


def _verify_inputs_unchanged(inputs: dict[str, ResolvedInput]) -> None:
    for name, item in sorted(inputs.items()):
        if item.acquisition == "authenticated-common-preprocessing-cache":
            continue
        path = item.resolved_path
        if not path.is_file():
            raise ValueError(f"Input changed during the run: {name} disappeared from {path}")
        if path.stat().st_size != item.size_bytes or sha256_file(path) != item.sha256:
            raise ValueError(f"Input changed during the run after initial hashing: {name} ({path})")


def _preserve_failure_logs(config: BuildConfig, work: Path, error: BaseException) -> Path | None:
    log_source = work / "logs"
    if not log_source.exists():
        return None
    destination = (
        config.output_dir.expanduser().resolve()
        / "_failed_runs"
        / (
            f"{config.identity_directory}-seed-{config.seed}-min-count-{config.min_count}-"
            f"{uuid.uuid4().hex}"
        )
    )
    try:
        destination.mkdir(parents=True, exist_ok=False)
        shutil.copytree(log_source, destination / "logs")
        _json(destination / "FAILURE.json", {
            "failed_at": utc_now(),
            "identity_percent": int(config.identity * 100),
            "fixture_mode": config.fixture_mode,
            "error_type": type(error).__name__,
            "error": str(error),
            "completion_marker_written": False,
        })
    except OSError:
        return None
    return destination


def _decision_index(decisions: list[MappingDecision]) -> tuple[dict[str, MappingDecision], dict[str, MappingDecision]]:
    by_raw = {decision.raw_accession: decision for decision in decisions}
    by_protein: dict[str, MappingDecision] = {}
    for decision in decisions:
        existing = by_protein.get(decision.protein_id)
        if existing is None or (not existing.mmseqs_cluster_id and decision.mmseqs_cluster_id):
            by_protein[decision.protein_id] = decision
    return by_raw, by_protein


def _write_mapping_manifests(stage: Path, decisions: list[MappingDecision]) -> None:
    fields = [
        "raw_uniprot_accession", "uniprot_accession", "accession_action", "uniref90_id",
        "mapping_status", "mapping_detail", "exists_in_uniref90_fasta",
        "canonical_sequence_available", "accession_lifecycle_status", "source_population",
    ]
    _tsv(stage / "uniprot_to_uniref90.tsv", fields, (
        {
            "raw_uniprot_accession": item.raw_accession,
            "uniprot_accession": item.protein_id,
            "accession_action": item.accession_action,
            "uniref90_id": item.uniref90_id,
            "mapping_status": item.status,
            "mapping_detail": item.detail,
            "exists_in_uniref90_fasta": (
                "" if item.exists_in_fasta is None else int(item.exists_in_fasta)
            ),
            "canonical_sequence_available": int(item.canonical_sequence_available),
            "accession_lifecycle_status": item.accession_lifecycle_status,
            "source_population": item.source_population,
        }
        for item in sorted(decisions, key=lambda row: (row.protein_id, row.raw_accession))
    ))
    assignment_fields = fields + ["mmseqs_cluster_id", "split"]
    _tsv(stage / "protein_cluster_assignments.tsv", assignment_fields, (
        {
            "raw_uniprot_accession": item.raw_accession,
            "uniprot_accession": item.protein_id,
            "accession_action": item.accession_action,
            "uniref90_id": item.uniref90_id,
            "mapping_status": item.status,
            "mapping_detail": item.detail,
            "exists_in_uniref90_fasta": (
                "" if item.exists_in_fasta is None else int(item.exists_in_fasta)
            ),
            "canonical_sequence_available": int(item.canonical_sequence_available),
            "accession_lifecycle_status": item.accession_lifecycle_status,
            "source_population": item.source_population,
            "mmseqs_cluster_id": item.mmseqs_cluster_id,
            "split": item.split,
        }
        for item in sorted(decisions, key=lambda row: (row.protein_id, row.raw_accession))
    ))


def _write_cluster_manifests(
    stage: Path,
    cluster_index: ClusterIndex,
    uniref: UniRefIndex,
    retained: dict,
    assignments: dict,
    decisions: list[MappingDecision],
    giant_threshold: int,
) -> tuple[int, int]:
    retained_ids = set(retained)
    annotated_uniref = {
        item.uniref90_id for item in decisions
        if item.mmseqs_cluster_id in retained_ids and item.uniref90_id
    }
    _tsv(stage / "mmseqs_cluster_membership.tsv.gz", ["mmseqs_cluster_id", "uniref90_id"], (
        {"mmseqs_cluster_id": cluster, "uniref90_id": member}
        for cluster, member in cluster_index.iter_assignments()
    ))
    _tsv(stage / "cluster_split_assignments.tsv", [
        "mmseqs_cluster_id", "split", "uniref90_member_count",
        "qualifying_uniprot_count", "assignment_stage",
    ], (
        {
            "mmseqs_cluster_id": item.cluster_id,
            "split": item.split,
            "uniref90_member_count": item.member_count,
            "qualifying_uniprot_count": item.labelled_protein_count,
            "assignment_stage": item.stage,
        }
        for item in (assignments[key] for key in sorted(assignments))
    ))
    _tsv(stage / "retained_clusters.tsv", [
        "mmseqs_cluster_id", "split", "uniref90_member_count",
        "qualifying_uniprot_count", "singleton", "giant_cluster",
    ], (
        {
            "mmseqs_cluster_id": cluster_id,
            "split": assignments[cluster_id].split,
            "uniref90_member_count": info.member_count,
            "qualifying_uniprot_count": info.labelled_protein_count,
            "singleton": int(info.member_count == 1),
            "giant_cluster": int(info.member_count >= giant_threshold),
        }
        for cluster_id, info in sorted(retained.items())
    ))

    retained_members = 0
    retained_unannotated = 0
    with _output_text(stage / "retained_cluster_members.tsv.gz") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mmseqs_cluster_id", "split", "uniref90_id", "sequence_sha256",
                "sequence_length", "connected_to_qualifying_uniprot",
            ],
            delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        for cluster_id, member_id, digest, length in (
            cluster_index.iter_assignments_with_metadata_for_clusters(uniref, retained_ids)
        ):
            retained_members += 1
            connected = member_id in annotated_uniref
            retained_unannotated += int(not connected)
            writer.writerow({
                "mmseqs_cluster_id": cluster_id,
                "split": assignments[cluster_id].split,
                "uniref90_id": member_id,
                "sequence_sha256": digest,
                "sequence_length": length,
                "connected_to_qualifying_uniprot": int(connected),
            })
    return retained_members, retained_unannotated


def _write_annotation_manifests(stage: Path, goa, decisions: list[MappingDecision]) -> None:
    by_raw, by_protein = _decision_index(decisions)
    _tsv(stage / "qualifying_annotations.tsv.gz", [
        "database", "raw_uniprot_accession", "uniprot_accession", "symbol",
        "raw_go_id", "canonical_go_id", "namespace", "aspect", "evidence_code", "qualifier",
        "reference", "with_from", "taxon_id", "assigned_date", "assigned_by",
        "annotation_extension", "gene_product_form", "gaf_line_number",
        "go_term_action", "accession_action", "uniref90_id", "mmseqs_cluster_id", "split",
    ], (
        {
            "database": record.database,
            "raw_uniprot_accession": record.raw_accession,
            "uniprot_accession": record.protein_id,
            "symbol": record.symbol,
            "raw_go_id": record.raw_go_id,
            "canonical_go_id": record.go_id,
            "namespace": record.namespace,
            "aspect": record.aspect,
            "evidence_code": record.evidence,
            "qualifier": record.qualifier,
            "reference": record.reference,
            "with_from": record.with_from,
            "taxon_id": record.taxon_id,
            "assigned_date": record.assigned_date,
            "assigned_by": record.assigned_by,
            "annotation_extension": record.annotation_extension,
            "gene_product_form": record.gene_product_form,
            "gaf_line_number": record.line_number,
            "go_term_action": record.term_action,
            "accession_action": record.accession_action,
            "uniref90_id": (by_raw.get(record.raw_accession) or by_protein.get(record.protein_id)).uniref90_id
            if (by_raw.get(record.raw_accession) or by_protein.get(record.protein_id)) else "",
            "mmseqs_cluster_id": (by_raw.get(record.raw_accession) or by_protein.get(record.protein_id)).mmseqs_cluster_id
            if (by_raw.get(record.raw_accession) or by_protein.get(record.protein_id)) else "",
            "split": (by_raw.get(record.raw_accession) or by_protein.get(record.protein_id)).split
            if (by_raw.get(record.raw_accession) or by_protein.get(record.protein_id)) else "",
        }
        for record in iter_annotation_records(goa)
    ))
    _tsv(stage / "excluded_annotations_sample.tsv.gz", [
        "database", "raw_uniprot_accession", "symbol", "raw_go_id", "evidence_code",
        "qualifier", "reference", "with_from", "aspect", "object_type", "taxon_id",
        "assigned_date", "assigned_by", "annotation_extension", "gene_product_form",
        "gaf_line_number", "rejection_reason", "detail",
    ], (
        {
            "database": item.database,
            "raw_uniprot_accession": item.raw_accession,
            "symbol": item.symbol,
            "raw_go_id": item.raw_go_id,
            "evidence_code": item.evidence,
            "qualifier": item.qualifier,
            "reference": item.reference,
            "with_from": item.with_from,
            "aspect": item.aspect,
            "object_type": item.object_type,
            "taxon_id": item.taxon_id,
            "assigned_date": item.assigned_date,
            "assigned_by": item.assigned_by,
            "annotation_extension": item.annotation_extension,
            "gene_product_form": item.gene_product_form,
            "gaf_line_number": item.line_number,
            "rejection_reason": item.rejection_reason,
            "detail": item.detail,
        }
        for item in iter_excluded_annotations(goa)
    ))


def _write_evidence_summary(stage: Path, goa) -> None:
    rejected = goa.rejected_evidence_reason_counts
    rows = [
        {"evidence_code": code, "disposition": "accepted", "reason": "qualifying", "rows": count}
        for code, count in sorted(goa.evidence_counts.items())
    ]
    rows.extend(
        {"evidence_code": code, "disposition": "rejected", "reason": reason, "rows": count}
        for (code, reason), count in sorted(rejected.items())
    )
    _tsv(stage / "evidence_summary.tsv", ["evidence_code", "disposition", "reason", "rows"], rows)
    _tsv(stage / "annotation_decision_counts.tsv", [
        "disposition", "reason", "evidence_code", "database", "object_type", "aspect", "rows",
    ], (
        {
            "disposition": key[0],
            "reason": key[1],
            "evidence_code": key[2],
            "database": key[3],
            "object_type": key[4],
            "aspect": key[5],
            "rows": count,
        }
        for key, count in sorted(goa.annotation_decision_counts.items())
    ))


def _write_attrition_summary(stage: Path, labels) -> None:
    payload = {
        "schema_version": 1,
        "row_population": "qualifying annotation rows entering the mapping/label chain",
        "row_denominator": labels.intended_annotation_rows,
        "row_terminal_counts": dict(sorted(labels.row_attrition_counts.items())),
        "protein_population": "unique pre-ontology-policy GOA raw protein accessions",
        "protein_denominator": labels.intended_accessions,
        "protein_terminal_counts": dict(sorted(labels.protein_attrition_counts.items())),
        "row_reconciles": sum(labels.row_attrition_counts.values()) == labels.intended_annotation_rows,
        "protein_reconciles": (
            sum(labels.protein_attrition_counts.values()) == labels.intended_accessions
        ),
        "retained_unannotated_members_policy": (
            "reported separately; never treated as supervised negatives or label-transfer recipients"
        ),
    }
    _json(stage / "attrition_summary.json", payload)
    _tsv(stage / "attrition_summary.tsv", ["unit", "terminal_status", "count", "denominator"], (
        {"unit": unit, "terminal_status": status, "count": count, "denominator": denominator}
        for unit, counts, denominator in (
            ("annotation_row", labels.row_attrition_counts, labels.intended_annotation_rows),
            ("protein_accession", labels.protein_attrition_counts, labels.intended_accessions),
        )
        for status, count in sorted(counts.items())
    ))


def _write_nonproduction_attrition_policy(
    path: Path,
    config: BuildConfig,
    repository_commit: str,
    frozen_manifest_sha256: str,
) -> None:
    """Write a non-authorizing measurement policy for fixtures and diagnostic pilots."""
    metrics = {}
    for name, definition in METRIC_DEFINITIONS.items():
        bound_key = f"allowed_{definition.bound}_ratio"
        metrics[name] = {
            "numerator_definition": definition.numerator,
            "denominator_definition": definition.denominator,
            bound_key: 0.0 if definition.bound == "minimum" else 1.0,
            "rationale": "Non-production measurement boundary; not a biological approval.",
            "evidence_source": "synthetic fixture or diagnostic pilot measurement only",
        }
    _json(path, {
        "schema_name": "homology-cluster-attrition-policy",
        "schema_version": 1,
        "uniprot_source_scope": config.uniprot_source_scope,
        "expected_releases": {
            "uniprot_uniref": config.release_uniprot,
            "goa": config.release_goa,
            "ontology": config.release_ontology,
        },
        "metrics": metrics,
        "rationale": "Measure software behavior without authorizing dissertation production.",
        "evidence_source": "non-production run",
        "author": "software-generated-non-production",
        "reviewer": "not-reviewed-for-production",
        "review_date": "2026-07-14",
        "framework_commit": repository_commit,
        "frozen_input_manifest_sha256": frozen_manifest_sha256,
    })


def _attrition_observations(
    config: BuildConfig,
    ontology: Ontology,
    labels,
    decisions: list[MappingDecision],
    assignments: dict,
    retained_members: int,
    uniref_count: int,
) -> dict[str, dict[str, object]]:
    selected = sum(item.canonical_sequence_available for item in decisions)
    mapped = sum(item.status == "mapped" for item in decisions)
    eligible_rows = int(labels.row_attrition_counts.get("eligible_annotation_row", 0))
    evaluable_raw = int(labels.protein_attrition_counts.get("evaluable_pfp", 0))
    intermediate = sum(len(labels.frames[split]) for split in SPLITS)
    evaluable = sum(
        bool(labels.restricted_annotations.get(str(row.proteins), ()))
        for split in SPLITS
        for row in labels.frames[split].itertuples(index=False)
    )
    namespace_prefix = {
        "bp": "biological_process",
        "cc": "cellular_component",
        "mf": "molecular_function",
    }
    ontology_evaluable = {}
    for prefix, namespace in namespace_prefix.items():
        terms = {
            term for term in labels.term_universe if ontology.namespace(term) == namespace
        }
        ontology_evaluable[prefix] = sum(
            bool(set(row.annotations) & terms)
            for split in SPLITS
            for row in labels.frames[split].itertuples(index=False)
        )
    balance = split_balance_metrics(assignments, config)
    return {
        "goa_to_selected_uniprot_mapping_ratio": observation(
            "goa_to_selected_uniprot_mapping_ratio", selected, len(decisions)
        ),
        "selected_uniprot_to_uniref90_mapping_ratio": observation(
            "selected_uniprot_to_uniref90_mapping_ratio", mapped, selected
        ),
        "qualifying_annotation_retention_ratio": observation(
            "qualifying_annotation_retention_ratio",
            eligible_rows,
            labels.intended_annotation_rows,
        ),
        "retained_cluster_member_ratio": observation(
            "retained_cluster_member_ratio", retained_members, uniref_count
        ),
        "evaluable_protein_ratio": observation(
            "evaluable_protein_ratio", evaluable_raw, labels.intended_accessions
        ),
        "propagated_term_evaluable_ratio": observation(
            "propagated_term_evaluable_ratio", evaluable, intermediate
        ),
        "bp_evaluable_ratio": observation(
            "bp_evaluable_ratio", ontology_evaluable["bp"], intermediate
        ),
        "cc_evaluable_ratio": observation(
            "cc_evaluable_ratio", ontology_evaluable["cc"], intermediate
        ),
        "mf_evaluable_ratio": observation(
            "mf_evaluable_ratio", ontology_evaluable["mf"], intermediate
        ),
        "development_split_deviation": observation(
            "development_split_deviation", balance["development_deviation"], 1
        ),
        "training_split_deviation": observation(
            "training_split_deviation", balance["training_deviation"], 1
        ),
    }


def _write_go_term_summary(stage: Path, goa, labels, ontology: Ontology) -> None:
    direct = goa.direct_go_counts
    unrestricted: dict[str, Counter] = {split: Counter() for split in SPLITS}
    restricted: dict[str, Counter] = {split: Counter() for split in SPLITS}
    for split in SPLITS:
        for row in labels.frames[split].itertuples(index=False):
            restricted[split].update(
                labels.restricted_annotations.get(str(row.proteins), ())
            )
            unrestricted[split].update(labels.unrestricted_annotations.get(str(row.proteins), ()))
    terms = sorted(set(direct) | set().union(*(set(counter) for counter in unrestricted.values())))
    _tsv(stage / "go_term_summary.tsv", [
        "go_id", "namespace", "direct_annotation_rows", "in_development_universe",
        "training_unrestricted_proteins", "validation_unrestricted_proteins", "test_unrestricted_proteins",
        "training_evaluable_proteins", "validation_evaluable_proteins", "test_evaluable_proteins",
    ], (
        {
            "go_id": term,
            "namespace": ontology.namespace(term),
            "direct_annotation_rows": direct[term],
            "in_development_universe": int(term in labels.term_universe),
            "training_unrestricted_proteins": unrestricted["training"][term],
            "validation_unrestricted_proteins": unrestricted["validation"][term],
            "test_unrestricted_proteins": unrestricted["test"][term],
            "training_evaluable_proteins": restricted["training"][term],
            "validation_evaluable_proteins": restricted["validation"][term],
            "test_evaluable_proteins": restricted["test"][term],
        }
        for term in terms
    ))


def _write_taxonomy_summary(stage: Path, goa, labels, catalog) -> None:
    annotation_taxa = goa.taxonomy_counts
    rows = [
        {"stage": "qualifying_annotation_rows", "split": "", "taxon_id": taxon, "count": count}
        for taxon, count in sorted(annotation_taxa.items())
    ]
    for split in SPLITS:
        taxa = Counter(
            (catalog.records.get(str(row.proteins)).taxon_id or "unknown")
            if catalog.records.get(str(row.proteins)) else "unknown"
            for row in labels.frames[split].itertuples(index=False)
        )
        rows.extend(
            {"stage": "supervised_proteins", "split": split, "taxon_id": taxon, "count": count}
            for taxon, count in sorted(taxa.items())
        )
    _tsv(stage / "taxonomy_summary.tsv", ["stage", "split", "taxon_id", "count"], rows)


def _write_split_summary(stage: Path, assignments: dict, labels) -> dict[str, dict[str, float | int]]:
    totals = {
        "clusters": len(assignments),
        "members": sum(item.member_count for item in assignments.values()),
        "labelled": sum(item.labelled_protein_count for item in assignments.values()),
    }
    summary: dict[str, dict[str, float | int]] = {}
    rows = []
    for split in SPLITS:
        selected = [item for item in assignments.values() if item.split == split]
        clusters = len(selected)
        members = sum(item.member_count for item in selected)
        labelled = sum(item.labelled_protein_count for item in selected)
        intermediate = len(labels.frames[split])
        evaluable = sum(
            bool(labels.restricted_annotations.get(str(row.proteins), ()))
            for row in labels.frames[split].itertuples(index=False)
        )
        summary[split] = {
            "clusters": clusters,
            "cluster_ratio": clusters / totals["clusters"] if totals["clusters"] else 0,
            "uniref90_members": members,
            "uniref90_member_ratio": members / totals["members"] if totals["members"] else 0,
            "qualifying_uniprot_mappings": labelled,
            "qualifying_uniprot_ratio": labelled / totals["labelled"] if totals["labelled"] else 0,
            "label_intermediate_proteins": intermediate,
            "evaluable_pfp_proteins": evaluable,
        }
        rows.append({"split": split, **summary[split]})
    _tsv(stage / "split_summary.tsv", [
        "split", "clusters", "cluster_ratio", "uniref90_members", "uniref90_member_ratio",
        "qualifying_uniprot_mappings", "qualifying_uniprot_ratio",
        "label_intermediate_proteins", "evaluable_pfp_proteins",
    ], rows)
    return summary


def _write_split_balance_summary(stage: Path, assignments: dict, config: BuildConfig) -> None:
    total_members = sum(item.member_count for item in assignments.values())
    development = [item for item in assignments.values() if item.split != "test"]
    development_members = sum(item.member_count for item in development)
    training_members = sum(
        item.member_count for item in assignments.values() if item.split == "training"
    )
    largest = max(assignments.values(), key=lambda item: (item.member_count, item.cluster_id))
    payload = {
        "schema_version": 1,
        "objective": (
            "uniref90_member_count" if config.split_policy == "sequence-balanced"
            else "cluster_count"
        ),
        "development_vs_test": {
            "requested_first_fraction": config.development_fraction,
            "achieved_first_fraction_by_members": development_members / total_members,
            "absolute_percentage_point_deviation_by_members": abs(
                development_members / total_members - config.development_fraction
            ) * 100,
            "denominator_members": total_members,
        },
        "training_vs_validation_within_development": {
            "requested_first_fraction": config.training_fraction_within_development,
            "achieved_first_fraction_by_members": training_members / development_members,
            "absolute_percentage_point_deviation_by_members": abs(
                training_members / development_members
                - config.training_fraction_within_development
            ) * 100,
            "denominator_members": development_members,
        },
        "largest_indivisible_cluster": {
            "cluster_id": largest.cluster_id,
            "member_count": largest.member_count,
            "share_of_retained_members": largest.member_count / total_members,
        },
        "optimization_boundary": (
            "deterministic bounded multi-candidate/local-swap heuristic; whole clusters remain "
            "indivisible and global subset-sum optimality is not claimed"
        ),
    }
    _json(stage / "split_balance_summary.json", payload)


def _percentile(values: list[int], fraction: float) -> int:
    if not values:
        return 0
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, max(0, math.ceil(fraction * len(ordered)) - 1))]


def _write_cluster_size_summary(stage: Path, assignments: dict, giant_threshold: int) -> None:
    rows = []
    for split in ("all", *SPLITS):
        sizes = [
            item.member_count for item in assignments.values()
            if split == "all" or item.split == split
        ]
        rows.append({
            "split": split,
            "cluster_count": len(sizes),
            "member_count": sum(sizes),
            "singleton_clusters": sum(value == 1 for value in sizes),
            "giant_clusters": sum(value >= giant_threshold for value in sizes),
            "minimum": min(sizes) if sizes else 0,
            "median": statistics.median(sizes) if sizes else 0,
            "mean": statistics.mean(sizes) if sizes else 0,
            "p95": _percentile(sizes, 0.95),
            "maximum": max(sizes) if sizes else 0,
        })
    _tsv(stage / "cluster_size_summary.tsv", [
        "split", "cluster_count", "member_count", "singleton_clusters", "giant_clusters",
        "minimum", "median", "mean", "p95", "maximum",
    ], rows)


def _write_benchmark_summary(stage: Path, summary: dict[str, object]) -> None:
    _json(stage / "benchmark_summary.json", summary)
    counts = summary["counts"]
    lines = [
        "# Homology-cluster benchmark summary", "",
        f"- Benchmark scope: `{summary['benchmark_scope']}`",
        f"- Selected UniProt source scope: `{summary['uniprot_source_scope']}`",
        f"- Identity threshold: **{summary['identity_percent']}%**",
        f"- Coverage: **{summary['coverage']}** using MMseqs2 cov-mode 0",
        f"- Split policy: `{summary['split_policy']}`",
        f"- Retained clusters: {counts['retained_mmseqs_clusters']}",
        f"- Retained UniRef90 entries: {counts['retained_uniref90_entries']}",
        f"- Evaluable PFP proteins: {counts['evaluable_pfp_proteins']}",
        f"- Development-defined terms (roots retained): {counts['development_defined_terms']}",
        "", "## Scientific boundary", "",
        "The software validates the recorded population, mappings, whole-cluster splits, schemas, "
        "and leakage constraints. It does not prove biological optimality, exhaustive low-identity "
        "edge discovery, or downstream model performance.", "",
    ]
    (stage / "benchmark_summary.md").write_text("\n".join(lines), encoding="utf-8")


def _validate_publication_marker(directory: Path) -> None:
    marker_path = directory / "RUN_COMPLETE.json"
    if not marker_path.is_file():
        raise ValueError(f"Completion marker is missing: {marker_path}")
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    if marker.get("complete") is not True:
        raise ValueError("Completion marker does not declare complete=true")
    if marker.get("post_publication_hash_verification") is not True:
        raise ValueError("Completion marker does not record final-path hash verification")
    if marker.get("manifest") != "output_manifest.json":
        raise ValueError("Completion marker must name literal output_manifest.json")
    if marker.get("publication_metadata") != "publication_metadata.json":
        raise ValueError("Completion marker must name literal publication_metadata.json")
    manifest_path = directory / "output_manifest.json"
    if sha256_file(manifest_path) != marker.get("manifest_sha256"):
        raise ValueError("Completion marker does not agree with output_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if marker.get("payload_file_count") != manifest.get("payload_file_count"):
        raise ValueError("Completion marker payload count does not agree with the manifest")
    verify_output_manifest(directory)
    publication_path = directory / "publication_metadata.json"
    if sha256_file(publication_path) != marker.get("publication_metadata_sha256"):
        raise ValueError("Completion marker does not agree with publication_metadata.json")
    publication = json.loads(publication_path.read_text(encoding="utf-8"))
    required_publication_keys = {
        "schema_version", "scientific_fingerprint_payload", *PUBLICATION_MARKER_KEYS,
    }
    missing_publication_keys = sorted(required_publication_keys - set(publication))
    if missing_publication_keys:
        raise ValueError(
            "Publication metadata is missing required keys: "
            + ", ".join(missing_publication_keys)
        )
    if publication.get("schema_version") != 1:
        raise ValueError("Publication metadata has an unsupported schema_version")
    if not isinstance(publication.get("fixture_mode"), bool) or not isinstance(
        publication.get("production_eligible"), bool
    ):
        raise ValueError("Publication fixture/eligibility fields must be booleans")
    for key in ("benchmark_scope", "split_policy", "training_population"):
        if not isinstance(publication.get(key), str) or not publication[key]:
            raise ValueError(f"Publication metadata field {key} must be a nonempty string")
    for key in ("seed", "min_count"):
        if (
            not isinstance(publication.get(key), int)
            or isinstance(publication.get(key), bool)
            or (key == "min_count" and publication[key] < 1)
        ):
            raise ValueError(f"Publication metadata field {key} must be a valid integer")
    for key in ("run_input_manifest_sha256", "frozen_input_manifest_sha256"):
        value = publication.get(key)
        if not isinstance(value, str) or not SHA256_PATTERN.fullmatch(value):
            raise ValueError(f"Publication metadata field {key} must be one lowercase SHA-256")
    executable_hash = publication.get("mmseqs_executable_sha256")
    if executable_hash is not None and (
        not isinstance(executable_hash, str)
        or not SHA256_PATTERN.fullmatch(executable_hash)
    ):
        raise ValueError("Publication MMseqs2 executable hash must be null or one SHA-256")
    for key in PUBLICATION_MARKER_KEYS:
        if marker.get(key) != publication.get(key):
            raise ValueError(f"Completion marker/publication metadata mismatch for {key}")
    if publication.get("benchmark_scope") not in {
        "dissertation-production", "diagnostic-pilot", "fixture-only",
        "all-thresholds-summary"
    }:
        raise ValueError("Publication metadata has an unknown benchmark_scope")
    if publication.get("benchmark_scope") == "all-thresholds-summary":
        expected_identities = [int(value * 100) for value in SUPPORTED_IDENTITIES]
        if (
            publication.get("identity_percent") is not None
            or publication.get("identities") != expected_identities
        ):
            raise ValueError("Aggregate publication identity scope is inconsistent")
    else:
        scope = publication.get("benchmark_scope")
        if (scope == "fixture-only") != publication["fixture_mode"]:
            raise ValueError("Publication benchmark_scope is inconsistent with fixture_mode")
        if scope in {"fixture-only", "diagnostic-pilot"} and publication.get(
            "production_eligible"
        ) is not False:
            raise ValueError("Non-production publication cannot be production eligible")
        if scope == "dissertation-production" and publication.get(
            "production_eligible"
        ) is not True:
            raise ValueError("Dissertation production publication lacks production authorization")
        if publication.get("identities") is not None or publication.get(
            "identity_percent"
        ) not in {30, 25, 20, 15, 10, 5}:
            raise ValueError("Child publication identity scope is inconsistent")
        if publication.get("benchmark_scope") in {
            "dissertation-production", "diagnostic-pilot"
        }:
            for key in (
                "expected_mmseqs_version", "observed_mmseqs_version",
                "mmseqs_resolved_executable", "repository_commit",
            ):
                if not isinstance(publication.get(key), str) or not publication[key].strip():
                    raise ValueError(
                        f"Production publication requires nonempty provenance field {key}"
                    )
    input_manifest_path = directory / "input_manifest.json"
    if sha256_file(input_manifest_path) != publication.get("run_input_manifest_sha256"):
        raise ValueError("Publication metadata does not bind input_manifest.json")
    frozen_manifest_path = directory / "frozen_input_manifest.json"
    if sha256_file(frozen_manifest_path) != publication.get("frozen_input_manifest_sha256"):
        raise ValueError("Publication metadata does not bind frozen_input_manifest.json")
    parameters = json.loads((directory / "parameters.json").read_text(encoding="utf-8"))
    for key in ("split_policy", "training_population", "seed", "min_count"):
        if publication.get(key) != parameters.get(key):
            raise ValueError(f"Publication metadata/parameters mismatch for {key}")
    if (
        publication.get("benchmark_scope") != "all-thresholds-summary"
        and publication.get("identity_percent") != parameters.get("identity_percent")
    ):
        raise ValueError("Publication metadata/parameters mismatch for identity_percent")
    input_manifest = json.loads(input_manifest_path.read_text(encoding="utf-8"))
    if publication.get("benchmark_scope") != "all-thresholds-summary":
        if publication.get("fixture_mode") != input_manifest.get("fixture_mode"):
            raise ValueError("Publication metadata/input manifest fixture-mode mismatch")
        if publication.get("benchmark_scope") != input_manifest.get("benchmark_scope"):
            raise ValueError("Publication metadata/input manifest benchmark-scope mismatch")
        if publication.get("uniprot_source_scope") != input_manifest.get(
            "uniprot_source_scope"
        ):
            raise ValueError("Publication metadata/input manifest source-scope mismatch")
        if (
            publication.get("frozen_input_manifest_sha256")
            != input_manifest.get("frozen_input_manifest", {}).get("sha256")
        ):
            raise ValueError("Publication metadata/input manifest frozen-manifest mismatch")
    provenance = json.loads((directory / "run_provenance.json").read_text(encoding="utf-8"))
    if publication.get("repository_commit") != provenance.get("repository", {}).get("commit"):
        raise ValueError("Publication metadata/repository commit mismatch")
    if publication.get("framework_revision") != publication.get("repository_commit"):
        raise ValueError("Publication requested/observed framework revision mismatch")
    runtime = provenance.get("runtime", {}).get("mmseqs2", {})
    if publication.get("expected_mmseqs_version") != runtime.get("expected_version"):
        raise ValueError("Publication metadata/MMseqs expected-version mismatch")
    if publication.get("observed_mmseqs_version") != runtime.get("observed_version_token"):
        raise ValueError("Publication metadata/MMseqs observed-version mismatch")
    if publication.get("mmseqs_resolved_executable") != runtime.get("resolved_executable"):
        raise ValueError("Publication metadata/MMseqs resolved-executable mismatch")
    if publication.get("mmseqs_executable_sha256") != runtime.get("executable_sha256"):
        raise ValueError("Publication metadata/MMseqs executable-hash mismatch")
    if publication["benchmark_scope"] in {
        "dissertation-production", "diagnostic-pilot"
    }:
        validate_recorded_exact_mmseqs_version(
            str(publication["expected_mmseqs_version"]),
            str(publication["observed_mmseqs_version"]),
        )
        if provenance.get("repository", {}).get("dirty") is not False:
            raise ValueError("Production publication provenance must record a clean repository")

    fingerprint = publication.get("scientific_fingerprint")
    payload = publication.get("scientific_fingerprint_payload")
    if not isinstance(fingerprint, str) or not SHA256_PATTERN.fullmatch(fingerprint):
        raise ValueError("Publication scientific_fingerprint must be one lowercase SHA-256")
    if not isinstance(payload, dict) or _scientific_fingerprint(payload) != fingerprint:
        raise ValueError("Publication scientific fingerprint payload/hash mismatch")
    frozen_manifest = load_frozen_input_manifest(
        frozen_manifest_path,
        uniprot_source_scope=str(publication["uniprot_source_scope"]),
        fixture_mode=publication["fixture_mode"],
    )
    expected_parameter_fields = (
        "split_policy", "development_fraction", "training_fraction_within_development",
        "training_population", "seed", "min_count", "evidence_codes",
        "include_relationships", "root_policy", "coverage", "cov_mode", "cluster_mode",
        "alignment_mode", "seq_id_mode", "createdb_shuffle", "cluster_reassign",
        "sensitivity", "evalue", "uniprot_release", "goa_release", "ontology_release",
        "uniprot_source_scope",
    )
    for key in expected_parameter_fields:
        if key not in parameters or payload.get(key) != parameters.get(key):
            raise ValueError(f"Scientific fingerprint/parameters mismatch for {key}")
    bound_values = {
        "frozen_input_source_fingerprint": frozen_manifest.source_fingerprint,
        "expected_mmseqs_version": publication.get("expected_mmseqs_version"),
        "observed_mmseqs_version": publication.get("observed_mmseqs_version"),
        "mmseqs_executable_sha256": publication.get("mmseqs_executable_sha256"),
        "repository_commit": publication.get("repository_commit"),
        "framework_revision": publication.get("framework_revision"),
        "uniprot_source_scope": publication.get("uniprot_source_scope"),
        "attrition_policy_sha256": publication.get("attrition_policy_sha256"),
        "requested_slots": publication.get("requested_slots"),
        "allocated_slots": publication.get("allocated_slots"),
        "mmseqs_threads": publication.get("mmseqs_threads"),
    }
    for key, expected in bound_values.items():
        if payload.get(key) != expected:
            raise ValueError(f"Scientific fingerprint binding mismatch for {key}")
    if publication.get("benchmark_scope") != "all-thresholds-summary":
        attrition_policy_path = directory / "attrition_policy.json"
        attrition_report_path = directory / "attrition_report.json"
        if sha256_file(attrition_policy_path) != publication.get("attrition_policy_sha256"):
            raise ValueError("Publication attrition-policy hash mismatch")
        if sha256_file(attrition_report_path) != publication.get("attrition_report_sha256"):
            raise ValueError("Publication attrition-report hash mismatch")
        attrition_report = json.loads(attrition_report_path.read_text(encoding="utf-8"))
        if attrition_report.get("uniprot_source_scope") != publication.get(
            "uniprot_source_scope"
        ):
            raise ValueError("Publication attrition report source-scope mismatch")
        if attrition_report.get("policy_sha256") != publication.get("attrition_policy_sha256"):
            raise ValueError("Publication attrition report policy binding mismatch")
        if attrition_report.get("input_manifest_sha256") != publication.get(
            "run_input_manifest_sha256"
        ):
            raise ValueError("Publication attrition report run-input-manifest binding mismatch")
        if publication.get("attrition_policy_passed") != attrition_report.get("policy_passed"):
            raise ValueError("Publication attrition policy outcome mismatch")
        if publication.get("attrition_override_valid") != attrition_report.get("override_valid"):
            raise ValueError("Publication attrition override outcome mismatch")
        expected_diagnostic = publication.get("benchmark_scope") == "diagnostic-pilot"
        if attrition_report.get("diagnostic") is not expected_diagnostic:
            raise ValueError("Publication attrition report diagnostic scope mismatch")
        if expected_diagnostic and attrition_report.get("production_authorized") is not False:
            raise ValueError("Diagnostic pilot attrition report cannot authorize production")
        if publication.get("production_eligible") and not attrition_report.get(
            "production_authorized"
        ):
            raise ValueError("Production publication lacks attrition authorization")


def _mapped_protein_sequence_rows(decisions, catalog):
    seen: set[str] = set()
    for decision in sorted(decisions, key=lambda item: (item.protein_id, item.raw_accession)):
        if (
            not decision.split
            or decision.protein_id in seen
            or decision.protein_id not in catalog.records
        ):
            continue
        seen.add(decision.protein_id)
        sequence = catalog.records[decision.protein_id].sequence
        yield (
            hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            decision.protein_id,
            decision.split,
        )


def _add_risk_warnings(report, goa, labels, decisions, catalog, config, ontology) -> None:
    rejected_term_rows = sum(
        count for key, count in goa.counters.items()
        if key in {
            "rejected_unknown", "rejected_obsolete_unresolved", "rejected_namespace_mismatch"
        }
    )
    if rejected_term_rows:
        report.add_warning(
            "go_canonicalization_loss",
            "Qualifying-policy rows were lost during GO canonicalization/namespace validation",
            rows=rejected_term_rows,
        )
    total = len(decisions)
    mapped = sum(item.status == "mapped" for item in decisions)
    if total and mapped / total < 0.90:
        report.add_warning(
            "low_mapping_coverage",
            "Less than 90% of qualifying GOA accessions completed the UniRef90/MMseqs mapping chain",
            mapped=mapped, total=total, ratio=mapped / total,
        )
    if labels.no_evaluable_term:
        report.add_warning(
            "proteins_without_evaluable_terms",
            "Some qualifying proteins cannot contribute evaluable PFP labels",
            counts=dict(labels.no_evaluable_term),
        )
    if labels.annotation_exclusion_counts:
        report.add_warning(
            "annotation_mapping_exclusions",
            "Some qualifying annotation rows were excluded because their own accession did not "
            "complete an unambiguous mapping chain",
            counts=dict(labels.annotation_exclusion_counts),
        )
    terms_by_namespace = {
        namespace: {term for term in labels.term_universe if ontology.namespace(term) == namespace}
        for namespace in PREFIX_TO_NAMESPACE.values()
    }
    ontology_counts = {
        namespace: sum(
            bool(set(row.annotations) & terms)
            for split in SPLITS for row in labels.frames[split].itertuples(index=False)
        )
        for namespace, terms in terms_by_namespace.items()
    }
    positive = [value for value in ontology_counts.values() if value]
    if len(positive) >= 2 and max(positive) / min(positive) >= 5:
        report.add_warning(
            "ontology_imbalance", "Supervised ontology populations differ by at least five-fold",
            counts=ontology_counts,
        )
    taxon_counts = Counter(
        catalog.records[str(row.proteins)].taxon_id or "unknown"
        for split in SPLITS for row in labels.frames[split].itertuples(index=False)
        if str(row.proteins) in catalog.records
    )
    if sum(taxon_counts.values()) >= 10 and taxon_counts:
        largest = max(taxon_counts.values()) / sum(taxon_counts.values())
        if largest >= 0.80:
            report.add_warning(
                "taxonomic_imbalance", "One taxon contributes at least 80% of supervised proteins",
                counts=dict(taxon_counts), largest_fraction=largest,
            )
    evidence_counts = {
        code: int(goa.evidence_counts.get(code, 0)) for code in sorted(config.evidence_codes)
    }
    positive_evidence = [count for count in evidence_counts.values() if count]
    absent_evidence = [code for code, count in evidence_counts.items() if count == 0]
    evidence_ratio = (
        max(positive_evidence) / min(positive_evidence)
        if len(positive_evidence) >= 2 else 1.0
    )
    if absent_evidence or evidence_ratio >= 10:
        report.add_warning(
            "evidence_code_imbalance",
            "One or more allowed evidence codes are absent or positive row counts differ by at least ten-fold",
            counts=evidence_counts, absent_codes=absent_evidence, positive_max_to_min=evidence_ratio,
        )
    development_support = Counter(
        term
        for split in ("training", "validation")
        for row in labels.frames[split].itertuples(index=False)
        for term in row.annotations
    )
    near_threshold = sorted(
        term for term in labels.term_universe
        if development_support[term] <= max(config.min_count, math.ceil(config.min_count * 1.10))
    )
    if near_threshold:
        report.add_warning(
            "low_effective_development_support",
            "Some retained terms are at or within 10% above the development min_count threshold",
            term_count=len(near_threshold), preview=near_threshold[:20],
        )


def validate_publication(directory: Path) -> None:
    _validate_publication_marker(directory.resolve())
    validation = json.loads((directory / "validation_report.json").read_text(encoding="utf-8"))
    if validation.get("valid") is not True:
        raise ValueError("Published validation_report.json does not declare valid=true")
    publication = json.loads(
        (directory / "publication_metadata.json").read_text(encoding="utf-8")
    )
    if publication.get("benchmark_scope") != "all-thresholds-summary":
        errors, _ = validate_term_universe_artifacts(directory, int(publication["min_count"]))
        if errors:
            raise ValueError(
                "Published terms.pkl fails independent development recount: " + "; ".join(errors)
            )


def build_benchmark(config: BuildConfig) -> BuildResult:
    config.validate()
    repository = Path(__file__).resolve().parents[4]
    repository_state = git_state(repository)
    if not config.fixture_mode and (
        repository_state["commit"] is None or repository_state["dirty"] is not False
    ):
        raise ValueError(
            "Production publication requires a clean, commit-addressable framework checkout; "
            f"observed repository state={repository_state}"
        )
    if not config.fixture_mode and repository_state["commit"] != config.framework_revision:
        raise ValueError(
            "Configured full framework revision does not match checked-out HEAD: "
            f"expected={config.framework_revision}, observed={repository_state['commit']}"
        )
    final_dir = config.output_dir.expanduser().resolve() / config.publication_relative_path
    work = config.temp_dir.expanduser().resolve() / (
        f"homology-{config.uniprot_source_scope}-{config.identity_directory}-"
        f"{config.split_policy}-{config.training_population}-seed-{config.seed}"
        f"-min-count-{config.min_count}-{config.run_id}"
    )
    if (
        work == final_dir
        or work.is_relative_to(final_dir)
        or final_dir.is_relative_to(work)
    ):
        raise ValueError(
            f"Scratch and final output paths must be disjoint: work={work}, final={final_dir}"
        )
    work.mkdir(parents=True, exist_ok=False)
    (work / "logs").mkdir()
    _json(work / "logs" / "run_started.json", {
        "started_at": utc_now(),
        "identity_percent": int(config.identity * 100),
        "fixture_mode": config.fixture_mode,
    })
    LOGGER.info(
        "Run started: identity=%d%% policy=%s fixture=%s scratch=%s final=%s",
        int(config.identity * 100), config.split_policy, config.fixture_mode, work, final_dir,
    )
    try:
        mmseqs_runtime = resolve_mmseqs_runtime(config.mmseqs_bin)
        if not config.fixture_mode:
            validate_exact_mmseqs_version(
                str(config.expected_mmseqs_version), mmseqs_runtime
            )
        elif config.cluster_assignments is None:
            if mmseqs_runtime.resolved_executable is None:
                raise FileNotFoundError(
                    f"MMseqs2 executable is unavailable: {config.mmseqs_bin}"
                )
            if mmseqs_runtime.observed_version is None:
                raise ValueError("Fixture MMseqs2 version output is empty")
            validate_mmseqs_version(mmseqs_runtime.observed_version)
        runtime_config = replace(
            config,
            mmseqs_bin=mmseqs_runtime.resolved_executable or config.mmseqs_bin,
        )
        initial_free = shutil.disk_usage(work).free
        if initial_free < config.minimum_free_disk_bytes:
            raise OSError(
                "Initial scratch preflight failed before input acquisition: "
                f"free={initial_free}, required={config.minimum_free_disk_bytes}"
            )
        expected_releases = {
            "uniprot_uniref": config.release_uniprot,
            "goa": config.release_goa,
            "ontology": config.release_ontology,
        }
        repository_commit = str(repository_state["commit"])
        policy_source = config.attrition_policy
        frozen_manifest = None
        policy = None
        attrition_policy_sha256 = None
        stage_started = time.monotonic()
        LOGGER.info("Stage started: validate frozen manifest and attrition policy contracts")
        if config.frozen_input_manifest is not None:
            frozen_manifest = load_frozen_input_manifest(
                config.frozen_input_manifest,
                uniprot_source_scope=config.uniprot_source_scope,
                fixture_mode=config.fixture_mode,
            )
            if policy_source is None:
                policy_source = work / "nonproduction_attrition_policy.json"
                _write_nonproduction_attrition_policy(
                    policy_source,
                    config,
                    repository_commit,
                    frozen_manifest.sha256,
                )
            policy, attrition_policy_sha256 = load_attrition_policy(
                policy_source,
                source_scope=config.uniprot_source_scope,
                expected_releases=expected_releases,
                framework_commit=repository_commit,
                frozen_input_manifest_sha256=frozen_manifest.sha256,
            )
        elif not config.fixture_mode:  # guarded by BuildConfig.validate
            raise ValueError("Production frozen-input manifest is unavailable")
        LOGGER.info(
            "Stage completed: validate frozen manifest and attrition policy contracts "
            "elapsed_seconds=%.1f",
            time.monotonic() - stage_started,
        )

        stage_started = time.monotonic()
        LOGGER.info("Stage started: resolve and hash frozen inputs")
        inputs = _resolved_inputs(config, work, frozen_manifest)
        if frozen_manifest is None:
            frozen_manifest = write_synthetic_fixture_manifest(
                work / "synthetic_frozen_input_manifest.json",
                _input_specs(config),
                inputs,
                config.uniprot_source_scope,
            )
            if policy_source is None:
                policy_source = work / "nonproduction_attrition_policy.json"
                _write_nonproduction_attrition_policy(
                    policy_source,
                    config,
                    repository_commit,
                    frozen_manifest.sha256,
                )
            policy, attrition_policy_sha256 = load_attrition_policy(
                policy_source,
                source_scope=config.uniprot_source_scope,
                expected_releases=expected_releases,
                framework_commit=repository_commit,
                frozen_input_manifest_sha256=frozen_manifest.sha256,
            )
        if policy is None or attrition_policy_sha256 is None or policy_source is None:
            raise RuntimeError("Attrition policy contract was not initialized")
        eligibility = bind_frozen_inputs(frozen_manifest, _input_specs(config), inputs)
        disk_preflight = _disk_preflight(config, work, inputs)
        loaded_common_cache = None
        LOGGER.info(
            "Stage completed: resolve and hash frozen inputs elapsed_seconds=%.1f",
            time.monotonic() - stage_started,
        )
        with staging_output(final_dir) as stage:
            parameters = _parameters(config)
            _json(stage / "parameters.json", parameters)
            _json(stage / "disk_preflight.json", disk_preflight)
            shutil.copyfile(frozen_manifest.path, stage / "frozen_input_manifest.json")

            ontology = Ontology(inputs["go_obo"].resolved_path, config.include_relationships)
            selected_sources = {
                {
                    "uniprot_sprot_sequences": "sprot",
                    "uniprot_trembl_sequences": "trembl",
                }[name]: inputs[name].resolved_path
                for name in config.selected_uniprot_input_names
            }
            if config.common_preprocessing_cache is None:
                stage_started = time.monotonic()
                LOGGER.info("Stage started: ontology and GOA parsing")
                goa = load_goa(
                    inputs["goa"].resolved_path,
                    ontology,
                    config.evidence_codes,
                    strict_malformed=config.strict_qc,
                    spool_dir=work / "goa_spool",
                    excluded_sample_per_reason=config.excluded_sample_per_reason,
                )
                _validate_frozen_metadata(config, ontology, goa.headers)
                _validate_manifest_embedded_metadata(frozen_manifest, ontology, goa.headers)
                LOGGER.info(
                    "Stage completed: ontology and GOA parsing qualifying_accessions=%d "
                    "elapsed_seconds=%.1f",
                    len(goa.qualifying_accessions), time.monotonic() - stage_started,
                )

                stage_started = time.monotonic()
                LOGGER.info(
                    "Stage started: UniRef90, UniProt, and accession mapping indexes"
                )
                uniref = UniRefIndex.build(
                    inputs["uniref90_fasta"].resolved_path, work / "uniref90.sqlite"
                )
                requested_raw = set(goa.qualifying_accessions or goa.annotations)
                catalog = load_requested_proteins_from_sources(
                    selected_sources,
                    requested_raw,
                    strict_collisions=not config.fixture_mode,
                    collision_database=work / "uniprot_accessions.sqlite",
                )
                goa = canonicalize_goa_accessions(goa, catalog)
                decisions = load_uniref90_mappings(
                    inputs["idmapping"].resolved_path, requested_raw, catalog, uniref
                )
                LOGGER.info(
                    "Stage completed: UniRef90, UniProt, and accession mapping indexes "
                    "uniref_entries=%d mapping_decisions=%d elapsed_seconds=%.1f",
                    uniref.count(), len(decisions), time.monotonic() - stage_started,
                )
            else:
                stage_started = time.monotonic()
                LOGGER.info("Stage started: load common preprocessing cache")
                loaded_common_cache = load_common_preprocessing_cache(
                    config.common_preprocessing_cache,
                    source_scope=config.uniprot_source_scope,
                    expected_input_sha256={
                        name: str(entry["sha256"])
                        for name, entry in frozen_manifest.entries.items()
                    },
                    include_relationships=config.include_relationships,
                    strict_qc=config.strict_qc,
                    excluded_sample_per_reason=config.excluded_sample_per_reason,
                    frozen_input_manifest_sha256=frozen_manifest.sha256,
                )
                goa = loaded_common_cache.goa
                catalog = loaded_common_cache.catalog
                decisions = loaded_common_cache.decisions
                requested_raw = loaded_common_cache.requested_raw
                uniref = loaded_common_cache.uniref
                _validate_frozen_metadata(config, ontology, goa.headers)
                _validate_manifest_embedded_metadata(frozen_manifest, ontology, goa.headers)
                shutil.copyfile(
                    loaded_common_cache.root / CACHE_MARKER,
                    stage / "common_preprocessing_cache_manifest.json",
                )
                LOGGER.info(
                    "Stage completed: load common preprocessing cache "
                    "uniref_entries=%d mapping_decisions=%d elapsed_seconds=%.1f",
                    uniref.count(), len(decisions), time.monotonic() - stage_started,
                )

            stage_started = time.monotonic()
            LOGGER.info("Stage started: MMseqs2 execution or fixture assignment validation")
            mmseqs_work = work / "mmseqs"
            mmseqs_work.mkdir()
            commands = build_mmseqs_commands(
                runtime_config, inputs["uniref90_fasta"].resolved_path, mmseqs_work
            )
            write_command_manifest(stage / "mmseqs_commands.tsv", commands)
            write_command_manifest(work / "logs" / "mmseqs_commands.tsv", commands)
            if config.cluster_assignments:
                source_clusters = config.cluster_assignments.expanduser().resolve()
                if not source_clusters.is_file():
                    raise FileNotFoundError(source_clusters)
                fixture_log_dir = stage / "logs" / "mmseqs"
                fixture_log_dir.mkdir(parents=True)
                _json(fixture_log_dir / "NOT_EXECUTED.json", {
                    "fixture_mode": True,
                    "benchmark_scope": "fixture-only",
                    "mmseqs_execution": "not executed by this run",
                    "cluster_assignments": str(source_clusters),
                    "cluster_assignments_sha256": sha256_file(source_clusters),
                    "runtime_mmseqs_probe": mmseqs_runtime.as_dict(
                        config.expected_mmseqs_version
                    ),
                })
            else:
                if mmseqs_runtime.resolved_executable is None:
                    raise FileNotFoundError(
                        f"MMseqs2 executable is unavailable: {config.mmseqs_bin}; "
                        "supply --mmseqs-bin or use --cluster-assignments only for a validated fixture"
                    )
                execute_commands(commands, work / "logs" / "mmseqs")
                shutil.copytree(work / "logs" / "mmseqs", stage / "logs" / "mmseqs")
                source_clusters = mmseqs_work / "uniref90_clusters.tsv"
            cluster_index = ClusterIndex.build(
                source_clusters, uniref, work / "clusters.sqlite"
            )
            LOGGER.info(
                "Stage completed: MMseqs2 execution or fixture assignment validation "
                "clusters=%d members=%d elapsed_seconds=%.1f",
                cluster_index.cluster_count(), cluster_index.member_count(),
                time.monotonic() - stage_started,
            )

            stage_started = time.monotonic()
            LOGGER.info("Stage started: cluster retention, splitting, labels, and PFP exports")
            decisions = connect_proteins_to_clusters(decisions, cluster_index)
            retained = retained_cluster_info(decisions, cluster_index)
            if len(retained) < 3:
                raise ValueError(
                    f"Only {len(retained)} retained clusters contain qualifying annotations; "
                    "at least three are required for train/validation/test"
                )
            development_test_assignments = assign_development_test(
                retained,
                config.split_policy,
                config.seed,
                config.development_fraction,
            )
            decisions = [
                replace(
                    item,
                    split=development_test_assignments[item.mmseqs_cluster_id].split,
                )
                if item.mmseqs_cluster_id in development_test_assignments else item
                for item in decisions
            ]
            labels = build_labels(
                ontology,
                goa,
                catalog,
                decisions,
                development_test_assignments,
                config.min_count,
                finalize_development=lambda: assign_training_validation(
                    retained,
                    development_test_assignments,
                    config.split_policy,
                    config.seed,
                    config.training_fraction_within_development,
                ),
            )
            assignments = labels.cluster_assignments
            decisions = [
                replace(item, split=assignments[item.mmseqs_cluster_id].split)
                if item.mmseqs_cluster_id in assignments else item
                for item in decisions
            ]
            export_all(labels, ontology, stage)
            LOGGER.info(
                "Stage completed: cluster retention, splitting, labels, and PFP exports "
                "retained_clusters=%d terms=%d elapsed_seconds=%.1f",
                len(retained), len(labels.term_universe), time.monotonic() - stage_started,
            )

            stage_started = time.monotonic()
            LOGGER.info("Stage started: audit manifests and summaries")
            _write_mapping_manifests(stage, decisions)
            retained_members, retained_unannotated = _write_cluster_manifests(
                stage, cluster_index, uniref, retained, assignments, decisions,
                config.giant_cluster_threshold,
            )
            _write_annotation_manifests(stage, goa, decisions)
            _write_evidence_summary(stage, goa)
            _write_attrition_summary(stage, labels)
            _write_go_term_summary(stage, goa, labels, ontology)
            _write_taxonomy_summary(stage, goa, labels, catalog)
            split_summary = _write_split_summary(stage, assignments, labels)
            _write_split_balance_summary(stage, assignments, config)
            _write_cluster_size_summary(stage, assignments, config.giant_cluster_threshold)

            mapping_summary = {
                "uniprot_source_scope": config.uniprot_source_scope,
                "qualifying_goa_accessions": len(requested_raw),
                "canonical_annotated_uniprot_proteins": len(goa.qualifying_accessions),
                "loaded_canonical_sequences": len(catalog.records),
                "ambiguous_sequence_aliases": sorted(catalog.ambiguous_aliases),
                "source_counts": catalog.source_counts,
                "collision_counts": dict(sorted(catalog.collision_counts.items())),
                "mapping_status_counts": mapping_counters(decisions),
                "mapping_counts_by_source": mapping_counters_by_source(
                    decisions, set(selected_sources), catalog.source_counts
                ),
            }
            _json(stage / "mapping_summary.json", mapping_summary)

            structural_input_eligible = not config.fixture_mode and all(eligibility.values())
            input_manifest = {
                "schema_version": 3,
                "fixture_mode": config.fixture_mode,
                "benchmark_scope": config.benchmark_scope,
                "uniprot_source_scope": config.uniprot_source_scope,
                "structural_input_eligible": structural_input_eligible,
                "eligibility": eligibility,
                "frozen_input_manifest": {
                    "published_path": "frozen_input_manifest.json",
                    "sha256": frozen_manifest.sha256,
                    "source_fingerprint": frozen_manifest.source_fingerprint,
                },
                "frozen_endpoint": {
                    "uniprot_uniref_release": config.release_uniprot,
                    "goa_release": config.release_goa,
                    "ontology_release": config.release_ontology,
                },
                "inputs": {name: item.as_dict() for name, item in sorted(inputs.items())},
                "common_preprocessing_cache": (
                    {
                        "used": True,
                        "published_manifest": "common_preprocessing_cache_manifest.json",
                        "marker_sha256": loaded_common_cache.marker_sha256,
                        "source_scope": config.uniprot_source_scope,
                        "cached_stages": [
                            "GOA parsing",
                            "UniRef90 sequence index",
                            "selected UniProt catalogue",
                            "GOA accession canonicalization",
                            "UniProt-to-UniRef90 mapping",
                        ],
                        "threshold_specific_data_cached": False,
                    }
                    if loaded_common_cache is not None else {"used": False}
                ),
                "selected_uniprot_population": {
                    "scope": config.uniprot_source_scope,
                    "required_roles": list(config.selected_uniprot_input_names),
                    "source_counts": catalog.source_counts,
                    "collision_counts": dict(sorted(catalog.collision_counts.items())),
                    "production_sequence_format": "DAT",
                    "fasta_limitation": (
                        "FASTA diagnostic fixtures cannot provide the secondary-accession "
                        "information available from authoritative UniProt DAT"
                    ),
                },
                "scratch_preflight": disk_preflight,
                "mmseqs_runtime": mmseqs_runtime.as_dict(config.expected_mmseqs_version),
                "embedded_metadata": {
                    "goa_headers": goa.headers,
                    "ontology_data_version": ontology.data_version,
                    "ontology_relationship_types_observed": sorted(ontology.relationship_types),
                    "uniref90_release_validation": (
                        "No release field is embedded in FASTA; configured release plus SHA-256 binds the input"
                    ),
                    "idmapping_schema": "headerless idmapping_selected; exactly 22 columns; UniRef90 column 9",
                },
            }
            if config.cluster_assignments:
                input_manifest["synthetic_or_external_cluster_assignments"] = {
                    "resolved_path": str(source_clusters),
                    "size_bytes": source_clusters.stat().st_size,
                    "sha256": sha256_file(source_clusters),
                    "mmseqs_execution": "not executed by this run",
                }
            _json(stage / "input_manifest.json", input_manifest)

            if sha256_file(policy_source) != attrition_policy_sha256:
                raise ValueError("Attrition policy changed before evaluation")
            shutil.copyfile(policy_source, stage / "attrition_policy.json")
            observations = _attrition_observations(
                config,
                ontology,
                labels,
                decisions,
                assignments,
                retained_members,
                uniref.count(),
            )
            attrition_report = evaluate_attrition(
                policy,
                attrition_policy_sha256,
                observations,
                source_scope=config.uniprot_source_scope,
                framework_commit=repository_commit,
                input_manifest_sha256=sha256_file(stage / "input_manifest.json"),
                override_path=config.attrition_override,
                diagnostic=config.diagnostic_pilot,
            )
            _json(stage / "attrition_report.json", attrition_report)
            if config.attrition_override is not None:
                shutil.copyfile(config.attrition_override, stage / "attrition_override.json")
            production_eligible = (
                config.benchmark_scope == "dissertation-production"
                and structural_input_eligible
                and attrition_report["production_authorized"] is True
            )

            label_intermediate_count = sum(len(labels.frames[split]) for split in SPLITS)
            evaluable_pfp_count = sum(
                bool(labels.restricted_annotations.get(str(row.proteins), ()))
                for split in SPLITS
                for row in labels.frames[split].itertuples(index=False)
            )
            summary = {
                "schema_version": 1,
                "fixture_mode": config.fixture_mode,
                "production_eligible": production_eligible,
                "benchmark_scope": config.benchmark_scope,
                "uniprot_source_scope": config.uniprot_source_scope,
                "framework_revision": config.framework_revision,
                "run_id": config.run_id,
                "identity_percent": int(config.identity * 100),
                "identity_fraction": config.identity,
                "coverage": config.coverage,
                "split_policy": config.split_policy,
                "training_population": config.training_population,
                "split_summary": split_summary,
                "counts": {
                    "total_uniref90_entries": uniref.count(),
                    "total_mmseqs_clusters": cluster_index.cluster_count(),
                    "qualifying_goa_accessions": len(requested_raw),
                    "qualifying_mapped_to_uniref90": sum(item.status == "mapped" for item in decisions),
                    "qualifying_unmapped_or_ambiguous": sum(item.status != "mapped" for item in decisions),
                    "retained_mmseqs_clusters": len(retained),
                    "retained_uniref90_entries": retained_members,
                    "retained_annotated_uniprot_proteins": len({
                        item.protein_id for item in decisions if item.mmseqs_cluster_id in retained
                    }),
                    "retained_unannotated_uniref90_entries": retained_unannotated,
                    "label_intermediate_proteins": label_intermediate_count,
                    "evaluable_pfp_proteins": evaluable_pfp_count,
                    "development_defined_terms": len(labels.term_universe),
                    "singleton_retained_clusters": sum(info.member_count == 1 for info in retained.values()),
                    "giant_retained_clusters": sum(
                        info.member_count >= config.giant_cluster_threshold for info in retained.values()
                    ),
                },
                "labels_removed_outside_development_universe": dict(labels.removed_term_counts),
                "annotation_rows_excluded_by_mapping_chain": dict(
                    labels.annotation_exclusion_counts
                ),
                "proteins_without_evaluable_terms": dict(labels.no_evaluable_term),
                "attrition_gate": {
                    "policy_sha256": attrition_policy_sha256,
                    "policy_passed": attrition_report["policy_passed"],
                    "override_valid": attrition_report["override_valid"],
                    "production_authorized": attrition_report["production_authorized"],
                    "failed_metrics": attrition_report["failed_metrics"],
                },
                "limitations": [
                    "MMseqs2 cluster mode 0 is greedy set cover, not an exhaustive equivalence relation.",
                    "At 5-15% identity, prefilter recall, E-value, database size, and MMseqs2 version strongly affect results.",
                    "Internal validation does not prove biological optimality or downstream PFP performance.",
                ],
            }
            _write_benchmark_summary(stage, summary)

            _json(
                stage / "run_provenance.json",
                runtime_provenance(
                    repository,
                    inputs,
                    parameters,
                    mmseqs_runtime.as_dict(config.expected_mmseqs_version),
                ),
            )
            fingerprint_payload = _scientific_fingerprint_payload(
                config,
                frozen_manifest,
                mmseqs_runtime,
                repository_state["commit"],
                attrition_policy_sha256,
            )
            scientific_fingerprint = _scientific_fingerprint(fingerprint_payload)
            publication_metadata = {
                "schema_version": 1,
                "fixture_mode": config.fixture_mode,
                "production_eligible": production_eligible,
                "benchmark_scope": config.benchmark_scope,
                "uniprot_source_scope": config.uniprot_source_scope,
                "framework_revision": config.framework_revision or repository_commit,
                "run_id": config.run_id,
                "identity_percent": int(config.identity * 100),
                "identities": None,
                "split_policy": config.split_policy,
                "training_population": config.training_population,
                "seed": config.seed,
                "min_count": config.min_count,
                "run_input_manifest_sha256": sha256_file(stage / "input_manifest.json"),
                "frozen_input_manifest_sha256": sha256_file(
                    stage / "frozen_input_manifest.json"
                ),
                "attrition_policy_sha256": sha256_file(stage / "attrition_policy.json"),
                "attrition_report_sha256": sha256_file(stage / "attrition_report.json"),
                "attrition_override_sha256": attrition_report["override_sha256"],
                "attrition_policy_passed": attrition_report["policy_passed"],
                "attrition_override_valid": attrition_report["override_valid"],
                "requested_slots": config.requested_slots,
                "allocated_slots": config.allocated_slots,
                "mmseqs_threads": config.threads,
                "expected_mmseqs_version": config.expected_mmseqs_version,
                "observed_mmseqs_version": mmseqs_runtime.version_token,
                "mmseqs_resolved_executable": mmseqs_runtime.resolved_executable,
                "mmseqs_executable_sha256": mmseqs_runtime.executable_sha256,
                "repository_commit": repository_state["commit"],
                "scientific_fingerprint_payload": fingerprint_payload,
                "scientific_fingerprint": scientific_fingerprint,
            }
            write_publication_metadata(stage, publication_metadata)
            LOGGER.info(
                "Stage completed: audit manifests and summaries elapsed_seconds=%.1f",
                time.monotonic() - stage_started,
            )

            stage_started = time.monotonic()
            LOGGER.info("Stage started: strict semantic validation")
            cluster_splits = {key: item.split for key, item in assignments.items()}
            validation = validate_outputs(
                stage, config, ontology, goa, labels, decisions, assignments, commands,
                cluster_index, uniref, uniref.count(), cluster_index.member_count(),
                cluster_index.global_exact_sequence_split_conflicts(
                    uniref, cluster_splits, _mapped_protein_sequence_rows(decisions, catalog)
                ),
            )
            _add_risk_warnings(
                validation, goa, labels, decisions, catalog, config, ontology
            )
            write_validation_report(validation, stage)
            if not validation.valid:
                shutil.copyfile(
                    stage / "validation_report.json", work / "logs" / "validation_report.json"
                )
                shutil.copyfile(
                    stage / "validation_report.md", work / "logs" / "validation_report.md"
                )
                failed = [item["name"] for item in validation.checks if not item["passed"]]
                raise ValueError("Benchmark validation failed: " + ", ".join(failed))
            if (
                config.benchmark_scope == "dissertation-production"
                and not attrition_report["production_authorized"]
            ):
                shutil.copyfile(
                    stage / "attrition_report.json", work / "logs" / "attrition_report.json"
                )
                failed_metrics = [
                    item["name"] for item in attrition_report["failed_metrics"]
                ]
                raise ValueError(
                    "Production attrition policy failed without a valid reviewed override: "
                    + ", ".join(failed_metrics)
                )
            LOGGER.info(
                "Stage completed: strict semantic validation checks=%d warnings=%d "
                "elapsed_seconds=%.1f",
                len(validation.checks), len(validation.warnings),
                time.monotonic() - stage_started,
            )

            stage_started = time.monotonic()
            LOGGER.info("Stage started: final input recheck and atomic publication")
            _verify_inputs_unchanged(inputs)
            if loaded_common_cache is not None:
                refreshed_cache = inspect_common_preprocessing_cache(
                    loaded_common_cache.root,
                    expected_source_scope=config.uniprot_source_scope,
                    expected_input_sha256={
                        name: str(entry["sha256"])
                        for name, entry in frozen_manifest.entries.items()
                    },
                    verify_file_hashes=True,
                )
                if sha256_file(loaded_common_cache.root / CACHE_MARKER) != (
                    loaded_common_cache.marker_sha256
                ):
                    raise ValueError(
                        "Common preprocessing cache marker changed while the run was in progress"
                    )
                if refreshed_cache != loaded_common_cache.payload:
                    raise ValueError(
                        "Common preprocessing cache contract changed while the run was in progress"
                    )
            verify_frozen_manifest_unchanged(frozen_manifest)
            if sha256_file(policy_source) != attrition_policy_sha256:
                raise ValueError("Attrition policy changed while the run was in progress")
            if config.attrition_override is not None and (
                sha256_file(config.attrition_override)
                != attrition_report["override_sha256"]
            ):
                raise ValueError("Attrition override changed while the run was in progress")
            verify_mmseqs_executable_unchanged(mmseqs_runtime)
            if not config.fixture_mode:
                final_repository_state = git_state(repository)
                if (
                    final_repository_state.get("commit") != repository_state.get("commit")
                    or final_repository_state.get("dirty") is not False
                ):
                    raise ValueError(
                        "Framework repository changed or became dirty while the production run "
                        f"was in progress: initial={repository_state}, final={final_repository_state}"
                    )
            if config.cluster_assignments:
                fixture_record = json.loads(
                    (stage / "logs" / "mmseqs" / "NOT_EXECUTED.json").read_text(encoding="utf-8")
                )
                expected_fixture_hash = fixture_record["cluster_assignments_sha256"]
                if sha256_file(source_clusters) != expected_fixture_hash:
                    raise ValueError("Fixture cluster assignments changed while the run was in progress")
            write_output_manifest(stage)
            publish(stage, final_dir)
            try:
                verify_output_manifest(final_dir)
                write_completion_marker(final_dir)
                validate_publication(final_dir)
            except BaseException:
                shutil.rmtree(final_dir, ignore_errors=True)
                raise
            LOGGER.info(
                "Stage completed: final input recheck and atomic publication elapsed_seconds=%.1f",
                time.monotonic() - stage_started,
            )

        files = tuple(sorted((path for path in final_dir.rglob("*") if path.is_file())))
        LOGGER.info("Run completed successfully: files=%d final=%s", len(files), final_dir)
        return BuildResult(final_dir, config.identity, files, True)
    except BaseException as exc:
        preserved = _preserve_failure_logs(config, work, exc)
        if preserved is not None:
            LOGGER.error("Preserved failure diagnostics at %s", preserved)
        raise
    finally:
        if not config.keep_temp:
            shutil.rmtree(work, ignore_errors=True)
