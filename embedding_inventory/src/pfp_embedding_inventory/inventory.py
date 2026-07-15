from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .benchmark import sequence_index, temporal_text_role
from .models import (
    AliasEntry,
    ArrayInfo,
    ArtifactVerification,
    BenchmarkData,
    InventoryRecord,
    InventoryResult,
    MODALITIES,
    ModalitySpec,
    PlannerConfig,
    ProteinRecord,
)
from .paths import PathSafetyError, require_resolved_within, resolve_within


class InventoryError(ValueError):
    pass


ArrayCache = Dict[Tuple[str, int], ArrayInfo]


@dataclass(frozen=True)
class Candidate:
    source_id: str
    route: str
    source_file: Path
    ambiguity: str = ""
    note: str = ""
    alias_source_identity: str = ""
    alias_mapping_evidence: str = ""


def build_inventory(
    benchmark: BenchmarkData,
    source_benchmark: BenchmarkData,
    embedding_cache: Path,
    config: PlannerConfig,
    policy: str,
    aliases: Optional[Dict[Tuple[str, str], List[AliasEntry]]] = None,
    array_cache: Optional[ArrayCache] = None,
    artifact_verification: Optional[ArtifactVerification] = None,
) -> InventoryResult:
    if policy not in {"paper-faithful", "maximize-coverage"}:
        raise InventoryError("Unsupported action policy: %s" % policy)
    aliases = aliases or {}
    artifact_verification = artifact_verification or ArtifactVerification(
        configured=False,
        verified=False,
        artifact_id="",
        checks={},
        reasons=["artifact verification was not performed"],
        expected={},
        observed={},
    )
    _validate_artifact_verification_binding(
        artifact_verification, embedding_cache, source_benchmark, config
    )
    array_cache = array_cache if array_cache is not None else {}
    embedding_cache = embedding_cache.resolve()
    unknown_alias_targets = sorted(
        {protein_id for protein_id, _ in aliases if protein_id not in benchmark.proteins}
    )
    if unknown_alias_targets:
        raise InventoryError(
            "Alias mappings target proteins absent from the benchmark: %s"
            % ", ".join(unknown_alias_targets[:5])
        )

    cache_ids: Dict[str, Set[str]] = {}
    for modality in MODALITIES:
        try:
            source_dir = resolve_within(
                embedding_cache,
                Path(config.modalities[modality].directory),
                "%s modality directory" % modality,
            )
        except PathSafetyError as exc:
            raise InventoryError(str(exc)) from exc
        cache_ids[modality] = {
            path.stem for path in source_dir.glob("*.npy") if path.is_file()
        } if source_dir.is_dir() else set()

    source_by_sequence = sequence_index(source_benchmark)
    records: List[InventoryRecord] = []
    used_source_ids: Dict[str, Set[str]] = {modality: set() for modality in MODALITIES}
    for protein_id in sorted(benchmark.proteins):
        target = benchmark.proteins[protein_id]
        for modality in MODALITIES:
            spec = config.modalities[modality]
            try:
                source_dir = resolve_within(
                    embedding_cache, Path(spec.directory),
                    "%s modality directory" % modality,
                )
            except PathSafetyError as exc:
                raise InventoryError(str(exc)) from exc
            record = _inventory_one(
                target=target,
                modality=modality,
                spec=spec,
                source_dir=source_dir,
                source_benchmark=source_benchmark,
                source_by_sequence=source_by_sequence,
                alias_entries=aliases.get((protein_id, modality), []),
                array_cache=array_cache,
                policy=policy,
                artifact_verification=artifact_verification,
                target_matches_artifact=(
                    benchmark.fingerprint
                    == config.artifact_scope.expected_benchmark_fingerprint
                ),
            )
            records.append(record)
            if record.requested_action == "reuse" and record.source_protein_id:
                used_source_ids[modality].add(record.source_protein_id)

    return InventoryResult(
        benchmark=benchmark,
        source_benchmark=source_benchmark,
        records=records,
        cache_ids=cache_ids,
        used_source_ids=used_source_ids,
        policy=policy,
        config=config,
        artifact_verification=artifact_verification,
    )


def _validate_artifact_verification_binding(
    verification: ArtifactVerification,
    embedding_cache: Path,
    source_benchmark: BenchmarkData,
    config: PlannerConfig,
) -> None:
    """Refuse a successful proof for a different published cache or source."""
    if not verification.verified:
        return
    observed_catalog = verification.observed.get("cache_catalog", {})
    bound = (
        config.artifact_scope.mode == "verified-published-cache"
        and verification.artifact_id == config.artifact_scope.artifact_id
        and bool(verification.checks)
        and all(verification.checks.values())
        and verification.observed.get("source_benchmark_fingerprint")
        == source_benchmark.fingerprint
        and verification.observed.get("embedding_cache_root")
        == str(embedding_cache.resolve())
        and observed_catalog.get("fingerprint")
        == config.artifact_scope.expected_cache_catalog_fingerprint
    )
    if not bound:
        raise InventoryError(
            "successful artifact verification is not bound to this source benchmark, "
            "cache catalog, and published-artifact configuration"
        )


def validate_array(path: Path, expected_dim: int, cache: ArrayCache) -> ArrayInfo:
    key = (str(path), expected_dim)
    cached = cache.get(key)
    if cached is not None:
        return cached
    if not path.exists():
        info = ArrayInfo(exists=False)
        cache[key] = info
        return info
    if not path.is_file():
        info = ArrayInfo(exists=True, error="path is not a regular file")
        cache[key] = info
        return info
    try:
        require_resolved_within(path.parent, path, "embedding array")
    except PathSafetyError as exc:
        info = ArrayInfo(exists=True, error=str(exc))
        cache[key] = info
        return info
    try:
        array = np.load(path, mmap_mode="r", allow_pickle=False)
        shape = str(tuple(array.shape))
        dtype = str(array.dtype)
        dtype_supported = bool(np.issubdtype(array.dtype, np.floating))
        if not dtype_supported:
            info = ArrayInfo(
                exists=True,
                observed_shape=shape,
                dtype=dtype,
                dtype_supported=False,
            )
        else:
            finite = bool(np.isfinite(array).all())
            valid = tuple(array.shape) == (expected_dim,) and finite
            info = ArrayInfo(
                exists=True,
                observed_shape=shape,
                dtype=dtype,
                finite=finite,
                dtype_supported=True,
                valid=valid,
            )
    except Exception as exc:  # report every NumPy/filesystem decoding failure
        info = ArrayInfo(exists=True, error="%s: %s" % (type(exc).__name__, exc))
    cache[key] = info
    return info


def _inventory_one(
    target: ProteinRecord,
    modality: str,
    spec: ModalitySpec,
    source_dir: Path,
    source_benchmark: BenchmarkData,
    source_by_sequence: Dict[str, List[str]],
    alias_entries: List[AliasEntry],
    array_cache: ArrayCache,
    policy: str,
    artifact_verification: ArtifactVerification,
    target_matches_artifact: bool,
) -> InventoryRecord:
    candidate = _select_candidate(
        target,
        modality,
        spec,
        source_dir,
        source_benchmark,
        source_by_sequence,
        alias_entries,
        array_cache,
        artifact_verification.verified,
    )
    if candidate.ambiguity:
        candidate_paths = sorted(
            {source_dir / (entry.source_protein_id + ".npy") for entry in alias_entries}
        )
        direct_path = source_dir / (target.protein_id + ".npy")
        paths = sorted(set(candidate_paths + [direct_path]))
        infos = [(path, validate_array(path, spec.expected_dim, array_cache)) for path in paths]
        existing = [(path, info) for path, info in infos if info.exists]
        exists = bool(existing)
        finite_values = [info.finite for _, info in existing if info.finite is not None]
        finite = all(finite_values) if finite_values else None
        valid = bool(existing) and all(info.valid for _, info in existing)
        observed_shapes = sorted({info.observed_shape for _, info in existing if info.observed_shape})
        dtypes = sorted({info.dtype for _, info in existing if info.dtype})
        source_files = [str(Path(spec.directory) / path.name) for path, _ in existing]
        status = "provenance-unknown"
        action = "regenerate"
        reason = "ambiguous explicit alias mapping: %s" % candidate.ambiguity
        return InventoryRecord(
            protein_id=target.protein_id,
            sequence_sha256=target.sequence_sha256,
            modality=modality,
            source_directory=spec.directory,
            source_file=";".join(source_files),
            exists=exists,
            observed_shape=";".join(observed_shapes),
            expected_shape=str((spec.expected_dim,)),
            dtype=";".join(dtypes),
            finite=finite,
            valid=valid,
            scientifically_eligible=False,
            source_protein_id="",
            match_route="ambiguous-explicit-alias",
            sequence_match="unknown" if spec.sequence_dependent else "not-required",
            provenance="unknown",
            factual_status=status,
            requested_action=action,
            reason=_with_action_reason(reason, action),
        )

    source_id = candidate.source_id
    route = candidate.route
    source_file = candidate.source_file
    array = validate_array(source_file, spec.expected_dim, array_cache)
    source = source_benchmark.proteins.get(source_id)
    sequence_match = "not-required"
    if spec.sequence_dependent:
        if source is None:
            sequence_match = "unknown"
        else:
            sequence_match = (
                "true" if source.sequence_sha256 == target.sequence_sha256 else "false"
            )
    eligible = False
    provenance = "%s:%s" % (spec.provenance.compatibility, spec.provenance.label)

    if not array.exists:
        status = "missing"
        factual_reason = "no source array at configured path"
    elif array.error:
        status = "unreadable"
        factual_reason = array.error
    elif tuple_text(spec.expected_dim) != array.observed_shape:
        status = "wrong-dimension"
        factual_reason = "observed %s; expected %s" % (array.observed_shape, tuple_text(spec.expected_dim))
    elif array.finite is False:
        status = "non-finite"
        factual_reason = "array contains NaN or Inf"
    elif array.dtype_supported is False:
        status = "unsupported-dtype"
        factual_reason = "array dtype %s is not a real floating type" % array.dtype
    elif source is None:
        status = "provenance-unknown"
        factual_reason = "source ID is absent from the explicit source benchmark"
    elif spec.sequence_dependent and sequence_match != "true":
        status = "sequence-mismatch"
        factual_reason = "complete source and target sequence SHA-256 values differ"
    elif modality == "text" and not _text_roles_compatible(target, source, spec):
        status = "provenance-incompatible"
        factual_reason = "source and target temporal text roles differ"
    elif not artifact_verification.verified:
        status = "provenance-unknown"
        factual_reason = (
            "published cache, archives, and PFP reference were not "
            "cryptographically verified"
        )
    elif spec.provenance.compatibility == "artifact-scoped":
        if (
            artifact_verification.verified
            and policy == "paper-faithful"
            and target_matches_artifact
            and route == "exact-id"
            and source_id == target.protein_id
        ):
            status = "present-valid"
            eligible = True
            provenance = "verified-artifact:%s" % artifact_verification.artifact_id
            factual_reason = (
                "valid direct-ID array is owned by the cryptographically verified "
                "published benchmark/cache artifact"
            )
        else:
            status = "provenance-unknown"
            factual_reason = (
                "artifact-scoped reuse requires the canonical target fingerprint, "
                "paper-faithful policy, and a direct-ID array from the authenticated cache"
            )
    elif spec.provenance.compatibility == "unknown":
        status = "provenance-unknown"
        factual_reason = "configuration does not establish compatible provenance"
    elif spec.provenance.compatibility == "incompatible":
        status = "provenance-incompatible"
        factual_reason = "configured source and target provenance are incompatible"
    else:
        direct_id_reuse = (
            spec.provenance.allow_direct_id_reuse
            and route == "exact-id"
            and source_id == target.protein_id
        )
        if route.startswith("explicit-alias:") and modality != "prott5":
            status = "provenance-unknown"
            factual_reason = (
                "%s alias reuse is unsupported without an authenticated external "
                "mapping/input artifact" % modality
            )
        elif modality in {"text", "structure", "ppi"} and source_id != target.protein_id:
            status = "provenance-unknown"
            factual_reason = (
                "%s cross-ID reuse is unsupported without an authenticated external "
                "mapping/input artifact" % modality
            )
        elif (
            spec.provenance.requires_mapping_evidence
            and not route.startswith("explicit-alias:")
            and not direct_id_reuse
        ):
            status = "provenance-unknown"
            factual_reason = "required per-protein mapping/source evidence is absent"
        else:
            status = "present-valid"
            eligible = True
            if direct_id_reuse:
                factual_reason = (
                    "valid direct-ID array is eligible under the fixed published "
                    "PPI source/extractor identity"
                )
            else:
                factual_reason = "valid array eligible under declared %s provenance evidence" % route

    action = "reuse" if eligible else "regenerate"
    reason_parts = [factual_reason]
    if candidate.note:
        reason_parts.append(candidate.note)
    if candidate.alias_mapping_evidence:
        reason_parts.append("mapping evidence: %s" % candidate.alias_mapping_evidence)
    reason = _with_action_reason("; ".join(reason_parts), action)
    return InventoryRecord(
        protein_id=target.protein_id,
        sequence_sha256=target.sequence_sha256,
        modality=modality,
        source_directory=spec.directory,
        source_file=str(Path(spec.directory) / source_file.name),
        exists=array.exists,
        observed_shape=array.observed_shape,
        expected_shape=tuple_text(spec.expected_dim),
        dtype=array.dtype,
        finite=array.finite,
        valid=array.valid,
        scientifically_eligible=eligible,
        source_protein_id=source_id,
        match_route=route,
        sequence_match=sequence_match,
        provenance=provenance,
        factual_status=status,
        requested_action=action,
        reason=reason,
    )


def _select_candidate(
    target: ProteinRecord,
    modality: str,
    spec: ModalitySpec,
    source_dir: Path,
    source_benchmark: BenchmarkData,
    source_by_sequence: Dict[str, List[str]],
    alias_entries: List[AliasEntry],
    array_cache: ArrayCache,
    verified_artifact: bool,
) -> Candidate:
    direct_file = source_dir / (target.protein_id + ".npy")
    direct_source = source_benchmark.proteins.get(target.protein_id)
    direct_sequence_ok = (
        direct_source is not None and direct_source.sequence_sha256 == target.sequence_sha256
    )

    if verified_artifact and direct_file.is_file() and direct_source is not None:
        if not spec.sequence_dependent or direct_sequence_ok:
            direct_info = validate_array(direct_file, spec.expected_dim, array_cache)
            if direct_info.valid:
                # Exact artifact ownership is deliberately stronger than an
                # alias. This route prevents aliases from manufacturing scope.
                return Candidate(target.protein_id, "exact-id", direct_file)

    if modality == "prott5" and direct_file.is_file() and direct_sequence_ok:
        direct_info = validate_array(direct_file, spec.expected_dim, array_cache)
        if direct_info.valid:
            return Candidate(target.protein_id, "exact-id", direct_file)

    ambiguous_alias: Optional[Candidate] = None
    if alias_entries:
        mappings = sorted(
            {
                (
                    entry.source_protein_id,
                    entry.mapping_route,
                    entry.source_identity,
                    entry.mapping_evidence,
                )
                for entry in alias_entries
            }
        )
        if len(mappings) != 1:
            detail = ",".join("%s via %s" % (item[0], item[1]) for item in mappings)
            ambiguous_alias = Candidate(
                "", "ambiguous-explicit-alias", source_dir, ambiguity=detail
            )
            if modality != "prott5":
                return ambiguous_alias
        else:
            source_id, mapping_route, source_identity, mapping_evidence = mappings[0]
            alias_candidate = Candidate(
                source_id,
                "explicit-alias:%s" % mapping_route,
                source_dir / (source_id + ".npy"),
                alias_source_identity=source_identity,
                alias_mapping_evidence=mapping_evidence,
            )
            if modality != "prott5" or validate_array(
                alias_candidate.source_file, spec.expected_dim, array_cache
            ).valid:
                return alias_candidate

    if spec.allow_sequence_hash_reuse:
        candidate_ids = [
            source_id
            for source_id in source_by_sequence.get(target.sequence_sha256, [])
            if (source_dir / (source_id + ".npy")).is_file()
        ]
        valid_ids = [
            source_id
            for source_id in candidate_ids
            if validate_array(
                source_dir / (source_id + ".npy"), spec.expected_dim, array_cache
            ).valid
        ]
        if valid_ids:
            selected = valid_ids[0]
            note = ""
            if direct_file.exists():
                direct_info = validate_array(direct_file, spec.expected_dim, array_cache)
                if not direct_sequence_ok:
                    note = "incompatible direct-ID sequence bypassed"
                elif not direct_info.valid:
                    note = "invalid direct-ID array bypassed"
            if len(valid_ids) > 1:
                suffix = "deterministically selected %s from %d valid exact-sequence sources" % (
                    selected,
                    len(valid_ids),
                )
                note = "%s; %s" % (note, suffix) if note else suffix
            return Candidate(
                selected, "sequence-sha256", source_dir / (selected + ".npy"), note=note
            )

    if ambiguous_alias is not None:
        return ambiguous_alias
    if direct_file.exists():
        return Candidate(target.protein_id, "exact-id", direct_file)
    if alias_entries and len(alias_entries) == 1:
        entry = alias_entries[0]
        return Candidate(
            entry.source_protein_id,
            "explicit-alias:%s" % entry.mapping_route,
            source_dir / (entry.source_protein_id + ".npy"),
            alias_source_identity=entry.source_identity,
            alias_mapping_evidence=entry.mapping_evidence,
        )
    return Candidate(target.protein_id, "exact-id-lookup", direct_file)


def _text_roles_compatible(
    target: ProteinRecord, source: ProteinRecord, spec: ModalitySpec
) -> bool:
    if spec.provenance.text_role_policy == "none":
        return True
    if spec.provenance.text_role_policy == "cafa3-mixed-temporal":
        target_role = temporal_text_role(target)
        source_role = temporal_text_role(source)
        allowed_roles = {"current-train-validation", "historical-test"}
        return target_role in allowed_roles and target_role == source_role
    return False


def _with_action_reason(reason: str, action: str) -> str:
    if action == "reuse":
        return "%s; reuse is positively proven" % reason
    return "%s; reuse is not positively proven, so regenerate" % reason


def tuple_text(expected_dim: int) -> str:
    return str((expected_dim,))
