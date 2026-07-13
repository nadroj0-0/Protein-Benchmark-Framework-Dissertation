from dataclasses import dataclass
from pathlib import Path
import re
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .benchmark import sequence_index, temporal_text_role
from .models import (
    AliasEntry,
    ArrayInfo,
    BenchmarkData,
    InventoryRecord,
    InventoryResult,
    MODALITIES,
    ModalitySpec,
    PlannerConfig,
    ProteinRecord,
)


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
) -> InventoryResult:
    if policy not in {"paper-faithful", "maximize-coverage"}:
        raise InventoryError("Unsupported action policy: %s" % policy)
    aliases = aliases or {}
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
        source_dir = embedding_cache / config.modalities[modality].directory
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
            source_dir = embedding_cache / spec.directory
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
        action = "manual-review"
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
            reason=_with_policy_reason(reason, action, policy),
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
    elif spec.provenance.compatibility == "unknown":
        status = "provenance-unknown"
        factual_reason = "configuration does not establish compatible provenance"
    elif spec.provenance.compatibility == "incompatible":
        status = "provenance-incompatible"
        factual_reason = "configured source and target provenance are incompatible"
    else:
        if spec.provenance.requires_mapping_evidence and not route.startswith("explicit-alias:"):
            status = "provenance-unknown"
            factual_reason = "required per-protein mapping/source evidence is absent"
        elif route.startswith("explicit-alias:") and (
            candidate.alias_source_identity != spec.provenance.source_identity
        ):
            status = "provenance-incompatible"
            factual_reason = "alias source identity does not match configured modality source identity"
        elif route.startswith("explicit-alias:") and not _mapping_evidence_compatible(
            candidate.alias_mapping_evidence, modality, target, spec
        )[0]:
            status = "provenance-unknown"
            factual_reason = _mapping_evidence_compatible(
                candidate.alias_mapping_evidence, modality, target, spec
            )[1]
        elif modality == "structure" and source_id != target.protein_id and not route.startswith("explicit-alias:"):
            status = "provenance-unknown"
            factual_reason = "structure cross-ID reuse requires an explicit compatible alias"
        elif modality in {"text", "structure", "ppi"} and (
            source_id != target.protein_id and not route.startswith("explicit-alias:")
        ):
            status = "provenance-unknown"
            factual_reason = "%s cross-ID reuse requires an explicit alias" % modality
        else:
            status = "present-valid"
            eligible = True
            factual_reason = "valid array eligible under declared %s provenance evidence" % route

    action = _requested_action(status, policy, spec)
    reason_parts = [factual_reason]
    if candidate.note:
        reason_parts.append(candidate.note)
    if candidate.alias_mapping_evidence:
        reason_parts.append("mapping evidence: %s" % candidate.alias_mapping_evidence)
    reason = _with_policy_reason("; ".join(reason_parts), action, policy)
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
) -> Candidate:
    direct_file = source_dir / (target.protein_id + ".npy")
    direct_source = source_benchmark.proteins.get(target.protein_id)
    direct_sequence_ok = (
        direct_source is not None and direct_source.sequence_sha256 == target.sequence_sha256
    )

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


def _mapping_evidence_compatible(
    evidence_text: str,
    modality: str,
    target: ProteinRecord,
    spec: ModalitySpec,
) -> Tuple[bool, str]:
    evidence: Dict[str, str] = {}
    for item in evidence_text.split(";"):
        key, separator, value = item.partition(":")
        if not separator or not key.strip() or not value.strip():
            return False, "mapping evidence must use semicolon-separated key:value fields"
        evidence[key.strip().lower()] = value.strip()

    required = {
        "prott5": {"sequence-sha256"},
        "text": {"description-sha256", "temporal-context"},
        "structure": {"structure-source", "structure-version"},
        "ppi": {"string-id", "string-release"},
    }[modality]
    missing = required - set(evidence)
    if missing:
        return False, "mapping evidence is missing fields: %s" % ",".join(sorted(missing))

    if modality == "prott5":
        if evidence["sequence-sha256"].lower() != target.sequence_sha256:
            return False, "alias sequence-sha256 evidence does not match the target sequence"
    elif modality == "text":
        if re.fullmatch(r"[0-9a-fA-F]{64}", evidence["description-sha256"]) is None:
            return False, "text description-sha256 evidence must be a complete SHA-256"
        identity_tokens = _identity_tokens(spec.provenance.source_identity)
        if evidence["temporal-context"].lower() not in identity_tokens:
            return False, "text temporal-context evidence is absent from the configured source identity"
    elif modality == "structure":
        identity_tokens = _identity_tokens(spec.provenance.source_identity)
        if any(
            evidence[key].lower() not in identity_tokens
            for key in ("structure-source", "structure-version")
        ):
            return False, "structure source/version evidence does not match the configured source identity"
    elif modality == "ppi":
        if evidence["string-release"].lower() not in _identity_tokens(
            spec.provenance.source_identity
        ):
            return False, "STRING release evidence does not match the configured source identity"
    return True, ""


def _identity_tokens(identity: str) -> Set[str]:
    return {token.strip().lower() for token in identity.split("|") if token.strip()}


def _requested_action(status: str, policy: str, spec: ModalitySpec) -> str:
    if status == "present-valid":
        return "reuse"
    if policy == "paper-faithful":
        if status in {"missing", "unreadable"}:
            return "leave-masked"
        return "manual-review"
    if spec.provenance.compatibility in {"unknown", "incompatible"}:
        return "manual-review"
    if status in {"provenance-unknown", "provenance-incompatible"}:
        return "manual-review"
    if status == "missing":
        return spec.missing_action
    return spec.invalid_action


def _with_policy_reason(reason: str, action: str, policy: str) -> str:
    if action == "leave-masked":
        return "%s; %s preserves PFP zero-vector/mask behavior" % (reason, policy)
    if action == "reuse":
        return "%s; reuse is eligible under the declared evidence" % reason
    return "%s; %s policy requests %s" % (reason, policy, action)


def tuple_text(expected_dim: int) -> str:
    return str((expected_dim,))
