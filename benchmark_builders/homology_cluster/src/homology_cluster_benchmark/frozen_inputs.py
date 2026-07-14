from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Any

from .models import InputSpec, ResolvedInput


SCHEMA_NAME = "homology-cluster-frozen-inputs"
SHARED_INPUT_NAMES = (
    "uniref90_fasta",
    "idmapping",
    "goa",
    "go_obo",
)
SOURCE_INPUT_NAMES = {
    "sprot-only": ("uniprot_sprot_sequences",),
    "trembl-only": ("uniprot_trembl_sequences",),
    "sprot-and-trembl": ("uniprot_sprot_sequences", "uniprot_trembl_sequences"),
}
SOURCE_POPULATIONS = {
    "uniref90_fasta": "uniref90-clustering-scaffold",
    "idmapping": "uniprotkb-shared-mapping",
    "uniprot_sprot_sequences": "sprot",
    "uniprot_trembl_sequences": "trembl",
    "goa": "uniprotkb-goa",
    "go_obo": "gene-ontology",
}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
PLACEHOLDER_TOKENS = ("replace_me", "placeholder", "example.invalid", "<", ">")


@dataclass(frozen=True)
class FrozenInputManifest:
    path: Path
    sha256: str
    payload: dict[str, Any]
    entries: dict[str, dict[str, Any]]
    source_fingerprint: str
    authoritative_origin_recorded: bool


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _non_placeholder(value: str) -> bool:
    lowered = value.strip().lower()
    return bool(lowered) and not any(token in lowered for token in PLACEHOLDER_TOKENS)


def expected_input_names(uniprot_source_scope: str) -> tuple[str, ...]:
    try:
        selected = SOURCE_INPUT_NAMES[uniprot_source_scope]
    except KeyError as exc:
        raise ValueError(f"Unsupported UniProt source scope: {uniprot_source_scope!r}") from exc
    return (*SHARED_INPUT_NAMES[:2], *selected, *SHARED_INPUT_NAMES[2:])


def load_frozen_input_manifest(
    path: Path, *, uniprot_source_scope: str | None = None, fixture_mode: bool = False
) -> FrozenInputManifest:
    resolved = path.expanduser().resolve()
    raw = resolved.read_bytes()
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Frozen-input manifest is not valid JSON: {resolved}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("Frozen-input manifest root must be an object")
    if payload.get("schema_name") != SCHEMA_NAME or payload.get("schema_version") != 2:
        raise ValueError("Frozen-input manifest has an unsupported schema_name/schema_version")
    manifest_scope = str(payload.get("uniprot_source_scope", ""))
    if manifest_scope not in SOURCE_INPUT_NAMES:
        raise ValueError("Frozen-input manifest has an invalid or missing uniprot_source_scope")
    if uniprot_source_scope is not None and manifest_scope != uniprot_source_scope:
        raise ValueError(
            "Frozen-input manifest source scope mismatch: "
            f"expected {uniprot_source_scope!r}, observed {manifest_scope!r}"
        )
    input_names = expected_input_names(manifest_scope)
    if not _non_placeholder(str(payload.get("created_at", ""))):
        raise ValueError("Frozen-input manifest created_at is missing or a placeholder")

    review = payload.get("review")
    if not isinstance(review, dict):
        raise ValueError("Frozen-input manifest review metadata must be an object")
    synthetic = review.get("status") == "synthetic-fixture"
    if synthetic and not fixture_mode:
        raise ValueError("A synthetic fixture manifest cannot authorize a production build")
    authoritative = (
        review.get("status") == "authoritative-source-reviewed"
        and review.get("authoritative_origin") is True
        and _non_placeholder(str(review.get("reviewed_by", "")))
        and _non_placeholder(str(review.get("evidence", "")))
    )
    if not fixture_mode and not authoritative:
        raise ValueError(
            "Production frozen-input manifest must record an authoritative-source review; "
            "file self-hashes alone are insufficient"
        )
    if fixture_mode and not (synthetic or authoritative):
        raise ValueError("Fixture manifest review status must be synthetic-fixture or reviewed")

    raw_entries = payload.get("inputs")
    if not isinstance(raw_entries, list) or len(raw_entries) != len(input_names):
        raise ValueError(
            f"Frozen-input manifest for {manifest_scope} must contain exactly "
            f"{len(input_names)} input entries"
        )
    entries: dict[str, dict[str, Any]] = {}
    required = {
        "name", "logical_role", "source_population", "release", "url",
        "local_filename", "size_bytes", "sha256", "acquisition",
        "embedded_metadata", "notes",
    }
    for index, entry in enumerate(raw_entries):
        if not isinstance(entry, dict) or not required.issubset(entry):
            raise ValueError(
                f"Frozen-input manifest entry {index} lacks required schema fields"
            )
        name = str(entry["name"])
        if name not in input_names:
            raise ValueError(f"Unknown frozen-input manifest source name: {name!r}")
        if entry["logical_role"] != name:
            raise ValueError(f"Frozen-input manifest logical_role mismatch for {name}")
        if entry["source_population"] != SOURCE_POPULATIONS[name]:
            raise ValueError(f"Frozen-input manifest source_population mismatch for {name}")
        if name in entries:
            raise ValueError(f"Duplicate frozen-input manifest source name: {name}")
        local_filename = str(entry["local_filename"])
        if Path(local_filename).name != local_filename or not _non_placeholder(local_filename):
            raise ValueError(f"Invalid or placeholder local_filename for {name}")
        release = str(entry["release"])
        url = str(entry["url"])
        acquisition = str(entry["acquisition"])
        notes = str(entry["notes"])
        digest = str(entry["sha256"]).lower()
        if not all(_non_placeholder(value) for value in (release, url, acquisition, notes)):
            raise ValueError(f"Frozen-input manifest entry {name} contains a placeholder")
        if not (url.startswith("https://") or url.startswith("ftp://")):
            raise ValueError(f"Frozen-input manifest URL for {name} must be an HTTPS/FTP source")
        if not isinstance(entry["size_bytes"], int) or entry["size_bytes"] <= 0:
            raise ValueError(f"Frozen-input manifest size_bytes for {name} must be positive")
        if not SHA256_RE.fullmatch(digest):
            raise ValueError(f"Frozen-input manifest SHA-256 for {name} is invalid")
        if not isinstance(entry["embedded_metadata"], dict):
            raise ValueError(f"Frozen-input manifest embedded_metadata for {name} must be an object")
        normalized = dict(entry)
        normalized["sha256"] = digest
        entries[name] = normalized
    if set(entries) != set(input_names):
        raise ValueError("Frozen-input manifest is missing one or more required sources")
    uniprot_releases = {
        entries[name]["release"]
        for name in ("uniref90_fasta", "idmapping", *SOURCE_INPUT_NAMES[manifest_scope])
    }
    if len(uniprot_releases) != 1:
        raise ValueError("Frozen-input manifest UniProt/UniRef source releases disagree")
    if not fixture_mode:
        goa_metadata = entries["goa"]["embedded_metadata"]
        ontology_metadata = entries["go_obo"]["embedded_metadata"]
        if not {"gaf_version", "date_generated", "go_version"}.issubset(goa_metadata):
            raise ValueError("Production GOA manifest entry lacks embedded GAF release metadata")
        if "data_version" not in ontology_metadata:
            raise ValueError("Production ontology manifest entry lacks embedded data_version")

    fingerprint_payload = {
        "uniprot_source_scope": manifest_scope,
        "inputs": [entries[name] for name in input_names],
    }
    fingerprint = _sha256_bytes(
        json.dumps(fingerprint_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )
    return FrozenInputManifest(
        path=resolved,
        sha256=_sha256_bytes(raw),
        payload=payload,
        entries=entries,
        source_fingerprint=fingerprint,
        authoritative_origin_recorded=authoritative,
    )


def bind_frozen_inputs(
    manifest: FrozenInputManifest,
    specs: dict[str, InputSpec],
    resolved_inputs: dict[str, ResolvedInput],
) -> dict[str, bool]:
    manifest_scope = str(manifest.payload["uniprot_source_scope"])
    input_names = expected_input_names(manifest_scope)
    if set(specs) != set(input_names) or set(resolved_inputs) != set(input_names):
        raise ValueError("Configured/resolved input roles do not exactly match the selected source scope")
    for name in input_names:
        entry = manifest.entries[name]
        spec = specs[name]
        observed = resolved_inputs[name]
        if spec.name != name or observed.name != name:
            raise ValueError(f"Frozen-input role mismatch for {name}")
        if spec.source_population != entry["source_population"]:
            raise ValueError(f"Configured source_population mismatch for {name}")
        if observed.source_population != entry["source_population"]:
            raise ValueError(f"Resolved source_population mismatch for {name}")
        if spec.release != entry["release"] or observed.release != entry["release"]:
            raise ValueError(f"Frozen-input manifest release mismatch for {name}")
        if spec.path is not None and spec.path.expanduser().resolve().name != entry["local_filename"]:
            raise ValueError(f"Frozen-input manifest local_filename mismatch for {name}")
        if spec.url is not None and spec.url != entry["url"]:
            raise ValueError(f"Frozen-input manifest URL mismatch for {name}")
        if spec.path is None and spec.url is not None:
            source_filename = Path(spec.url.split("?", 1)[0]).name
            if source_filename != entry["local_filename"]:
                raise ValueError(f"Frozen-input manifest source filename mismatch for {name}")
        if spec.expected_sha256 and spec.expected_sha256.lower() != entry["sha256"]:
            raise ValueError(f"Configured SHA-256 disagrees with frozen-input manifest for {name}")
        if observed.size_bytes != entry["size_bytes"]:
            raise ValueError(f"Frozen-input manifest byte-size mismatch for {name}")
        if observed.sha256 != entry["sha256"]:
            raise ValueError(f"Frozen-input manifest SHA-256 mismatch for {name}")
    return {
        "byte_reproducibility": True,
        "recorded_provenance": True,
        "authoritative_origin": manifest.authoritative_origin_recorded,
    }


def verify_frozen_manifest_unchanged(manifest: FrozenInputManifest) -> None:
    if _sha256_bytes(manifest.path.read_bytes()) != manifest.sha256:
        raise ValueError("Frozen-input manifest changed while the run was in progress")


def write_synthetic_fixture_manifest(
    path: Path,
    specs: dict[str, InputSpec],
    resolved_inputs: dict[str, ResolvedInput],
    uniprot_source_scope: str,
) -> FrozenInputManifest:
    input_names = expected_input_names(uniprot_source_scope)
    if set(specs) != set(input_names) or set(resolved_inputs) != set(input_names):
        raise ValueError("Synthetic fixture inputs do not match the selected source scope")
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": 2,
        "uniprot_source_scope": uniprot_source_scope,
        "created_at": "synthetic-fixture-generated",
        "review": {
            "status": "synthetic-fixture",
            "authoritative_origin": False,
            "reviewed_by": "not-applicable-fixture",
            "evidence": "synthetic fixture bytes generated or bundled for software tests",
        },
        "inputs": [
            {
                "name": name,
                "logical_role": name,
                "source_population": SOURCE_POPULATIONS[name],
                "release": specs[name].release,
                "url": specs[name].url
                or f"https://synthetic.invalid/{resolved_inputs[name].resolved_path.name}",
                "local_filename": resolved_inputs[name].resolved_path.name,
                "size_bytes": resolved_inputs[name].size_bytes,
                "sha256": resolved_inputs[name].sha256,
                "acquisition": "synthetic-fixture",
                "embedded_metadata": {},
                "notes": "Synthetic fixture only; never dissertation-production eligible.",
            }
            for name in input_names
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return load_frozen_input_manifest(
        path, uniprot_source_scope=uniprot_source_scope, fixture_mode=True
    )
