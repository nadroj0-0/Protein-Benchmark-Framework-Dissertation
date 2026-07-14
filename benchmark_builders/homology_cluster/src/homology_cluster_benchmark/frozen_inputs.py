from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Any

from .models import InputSpec, ResolvedInput


SCHEMA_NAME = "homology-cluster-frozen-inputs"
INPUT_NAMES = (
    "uniref90_fasta",
    "idmapping",
    "uniprot_sequences",
    "goa",
    "go_obo",
)
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


def load_frozen_input_manifest(
    path: Path, *, fixture_mode: bool = False
) -> FrozenInputManifest:
    resolved = path.expanduser().resolve()
    raw = resolved.read_bytes()
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Frozen-input manifest is not valid JSON: {resolved}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("Frozen-input manifest root must be an object")
    if payload.get("schema_name") != SCHEMA_NAME or payload.get("schema_version") != 1:
        raise ValueError("Frozen-input manifest has an unsupported schema_name/schema_version")
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
    if not isinstance(raw_entries, list) or len(raw_entries) != len(INPUT_NAMES):
        raise ValueError("Frozen-input manifest must contain exactly five input entries")
    entries: dict[str, dict[str, Any]] = {}
    required = {
        "name", "release", "url", "local_filename", "size_bytes", "sha256",
        "acquisition", "embedded_metadata", "notes",
    }
    for index, entry in enumerate(raw_entries):
        if not isinstance(entry, dict) or not required.issubset(entry):
            raise ValueError(
                f"Frozen-input manifest entry {index} lacks required schema fields"
            )
        name = str(entry["name"])
        if name not in INPUT_NAMES:
            raise ValueError(f"Unknown frozen-input manifest source name: {name!r}")
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
    if set(entries) != set(INPUT_NAMES):
        raise ValueError("Frozen-input manifest is missing one or more required sources")
    if not fixture_mode:
        goa_metadata = entries["goa"]["embedded_metadata"]
        ontology_metadata = entries["go_obo"]["embedded_metadata"]
        if not {"gaf_version", "date_generated", "go_version"}.issubset(goa_metadata):
            raise ValueError("Production GOA manifest entry lacks embedded GAF release metadata")
        if "data_version" not in ontology_metadata:
            raise ValueError("Production ontology manifest entry lacks embedded data_version")

    fingerprint_payload = [entries[name] for name in INPUT_NAMES]
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
    for name in INPUT_NAMES:
        entry = manifest.entries[name]
        spec = specs[name]
        observed = resolved_inputs[name]
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
) -> FrozenInputManifest:
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": 1,
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
            for name in INPUT_NAMES
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return load_frozen_input_manifest(path, fixture_mode=True)
