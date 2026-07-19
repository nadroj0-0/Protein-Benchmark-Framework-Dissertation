from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
import json
import logging
from pathlib import Path
import pickle
import shutil
import time
import uuid

from .config import SUPERVISOR_EVIDENCE_CODES
from .frozen_inputs import (
    FrozenInputManifest,
    SOURCE_INPUT_NAMES,
    expected_input_names,
    load_frozen_input_manifest,
)
from .goa import load_goa
from .idmapping import load_uniref90_mappings
from .inputs import sha256_file
from .mapping import canonicalize_goa_accessions, load_requested_proteins_from_sources
from .models import GoaLoadResult, MappingDecision, ProteinCatalog
from .ontology import Ontology
from .uniref import UniRefIndex


LOGGER = logging.getLogger(__name__)

CACHE_SCHEMA_NAME = "homology-cluster-common-preprocessing"
CACHE_SCHEMA_VERSION = 3
CACHE_MARKER = "CACHE_COMPLETE.json"
STATE_FILE = "preprocessing_state.pkl.gz"
STATE_FORMAT = "homology-common-preprocessing-state-v1"
UNIREF_INDEX_FILE = "uniref90.sqlite"
FROZEN_MANIFEST_FILE = "frozen_input_manifest.json"
GOA_RECORD_FILE = "goa/qualifying_annotations.raw.jsonl.gz"
GOA_EXCLUDED_FILE = "goa/excluded_annotations.sample.jsonl.gz"
COLLISION_REPORT_FILE = "uniprot_accession_collisions.tsv.gz"

# A cache is invalidated if any source file that defines the shared scientific
# preprocessing changes. Threshold-specific MMseqs/splitting code is deliberately absent.
PREPROCESSING_SOURCE_FILES = (
    "common_cache.py",
    "config.py",
    "goa.py",
    "idmapping.py",
    "inputs.py",
    "mapping.py",
    "models.py",
    "ontology.py",
    "uniref.py",
)

# Schema 2 was produced by the first SAN cache build. Its CLI executed this
# module with ``python -m``, which made the wrapper dataclass pickle as
# ``__main__.CommonPreprocessingState``. The complete producer fingerprint is
# pinned here so compatibility cannot silently extend to an unknown cache.
SCHEMA_V2_PREPROCESSING_SOURCE_SHA256 = {
    "common_cache.py": "07eb91fe7cfa8fd3bb8c23f62d633c56ebf5e4ce1905c755dd2e6006cf146994",
    "config.py": "f41fdfbe5c5a2c0f9288899f46d1e527f07428df00976b547607543eab29b0ed",
    "goa.py": "db0e4a7d0bc8124c1ea00c88c7a7694ac9ae52a76a6a031254fa562fb51bc055",
    "idmapping.py": "c23b8232292ff3e868364014692c2df3f73234444ed198ed50b10509210bbefd",
    "inputs.py": "7da70361cfc0db42bb529beabda64d01c09898b1e1d2476380de249e1e9983ad",
    "mapping.py": "4f91448aab74271cc4e1562de46ed7ed685e8b5b39bbf7403906066d2333f9ff",
    "models.py": "f8529658f3689a1499ad00b49d179c1f50ef4328c8824cbe54a7e7fcb19cee53",
    "ontology.py": "138374a3338c9f9e38411298284810d16a039b4fc0abaae04ab757ed0ae2439b",
    "uniref.py": "9a03bb24f3e6697730550e8420ab70adb9eb6dd344232c6a2c63a16ac98164c3",
}


@dataclass(frozen=True)
class CommonPreprocessingState:
    goa: GoaLoadResult
    catalog: ProteinCatalog
    decisions: list[MappingDecision]
    requested_raw: set[str]


@dataclass(frozen=True)
class LoadedCommonPreprocessing:
    root: Path
    marker_sha256: str
    payload: dict[str, object]
    uniref: UniRefIndex
    goa: GoaLoadResult
    catalog: ProteinCatalog
    decisions: list[MappingDecision]
    requested_raw: set[str]


class _SchemaV2StateUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str) -> object:
        if module == "__main__" and name == "CommonPreprocessingState":
            return CommonPreprocessingState
        return super().find_class(module, name)


def _load_common_preprocessing_state(
    path: Path, schema_version: int
) -> CommonPreprocessingState:
    with gzip.open(path, "rb") as handle:
        if schema_version == 2:
            state = _SchemaV2StateUnpickler(handle).load()
        else:
            payload = pickle.load(handle)
            if not isinstance(payload, dict) or payload.get("state_format") != STATE_FORMAT:
                raise ValueError("Common preprocessing cache state has an unsupported format")
            expected = {"state_format", "goa", "catalog", "decisions", "requested_raw"}
            if set(payload) != expected:
                raise ValueError("Common preprocessing cache state fields are malformed")
            state = CommonPreprocessingState(
                goa=payload["goa"],
                catalog=payload["catalog"],
                decisions=payload["decisions"],
                requested_raw=payload["requested_raw"],
            )
    if not isinstance(state, CommonPreprocessingState):
        raise ValueError("Common preprocessing cache state has an unexpected type")
    return state


def _json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _source_hashes() -> dict[str, str]:
    root = Path(__file__).resolve().parent
    return {name: sha256_file(root / name) for name in PREPROCESSING_SOURCE_FILES}


def _cache_root(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved.is_file() and resolved.name == CACHE_MARKER:
        resolved = resolved.parent
    return resolved


def common_cache_root(path: Path) -> Path:
    """Return the cache directory for either a directory or marker-file argument."""
    return _cache_root(path)


def _directory_files(root: Path) -> list[Path]:
    return sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.name != CACHE_MARKER
    )


def _cache_policy(
    *,
    source_scope: str,
    include_relationships: bool,
    strict_qc: bool,
    excluded_sample_per_reason: int,
) -> dict[str, object]:
    return {
        "uniprot_source_scope": source_scope,
        "evidence_codes": sorted(SUPERVISOR_EVIDENCE_CODES),
        "include_relationships": include_relationships,
        "strict_qc": strict_qc,
        "excluded_sample_per_reason": excluded_sample_per_reason,
        "uniprot_release": "2026_02",
        "goa_release": "234",
        "ontology_release": "releases/2026-06-15",
    }


def _validate_frozen_metadata(
    manifest: FrozenInputManifest,
    ontology: Ontology,
    headers: dict[str, str],
) -> None:
    if headers.get("gaf_version") != "2.2":
        raise ValueError("Common preprocessing requires GOA GAF 2.2")
    if "2026-06-17" not in headers.get("date_generated", ""):
        raise ValueError("Common preprocessing requires GOA 234 from 2026-06-17")
    if "2026-06-15" not in headers.get("go_version", ""):
        raise ValueError("Common preprocessing requires GOA's 2026-06-15 GO version")
    if ontology.data_version != "releases/2026-06-15":
        raise ValueError(
            "Common preprocessing requires ontology data-version releases/2026-06-15"
        )
    declared_goa = manifest.entries["goa"]["embedded_metadata"]
    for key, expected in declared_goa.items():
        if expected and str(expected) not in str(headers.get(key, "")):
            raise ValueError(
                f"Frozen manifest GOA metadata mismatch for {key}: "
                f"expected {expected!r}, observed {headers.get(key)!r}"
            )
    expected_ontology = manifest.entries["go_obo"]["embedded_metadata"].get(
        "data_version"
    )
    if expected_ontology and ontology.data_version != expected_ontology:
        raise ValueError(
            "Frozen manifest ontology metadata mismatch: "
            f"expected {expected_ontology!r}, observed {ontology.data_version!r}"
        )


def _validate_input_paths(
    manifest: FrozenInputManifest,
    paths: dict[str, Path],
) -> dict[str, dict[str, object]]:
    expected_names = set(expected_input_names(str(manifest.payload["uniprot_source_scope"])))
    if set(paths) != expected_names:
        raise ValueError(
            "Common-cache inputs do not exactly match the frozen source scope: "
            f"expected={sorted(expected_names)}, observed={sorted(paths)}"
        )
    bindings: dict[str, dict[str, object]] = {}
    for name in sorted(expected_names):
        path = paths[name].expanduser().resolve()
        entry = manifest.entries[name]
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"Common-cache input is missing or empty: {name}: {path}")
        if path.name != entry["local_filename"]:
            raise ValueError(
                f"Common-cache input filename mismatch for {name}: "
                f"expected={entry['local_filename']!r}, observed={path.name!r}"
            )
        observed_size = path.stat().st_size
        if observed_size != entry["size_bytes"]:
            raise ValueError(
                f"Common-cache input size mismatch for {name}: "
                f"expected={entry['size_bytes']}, observed={observed_size}"
            )
        observed_sha = sha256_file(path)
        if observed_sha != entry["sha256"]:
            raise ValueError(
                f"Common-cache input SHA-256 mismatch for {name}: "
                f"expected={entry['sha256']}, observed={observed_sha}"
            )
        bindings[name] = {
            "release": entry["release"],
            "size_bytes": observed_size,
            "sha256": observed_sha,
            "source_population": entry["source_population"],
            "url": entry["url"],
            "local_filename": entry["local_filename"],
        }
    return bindings


def inspect_common_preprocessing_cache(
    path: Path,
    *,
    expected_source_scope: str | None = None,
    expected_input_sha256: dict[str, str] | None = None,
    expected_policy: dict[str, object] | None = None,
    verify_file_hashes: bool = False,
) -> dict[str, object]:
    root = _cache_root(path)
    marker = root / CACHE_MARKER
    if not root.is_dir() or not marker.is_file():
        raise ValueError(f"Common preprocessing cache is incomplete: {root}")
    payload = json.loads(marker.read_text(encoding="utf-8"))
    schema_version = payload.get("schema_version")
    if (
        payload.get("schema_name") != CACHE_SCHEMA_NAME
        or schema_version not in {2, CACHE_SCHEMA_VERSION}
        or payload.get("complete") is not True
    ):
        raise ValueError("Common preprocessing cache marker has an unsupported contract")
    source_scope = str(payload.get("uniprot_source_scope", ""))
    if source_scope not in SOURCE_INPUT_NAMES:
        raise ValueError("Common preprocessing cache has an invalid source scope")
    if expected_source_scope is not None and source_scope != expected_source_scope:
        raise ValueError(
            "Common preprocessing cache source-scope mismatch: "
            f"expected={expected_source_scope}, observed={source_scope}"
        )
    observed_source_hashes = payload.get("preprocessing_source_sha256")
    expected_source_hashes = (
        SCHEMA_V2_PREPROCESSING_SOURCE_SHA256
        if schema_version == 2 else _source_hashes()
    )
    if observed_source_hashes != expected_source_hashes:
        raise ValueError(
            "Common preprocessing cache was produced by unsupported preprocessing code"
        )
    if expected_policy is not None and payload.get("policy") != expected_policy:
        raise ValueError("Common preprocessing cache policy does not match this run")
    bindings = payload.get("input_bindings")
    if not isinstance(bindings, dict):
        raise ValueError("Common preprocessing cache lacks input bindings")
    if expected_input_sha256 is not None:
        if set(bindings) != set(expected_input_sha256):
            raise ValueError("Common preprocessing cache input-role set does not match")
        for name, expected_sha in expected_input_sha256.items():
            observed = bindings.get(name)
            if not isinstance(observed, dict) or observed.get("sha256") != expected_sha:
                raise ValueError(
                    f"Common preprocessing cache input SHA-256 mismatch for {name}"
                )
    files = payload.get("files")
    if not isinstance(files, list) or not files:
        raise ValueError("Common preprocessing cache has no file manifest")
    expected_paths: set[str] = set()
    for entry in files:
        if not isinstance(entry, dict):
            raise ValueError("Common preprocessing cache file entry is malformed")
        relative = str(entry.get("path", ""))
        candidate = (root / relative).resolve()
        if not relative or root not in candidate.parents or not candidate.is_file():
            raise ValueError(f"Common preprocessing cache file is missing: {relative}")
        expected_paths.add(relative)
        if candidate.stat().st_size != entry.get("size_bytes"):
            raise ValueError(f"Common preprocessing cache file-size mismatch: {relative}")
        if verify_file_hashes and sha256_file(candidate) != entry.get("sha256"):
            raise ValueError(f"Common preprocessing cache file hash mismatch: {relative}")
    observed_paths = {
        path.relative_to(root).as_posix() for path in _directory_files(root)
    }
    if observed_paths != expected_paths:
        raise ValueError(
            "Common preprocessing cache files do not reconcile with its marker: "
            f"missing={sorted(expected_paths - observed_paths)}, "
            f"extra={sorted(observed_paths - expected_paths)}"
        )
    return payload


def build_common_preprocessing_cache(
    output_dir: Path,
    work_dir: Path,
    frozen_input_manifest: Path,
    input_paths: dict[str, Path],
    *,
    source_scope: str,
    include_relationships: bool = True,
    strict_qc: bool = True,
    excluded_sample_per_reason: int = 1000,
    fixture_mode: bool = False,
    replace_existing: bool = False,
) -> Path:
    output = output_dir.expanduser().resolve()
    work = work_dir.expanduser().resolve()
    manifest = load_frozen_input_manifest(
        frozen_input_manifest,
        uniprot_source_scope=source_scope,
        fixture_mode=fixture_mode,
    )
    bindings = _validate_input_paths(manifest, input_paths)
    policy = _cache_policy(
        source_scope=source_scope,
        include_relationships=include_relationships,
        strict_qc=strict_qc,
        excluded_sample_per_reason=excluded_sample_per_reason,
    )
    expected_hashes = {
        name: str(entry["sha256"]) for name, entry in bindings.items()
    }
    if output.exists():
        try:
            inspect_common_preprocessing_cache(
                output,
                expected_source_scope=source_scope,
                expected_input_sha256=expected_hashes,
                expected_policy=policy,
                verify_file_hashes=True,
            )
            LOGGER.info("Common preprocessing cache is already complete: %s", output)
            return output
        except (OSError, ValueError):
            if not replace_existing:
                raise

    work.mkdir(parents=True, exist_ok=True)
    run_work = work / f"common-cache-build-{uuid.uuid4().hex}"
    run_work.mkdir()
    stage = output.parent / f".{output.name}.staging-{uuid.uuid4().hex}"
    old = output.parent / f".{output.name}.obsolete-{uuid.uuid4().hex}"
    started = time.monotonic()
    try:
        ontology = Ontology(input_paths["go_obo"], include_relationships)
        goa = load_goa(
            input_paths["goa"],
            ontology,
            SUPERVISOR_EVIDENCE_CODES,
            strict_malformed=strict_qc,
            spool_dir=run_work / "goa",
            excluded_sample_per_reason=excluded_sample_per_reason,
        )
        _validate_frozen_metadata(manifest, ontology, goa.headers)
        uniref = UniRefIndex.build(
            input_paths["uniref90_fasta"], run_work / UNIREF_INDEX_FILE
        )
        requested_raw = set(goa.qualifying_accessions or goa.annotations)
        source_names = {
            "uniprot_sprot_sequences": "sprot",
            "uniprot_trembl_sequences": "trembl",
        }
        selected_sources = {
            source_names[name]: input_paths[name]
            for name in SOURCE_INPUT_NAMES[source_scope]
        }
        catalog = load_requested_proteins_from_sources(
            selected_sources,
            requested_raw,
            strict_collisions=not fixture_mode,
            collision_database=run_work / "uniprot_accessions.sqlite",
            collision_report=run_work / COLLISION_REPORT_FILE,
        )
        goa = canonicalize_goa_accessions(goa, catalog)
        decisions = load_uniref90_mappings(
            input_paths["idmapping"], requested_raw, catalog, uniref
        )

        stage.mkdir(parents=True)
        (stage / "goa").mkdir()
        shutil.copy2(run_work / UNIREF_INDEX_FILE, stage / UNIREF_INDEX_FILE)
        shutil.copy2(
            run_work / COLLISION_REPORT_FILE,
            stage / COLLISION_REPORT_FILE,
        )
        if goa.record_spool is not None:
            shutil.copy2(goa.record_spool, stage / GOA_RECORD_FILE)
            goa.record_spool = Path(GOA_RECORD_FILE)
        if goa.excluded_spool is not None:
            shutil.copy2(goa.excluded_spool, stage / GOA_EXCLUDED_FILE)
            goa.excluded_spool = Path(GOA_EXCLUDED_FILE)
        state = {
            "state_format": STATE_FORMAT,
            "goa": goa,
            "catalog": catalog,
            "decisions": decisions,
            "requested_raw": requested_raw,
        }
        with gzip.open(stage / STATE_FILE, "wb", compresslevel=1) as handle:
            pickle.dump(state, handle, protocol=pickle.HIGHEST_PROTOCOL)
        shutil.copy2(manifest.path, stage / FROZEN_MANIFEST_FILE)

        file_manifest = [
            {
                "path": path.relative_to(stage).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in _directory_files(stage)
        ]
        payload: dict[str, object] = {
            "schema_name": CACHE_SCHEMA_NAME,
            "schema_version": CACHE_SCHEMA_VERSION,
            "complete": True,
            "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "uniprot_source_scope": source_scope,
            "policy": policy,
            "frozen_input_manifest_sha256": manifest.sha256,
            "frozen_input_source_fingerprint": manifest.source_fingerprint,
            "input_bindings": bindings,
            "preprocessing_source_sha256": _source_hashes(),
            "counts": {
                "uniref90_entries": uniref.count(),
                "qualifying_raw_accessions": len(requested_raw),
                "canonical_qualifying_accessions": len(goa.qualifying_accessions),
                "loaded_canonical_sequences": len(catalog.records),
                "mapping_decisions": len(decisions),
            },
            "files": file_manifest,
            "note": (
                "This cache ends immediately before threshold-specific MMseqs2 execution. "
                "It does not contain clusters, splits, term universes, or PFP outputs."
            ),
        }
        _json(stage / CACHE_MARKER, payload)
        inspect_common_preprocessing_cache(
            stage,
            expected_source_scope=source_scope,
            expected_input_sha256=expected_hashes,
            expected_policy=policy,
            verify_file_hashes=True,
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists():
            output.rename(old)
        stage.rename(output)
        shutil.rmtree(old, ignore_errors=True)
        return output
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        if old.exists() and not output.exists():
            old.rename(output)
        raise
    finally:
        shutil.rmtree(run_work, ignore_errors=True)


def load_common_preprocessing_cache(
    path: Path,
    *,
    source_scope: str,
    expected_input_sha256: dict[str, str],
    include_relationships: bool,
    strict_qc: bool,
    excluded_sample_per_reason: int,
    frozen_input_manifest_sha256: str,
) -> LoadedCommonPreprocessing:
    root = _cache_root(path)
    policy = _cache_policy(
        source_scope=source_scope,
        include_relationships=include_relationships,
        strict_qc=strict_qc,
        excluded_sample_per_reason=excluded_sample_per_reason,
    )
    payload = inspect_common_preprocessing_cache(
        root,
        expected_source_scope=source_scope,
        expected_input_sha256=expected_input_sha256,
        expected_policy=policy,
        verify_file_hashes=True,
    )
    if payload.get("frozen_input_manifest_sha256") != frozen_input_manifest_sha256:
        raise ValueError(
            "Common preprocessing cache was not built from this frozen-input manifest"
        )
    state = _load_common_preprocessing_state(
        root / STATE_FILE, int(payload["schema_version"])
    )
    goa = state.goa
    if goa.record_spool is not None:
        goa.record_spool = root / goa.record_spool
    if goa.excluded_spool is not None:
        goa.excluded_spool = root / goa.excluded_spool
    uniref = UniRefIndex(root / UNIREF_INDEX_FILE)
    expected_count = int(payload["counts"]["uniref90_entries"])  # type: ignore[index]
    if uniref.count() != expected_count:
        raise ValueError("Common preprocessing cache UniRef index count changed")
    return LoadedCommonPreprocessing(
        root=root,
        marker_sha256=sha256_file(root / CACHE_MARKER),
        payload=payload,
        uniref=uniref,
        goa=goa,
        catalog=state.catalog,
        decisions=state.decisions,
        requested_raw=state.requested_raw,
    )


def _parse_expected_sha(values: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--expected-input-sha256 must be ROLE=SHA256")
        name, digest = value.split("=", 1)
        if not name or len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
            raise ValueError(f"Invalid expected input SHA-256 binding: {value!r}")
        result[name] = digest
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build or verify shared homology preprocessing")
    subparsers = parser.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("build")
    build.add_argument("--output-dir", type=Path, required=True)
    build.add_argument("--work-dir", type=Path, required=True)
    build.add_argument("--frozen-input-manifest", type=Path, required=True)
    build.add_argument("--source-scope", choices=tuple(SOURCE_INPUT_NAMES), required=True)
    build.add_argument("--uniref90-fasta", type=Path, required=True)
    build.add_argument("--idmapping", type=Path, required=True)
    build.add_argument("--uniprot-sprot-sequences", type=Path)
    build.add_argument("--uniprot-trembl-sequences", type=Path)
    build.add_argument("--goa", type=Path, required=True)
    build.add_argument("--go-obo", type=Path, required=True)
    build.add_argument("--excluded-sample-per-reason", type=int, default=1000)
    build.add_argument("--fixture-mode", action="store_true")
    build.add_argument("--replace-existing", action="store_true")
    verify = subparsers.add_parser("verify")
    verify.add_argument("--cache-dir", type=Path, required=True)
    verify.add_argument("--source-scope", choices=tuple(SOURCE_INPUT_NAMES))
    verify.add_argument("--expected-input-sha256", action="append", default=[])
    verify.add_argument("--full-hashes", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.command == "verify":
        payload = inspect_common_preprocessing_cache(
            args.cache_dir,
            expected_source_scope=args.source_scope,
            expected_input_sha256=(
                _parse_expected_sha(args.expected_input_sha256)
                if args.expected_input_sha256 else None
            ),
            expected_policy=(
                _cache_policy(
                    source_scope=args.source_scope,
                    include_relationships=True,
                    strict_qc=True,
                    excluded_sample_per_reason=1000,
                )
                if args.source_scope else None
            ),
            verify_file_hashes=args.full_hashes,
        )
        print(json.dumps({
            "status": "valid",
            "cache_dir": str(_cache_root(args.cache_dir)),
            "counts": payload["counts"],
        }, sort_keys=True))
        return 0
    source_inputs = {
        "uniref90_fasta": args.uniref90_fasta,
        "idmapping": args.idmapping,
        "goa": args.goa,
        "go_obo": args.go_obo,
    }
    if args.source_scope != "trembl-only":
        if args.uniprot_sprot_sequences is None:
            raise ValueError("Selected source scope requires --uniprot-sprot-sequences")
        source_inputs["uniprot_sprot_sequences"] = args.uniprot_sprot_sequences
    if args.source_scope != "sprot-only":
        if args.uniprot_trembl_sequences is None:
            raise ValueError("Selected source scope requires --uniprot-trembl-sequences")
        source_inputs["uniprot_trembl_sequences"] = args.uniprot_trembl_sequences
    output = build_common_preprocessing_cache(
        args.output_dir,
        args.work_dir,
        args.frozen_input_manifest,
        source_inputs,
        source_scope=args.source_scope,
        excluded_sample_per_reason=args.excluded_sample_per_reason,
        fixture_mode=args.fixture_mode,
        replace_existing=args.replace_existing,
    )
    print(json.dumps({"status": "complete", "cache_dir": str(output)}, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
