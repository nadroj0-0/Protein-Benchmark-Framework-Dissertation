from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
import logging
import os
from pathlib import Path
import shutil
import time
import uuid

from .config import BuildConfig
from .inputs import open_text, sha256_file
from .mmseqs import ClusterIndex, MMseqsRuntime
from .uniref import UniRefIndex


LOGGER = logging.getLogger(__name__)

CACHE_SCHEMA_NAME = "homology-mmseqs-cluster-assignments"
CACHE_SCHEMA_VERSION = 1
CACHE_MARKER = "CACHE_COMPLETE.json"
CACHE_ROOT_MARKER = "CLUSTER_CACHE_ROOT.json"
CACHE_ROOT_SCHEMA_NAME = "homology-mmseqs-cluster-cache-root"
ASSIGNMENTS_FILE = "cluster_assignments.tsv.gz"
COMMANDS_FILE = "mmseqs_commands.tsv"


@dataclass(frozen=True)
class LoadedClusterCache:
    root: Path
    marker_sha256: str
    payload: dict[str, object]
    assignments: Path


def _json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _canonical_hash(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _identity_label(identity: float) -> str:
    percent = identity * 100
    if abs(percent - round(percent)) < 1e-10:
        return f"identity_{int(round(percent)):02d}"
    return "identity_" + f"{percent:.10g}".replace(".", "p")


def cluster_cache_root(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved.is_file() and resolved.name == CACHE_ROOT_MARKER:
        return resolved.parent
    return resolved


def _root_marker_payload() -> dict[str, object]:
    return {
        "schema_name": CACHE_ROOT_SCHEMA_NAME,
        "schema_version": CACHE_SCHEMA_VERSION,
        "cache_schema_name": CACHE_SCHEMA_NAME,
        "role": "persistent validated MMseqs2 cluster-assignment cache root",
        "note": (
            "Each child cache is independently keyed by UniRef90 and the exact MMseqs2 "
            "scientific/runtime contract. Annotation and split policy are downstream."
        ),
    }


def initialize_cluster_cache_root(path: Path) -> Path:
    root = cluster_cache_root(path)
    root.mkdir(parents=True, exist_ok=True)
    marker = root / CACHE_ROOT_MARKER
    expected = _root_marker_payload()
    if marker.exists():
        observed = json.loads(marker.read_text(encoding="utf-8"))
        if observed != expected:
            raise ValueError(f"Cluster-cache root marker is incompatible: {marker}")
        return marker
    temporary = root / f".{CACHE_ROOT_MARKER}.partial-{uuid.uuid4().hex}"
    _json(temporary, expected)
    try:
        os.link(temporary, marker)
        temporary.unlink()
    except OSError:
        if not marker.exists():
            raise
        temporary.unlink(missing_ok=True)
    observed = json.loads(marker.read_text(encoding="utf-8"))
    if observed != expected:
        raise ValueError(f"Cluster-cache root marker is incompatible: {marker}")
    return marker


def inspect_cluster_cache_root(path: Path) -> dict[str, object]:
    root = cluster_cache_root(path)
    marker = root / CACHE_ROOT_MARKER
    if not marker.is_file():
        raise ValueError(f"Cluster-cache root marker is missing: {marker}")
    payload = json.loads(marker.read_text(encoding="utf-8"))
    if payload != _root_marker_payload():
        raise ValueError(f"Cluster-cache root marker is incompatible: {marker}")
    return payload


def cluster_cache_contract_from_values(
    *,
    uniref90_release: str,
    uniref90_sha256: str,
    identity: float,
    coverage: float,
    cov_mode: int,
    cluster_mode: int,
    alignment_mode: int,
    cluster_reassign: int,
    sensitivity: float,
    evalue: float,
    expected_mmseqs_version: str,
    observed_mmseqs_version: str,
    mmseqs_executable_sha256: str,
) -> dict[str, object]:
    if not uniref90_release:
        raise ValueError("Cluster cache requires a UniRef90 release")
    if len(uniref90_sha256) != 64 or any(
        char not in "0123456789abcdef" for char in uniref90_sha256
    ):
        raise ValueError("Cluster cache requires a lowercase UniRef90 SHA-256")
    if not expected_mmseqs_version or not observed_mmseqs_version:
        raise ValueError("Cluster cache requires exact expected and observed MMseqs2 versions")
    if len(mmseqs_executable_sha256) != 64 or any(
        char not in "0123456789abcdef" for char in mmseqs_executable_sha256
    ):
        raise ValueError("Cluster cache requires a lowercase MMseqs2 executable SHA-256")
    return {
        "schema_name": CACHE_SCHEMA_NAME,
        "schema_version": CACHE_SCHEMA_VERSION,
        "input": {
            "uniref90_release": uniref90_release,
            "uniref90_sha256": uniref90_sha256,
        },
        "method": {
            "workflow": "cluster",
            "identity_fraction": f"{identity:.12g}",
            "coverage": f"{coverage:.12g}",
            "cov_mode": cov_mode,
            "cluster_mode": cluster_mode,
            "alignment_mode": alignment_mode,
            "seq_id_mode": 0,
            "cluster_reassign": cluster_reassign,
            "sensitivity": f"{sensitivity:.12g}",
            "evalue": f"{evalue:.12g}",
            "createdb_dbtype": 1,
            "createdb_shuffle": 0,
        },
        "runtime": {
            "expected_mmseqs_version": expected_mmseqs_version,
            "observed_mmseqs_version": observed_mmseqs_version,
            "mmseqs_executable_sha256": mmseqs_executable_sha256,
        },
    }


def cluster_cache_contract(
    config: BuildConfig,
    runtime: MMseqsRuntime,
    uniref90_sha256: str,
) -> dict[str, object]:
    if runtime.version_token is None or runtime.executable_sha256 is None:
        raise ValueError("Cluster cache requires an exact MMseqs2 version and executable SHA-256")
    return cluster_cache_contract_from_values(
        uniref90_release=config.release_uniprot,
        uniref90_sha256=uniref90_sha256,
        identity=config.identity,
        coverage=config.coverage,
        cov_mode=config.cov_mode,
        cluster_mode=config.cluster_mode,
        alignment_mode=config.alignment_mode,
        cluster_reassign=config.cluster_reassign,
        sensitivity=config.sensitivity,
        evalue=config.evalue,
        expected_mmseqs_version=config.expected_mmseqs_version or "",
        observed_mmseqs_version=runtime.version_token,
        mmseqs_executable_sha256=runtime.executable_sha256,
    )


def cluster_cache_directory(cache_root: Path, contract: dict[str, object]) -> Path:
    method = contract.get("method")
    inputs = contract.get("input")
    if not isinstance(method, dict) or not isinstance(inputs, dict):
        raise ValueError("Cluster cache contract is malformed")
    identity = float(str(method["identity_fraction"]))
    release = str(inputs["uniref90_release"])
    digest = _canonical_hash(contract)
    return (
        cluster_cache_root(cache_root)
        / f"uniref90_{release}"
        / _identity_label(identity)
        / f"contract_{digest[:16]}"
    )


def prepare_cluster_cache_destination(
    cache_root_path: Path, contract: dict[str, object]
) -> Path:
    initialize_cluster_cache_root(cache_root_path)
    output = cluster_cache_directory(cache_root_path, contract)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        return output
    probe = output.parent / f".write-probe-{uuid.uuid4().hex}"
    try:
        with probe.open("x", encoding="ascii") as handle:
            handle.write("cluster-cache-write-probe\n")
            handle.flush()
            os.fsync(handle.fileno())
    finally:
        probe.unlink(missing_ok=True)
    return output


def _directory_files(root: Path) -> list[Path]:
    return sorted(
        path for path in root.rglob("*")
        if path.is_file() and path != root / CACHE_MARKER
    )


def inspect_cluster_cache(
    path: Path,
    *,
    expected_contract: dict[str, object] | None = None,
    verify_file_hashes: bool = True,
) -> dict[str, object]:
    root = path.expanduser().resolve()
    if root.is_file() and root.name == CACHE_MARKER:
        root = root.parent
    marker = root / CACHE_MARKER
    if not root.is_dir() or not marker.is_file():
        raise ValueError(f"Cluster cache is incomplete: {root}")
    payload = json.loads(marker.read_text(encoding="utf-8"))
    if (
        payload.get("schema_name") != CACHE_SCHEMA_NAME
        or payload.get("schema_version") != CACHE_SCHEMA_VERSION
        or payload.get("complete") is not True
    ):
        raise ValueError("Cluster cache marker has an unsupported contract")
    contract = payload.get("contract")
    if not isinstance(contract, dict):
        raise ValueError("Cluster cache marker lacks its scientific contract")
    if payload.get("contract_sha256") != _canonical_hash(contract):
        raise ValueError("Cluster cache contract digest is invalid")
    if expected_contract is not None and contract != expected_contract:
        raise ValueError("Cluster cache scientific contract does not match this run")
    counts = payload.get("counts")
    if (
        not isinstance(counts, dict)
        or not isinstance(counts.get("members"), int)
        or not isinstance(counts.get("clusters"), int)
        or int(counts["members"]) < 1
        or int(counts["clusters"]) < 1
    ):
        raise ValueError("Cluster cache has invalid member or cluster counts")
    files = payload.get("files")
    if not isinstance(files, list) or not files:
        raise ValueError("Cluster cache has no file manifest")
    expected_paths: set[str] = set()
    for entry in files:
        if not isinstance(entry, dict):
            raise ValueError("Cluster cache file entry is malformed")
        relative = str(entry.get("path", ""))
        candidate = (root / relative).resolve()
        if not relative or root not in candidate.parents or not candidate.is_file():
            raise ValueError(f"Cluster cache file is missing: {relative}")
        expected_paths.add(relative)
        if candidate.stat().st_size != entry.get("size_bytes"):
            raise ValueError(f"Cluster cache file-size mismatch: {relative}")
        if verify_file_hashes and sha256_file(candidate) != entry.get("sha256"):
            raise ValueError(f"Cluster cache file hash mismatch: {relative}")
    observed_paths = {
        item.relative_to(root).as_posix() for item in _directory_files(root)
    }
    if observed_paths != expected_paths:
        raise ValueError(
            "Cluster cache files do not reconcile with its marker: "
            f"missing={sorted(expected_paths - observed_paths)}, "
            f"extra={sorted(observed_paths - expected_paths)}"
        )
    if ASSIGNMENTS_FILE not in expected_paths or COMMANDS_FILE not in expected_paths:
        raise ValueError("Cluster cache lacks assignments or its command manifest")
    return payload


def load_cluster_cache(
    cache_root_path: Path,
    contract: dict[str, object],
) -> LoadedClusterCache:
    inspect_cluster_cache_root(cache_root_path)
    root = cluster_cache_directory(cache_root_path, contract)
    payload = inspect_cluster_cache(
        root, expected_contract=contract, verify_file_hashes=True
    )
    return LoadedClusterCache(
        root=root,
        marker_sha256=sha256_file(root / CACHE_MARKER),
        payload=payload,
        assignments=root / ASSIGNMENTS_FILE,
    )


def _copy_assignments(source: Path, destination: Path, *, has_header: bool) -> None:
    with open_text(source) as incoming, destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as outgoing:
                first_content = True
                for line_number, raw_line in enumerate(incoming, start=1):
                    if not raw_line.strip():
                        continue
                    columns = raw_line.rstrip("\r\n").split("\t")
                    if has_header and first_content:
                        first_content = False
                        if columns != ["mmseqs_cluster_id", "uniref90_id"]:
                            raise ValueError(
                                "Published cluster membership has an unexpected header"
                            )
                        continue
                    first_content = False
                    if len(columns) != 2 or not columns[0] or not columns[1]:
                        raise ValueError(
                            f"Malformed cluster assignment at {source}:{line_number}"
                        )
                    outgoing.write(f"{columns[0]}\t{columns[1]}\n")


def _write_canonical_assignments(
    cluster_index: ClusterIndex, destination: Path
) -> None:
    with destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as outgoing:
                outgoing.write("mmseqs_cluster_id\tuniref90_id\n")
                for cluster_id, member_id in cluster_index.iter_assignments():
                    outgoing.write(f"{cluster_id}\t{member_id}\n")


def publish_cluster_cache(
    cache_root_path: Path,
    contract: dict[str, object],
    cluster_index: ClusterIndex,
    command_manifest: Path,
    *,
    producer: dict[str, object],
    log_dir: Path | None = None,
) -> LoadedClusterCache:
    initialize_cluster_cache_root(cache_root_path)
    output = cluster_cache_directory(cache_root_path, contract)
    if output.exists():
        return load_cluster_cache(cache_root_path, contract)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = output.parent / f".{output.name}.staging-{uuid.uuid4().hex}"
    started = time.monotonic()
    try:
        stage.mkdir()
        _write_canonical_assignments(cluster_index, stage / ASSIGNMENTS_FILE)
        shutil.copy2(command_manifest, stage / COMMANDS_FILE)
        if log_dir is not None and log_dir.is_dir():
            shutil.copytree(log_dir, stage / "logs")
        files = [
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
            "contract": contract,
            "contract_sha256": _canonical_hash(contract),
            "assignment_sha256": sha256_file(stage / ASSIGNMENTS_FILE),
            "counts": {
                "members": cluster_index.member_count(),
                "clusters": cluster_index.cluster_count(),
            },
            "producer": producer,
            "files": files,
            "note": (
                "This cache ends immediately after MMseqs2 createtsv and complete assignment "
                "validation. GOA retention, split allocation, labels, term universes, training "
                "population, and PFP exports are deliberately excluded."
            ),
        }
        _json(stage / CACHE_MARKER, payload)
        inspect_cluster_cache(
            stage, expected_contract=contract, verify_file_hashes=True
        )
        try:
            stage.rename(output)
        except OSError:
            if not output.exists():
                raise
            existing = load_cluster_cache(cache_root_path, contract)
            if existing.payload.get("assignment_sha256") != payload["assignment_sha256"]:
                raise ValueError(
                    "Concurrent cluster-cache publication produced different assignments"
                )
            shutil.rmtree(stage, ignore_errors=True)
            return existing
        return load_cluster_cache(cache_root_path, contract)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise


def verify_loaded_cluster_cache(cache: LoadedClusterCache) -> None:
    payload = inspect_cluster_cache(
        cache.root,
        expected_contract=cache.payload["contract"],  # type: ignore[arg-type]
        verify_file_hashes=True,
    )
    if sha256_file(cache.root / CACHE_MARKER) != cache.marker_sha256:
        raise ValueError("Cluster cache marker changed while the run was in progress")
    if payload != cache.payload:
        raise ValueError("Cluster cache contract changed while the run was in progress")


def import_publication_cluster_cache(
    run_dir: Path,
    common_preprocessing_cache: Path,
    cache_root_path: Path,
    work_dir: Path,
) -> LoadedClusterCache:
    from .common_cache import (  # Imported lazily to keep the pipeline import acyclic.
        UNIREF_INDEX_FILE,
        inspect_common_preprocessing_cache,
    )
    from .pipeline import validate_publication

    publication_root = run_dir.expanduser().resolve()
    validate_publication(publication_root)
    parameters = json.loads(
        (publication_root / "parameters.json").read_text(encoding="utf-8")
    )
    publication = json.loads(
        (publication_root / "publication_metadata.json").read_text(encoding="utf-8")
    )
    input_manifest = json.loads(
        (publication_root / "input_manifest.json").read_text(encoding="utf-8")
    )
    if publication.get("fixture_mode") is True:
        raise ValueError("Fixture publications cannot seed a production cluster cache")
    inputs = input_manifest.get("inputs")
    if not isinstance(inputs, dict) or not isinstance(inputs.get("uniref90_fasta"), dict):
        raise ValueError("Publication lacks its resolved UniRef90 input binding")
    uniref_binding = inputs["uniref90_fasta"]
    uniref_sha256 = str(uniref_binding.get("sha256", ""))

    common_root = common_preprocessing_cache.expanduser().resolve()
    if common_root.is_file():
        common_root = common_root.parent
    common_payload = inspect_common_preprocessing_cache(
        common_root, verify_file_hashes=True
    )
    common_bindings = common_payload.get("input_bindings")
    if (
        not isinstance(common_bindings, dict)
        or not isinstance(common_bindings.get("uniref90_fasta"), dict)
        or common_bindings["uniref90_fasta"].get("sha256") != uniref_sha256
    ):
        raise ValueError(
            "Common preprocessing cache and publication use different UniRef90 inputs"
        )
    uniref = UniRefIndex(common_root / UNIREF_INDEX_FILE)
    expected_uniref_count = int(common_payload["counts"]["uniref90_entries"])  # type: ignore[index]
    if uniref.count() != expected_uniref_count:
        raise ValueError("Common preprocessing cache UniRef90 index count changed")

    contract = cluster_cache_contract_from_values(
        uniref90_release=str(parameters["uniprot_release"]),
        uniref90_sha256=uniref_sha256,
        identity=float(parameters["identity_fraction"]),
        coverage=float(parameters["coverage"]),
        cov_mode=int(parameters["cov_mode"]),
        cluster_mode=int(parameters["cluster_mode"]),
        alignment_mode=int(parameters["alignment_mode"]),
        cluster_reassign=int(parameters["cluster_reassign"]),
        sensitivity=float(parameters["sensitivity"]),
        evalue=float(parameters["evalue"]),
        expected_mmseqs_version=str(publication["expected_mmseqs_version"]),
        observed_mmseqs_version=str(publication["observed_mmseqs_version"]),
        mmseqs_executable_sha256=str(publication["mmseqs_executable_sha256"]),
    )
    work_root = work_dir.expanduser().resolve() / f"cluster-cache-import-{uuid.uuid4().hex}"
    work_root.mkdir(parents=True)
    try:
        normalized = work_root / ASSIGNMENTS_FILE
        _copy_assignments(
            publication_root / "mmseqs_cluster_membership.tsv.gz",
            normalized,
            has_header=True,
        )
        cluster_index = ClusterIndex.build(
            normalized, uniref, work_root / "clusters.sqlite"
        )
        return publish_cluster_cache(
            cache_root_path,
            contract,
            cluster_index,
            publication_root / COMMANDS_FILE,
            producer={
                "imported_from_publication": {
                    "run_id": publication.get("run_id"),
                    "identity_percent": publication.get("identity_percent"),
                    "framework_revision": publication.get("framework_revision"),
                    "repository_commit": publication.get("repository_commit"),
                },
                "publication_output_manifest_sha256": sha256_file(
                    publication_root / "output_manifest.json"
                ),
                "framework_revision": publication.get("framework_revision"),
                "repository_commit": publication.get("repository_commit"),
                "run_id": publication.get("run_id"),
                "benchmark_scope": publication.get("benchmark_scope"),
                "threads": publication.get("mmseqs_threads"),
                "requested_slots": publication.get("requested_slots"),
                "allocated_slots": publication.get("allocated_slots"),
            },
            log_dir=(
                publication_root / "logs" / "mmseqs"
                if (publication_root / "logs" / "mmseqs").is_dir() else None
            ),
        )
    finally:
        shutil.rmtree(work_root, ignore_errors=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Initialize, verify, or import validated MMseqs2 cluster caches"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    initialize = subparsers.add_parser("init-root")
    initialize.add_argument("--cache-root", type=Path, required=True)
    verify_root = subparsers.add_parser("verify-root")
    verify_root.add_argument("--cache-root", type=Path, required=True)
    verify = subparsers.add_parser("verify")
    verify.add_argument("--cache-dir", type=Path, required=True)
    verify.add_argument("--full-hashes", action="store_true")
    import_run = subparsers.add_parser("import-publication")
    import_run.add_argument("--run-dir", type=Path, required=True)
    import_run.add_argument("--common-preprocessing-cache", type=Path, required=True)
    import_run.add_argument("--cache-root", type=Path, required=True)
    import_run.add_argument("--work-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.command == "init-root":
        marker = initialize_cluster_cache_root(args.cache_root)
        result: object = {"cache_root_marker": str(marker), "status": "ready"}
    elif args.command == "verify-root":
        result = inspect_cluster_cache_root(args.cache_root)
    elif args.command == "verify":
        result = inspect_cluster_cache(
            args.cache_dir, verify_file_hashes=args.full_hashes
        )
    else:
        cache = import_publication_cluster_cache(
            args.run_dir,
            args.common_preprocessing_cache,
            args.cache_root,
            args.work_dir,
        )
        result = {
            "cache_directory": str(cache.root),
            "marker_sha256": cache.marker_sha256,
            "counts": cache.payload["counts"],
            "contract_sha256": cache.payload["contract_sha256"],
        }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
