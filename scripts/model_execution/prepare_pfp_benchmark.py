#!/usr/bin/env python3
"""Validate nine CAFA-style CSVs and materialize PFP training artifacts.

The upstream PFP preparer is used without modification. This wrapper supplies
the fail-closed input and post-materialization checks that it does not provide.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tarfile
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

import numpy as np
import scipy.sparse as ssp

from common import (
    ASPECTS,
    ASPECT_TO_CSV,
    ASPECT_TO_NAMESPACE,
    CSV_SPLITS,
    FRAMEWORK_ROOT,
    PFP_SPLITS,
    atomic_write_json,
    load_run_config,
    require_empty_output,
    sha256_file,
    sha256_lines,
)


INVENTORY_SRC = FRAMEWORK_ROOT / "embedding_inventory" / "src"
HOMOLOGY_SRC = FRAMEWORK_ROOT / "benchmark_builders" / "homology_cluster" / "src"
sys.path.insert(0, str(INVENTORY_SRC))
from pfp_embedding_inventory.benchmark import (  # noqa: E402
    BenchmarkError,
    parse_benchmark,
    required_csv_names,
)
from pfp_embedding_inventory.models import BenchmarkContract, BenchmarkData  # noqa: E402


def read_obo(path: Path) -> Tuple[Dict[str, str], set[str]]:
    """Return live GO term namespaces and obsolete term IDs from an OBO file."""
    namespaces: Dict[str, str] = {}
    obsolete: set[str] = set()
    stanza: Dict[str, Any] = {}

    def publish(value: Dict[str, Any]) -> None:
        term_id = value.get("id")
        namespace = value.get("namespace")
        if not term_id or not str(term_id).startswith("GO:"):
            return
        if value.get("is_obsolete") == "true":
            obsolete.add(str(term_id))
            return
        if namespace:
            namespaces[str(term_id)] = str(namespace)

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                publish(stanza)
                stanza = {"type": "Term"}
            elif line.startswith("["):
                publish(stanza)
                stanza = {}
            elif stanza.get("type") == "Term" and ": " in line:
                key, value = line.split(": ", 1)
                if key in {"id", "namespace", "is_obsolete"}:
                    stanza[key] = value
        publish(stanza)
    if not namespaces:
        raise ValueError(f"No live GO terms were parsed from {path}")
    return namespaces, obsolete


def csv_profile(
    benchmark_dir: Path,
    allow_all_zero_rows: bool,
) -> Tuple[Dict[str, Any], Dict[str, list[str]]]:
    profiles: Dict[str, Any] = {}
    ordered_terms: Dict[str, list[str]] = {}
    for aspect in ASPECTS:
        prefix = ASPECT_TO_CSV[aspect]
        aspect_headers: Dict[str, list[str]] = {}
        for split in CSV_SPLITS:
            path = benchmark_dir / f"{prefix}-{split}.csv"
            with path.open("r", newline="", encoding="utf-8-sig") as handle:
                reader = csv.reader(handle, strict=True)
                try:
                    header = next(reader)
                except StopIteration as exc:
                    raise ValueError(f"Required CSV is empty: {path}") from exc
                terms = header[2:]
                row_count = 0
                positive_rows = 0
                all_zero_rows = 0
                term_positives: Counter[str] = Counter()
                for line_number, row in enumerate(reader, start=2):
                    row_count += 1
                    labels = row[2:]
                    positive_indices = [index for index, value in enumerate(labels) if value == "1"]
                    if positive_indices:
                        positive_rows += 1
                        for index in positive_indices:
                            term_positives[terms[index]] += 1
                    else:
                        all_zero_rows += 1
                if row_count == 0:
                    raise ValueError(f"Benchmark split contains no rows: {path.name}")
                if all_zero_rows and not allow_all_zero_rows:
                    raise ValueError(
                        f"{path.name} contains {all_zero_rows} all-zero label rows; "
                        "this run policy does not permit unlabeled proteins"
                    )
                aspect_headers[split] = terms
                profiles[path.name] = {
                    "rows": row_count,
                    "positive_rows": positive_rows,
                    "all_zero_rows": all_zero_rows,
                    "go_terms": len(terms),
                    "ordered_term_sha256": sha256_lines(terms),
                    "sha256": sha256_file(path),
                    "zero_positive_terms": [
                        term for term in terms if term_positives[term] == 0
                    ],
                }
        training_terms = aspect_headers["training"]
        for split in ("validation", "test"):
            if aspect_headers[split] != training_terms:
                raise ValueError(
                    f"{aspect} GO columns differ in identity or order between "
                    f"training and {split}"
                )
        ordered_terms[aspect] = training_terms
    return profiles, ordered_terms


def global_cross_split_diagnostics(benchmark: BenchmarkData) -> Dict[str, int]:
    """Count global split overlap without changing the configured fail policy."""
    values: Dict[str, Tuple[set[str], set[str]]] = {}
    for split in ("training", "validation", "test"):
        protein_ids = set().union(
            *(benchmark.file_members[(ontology, split)] for ontology in ("BP", "CC", "MF"))
        )
        sequence_hashes = {
            benchmark.proteins[protein_id].sequence_sha256
            for protein_id in protein_ids
        }
        values[split] = (protein_ids, sequence_hashes)

    diagnostics: Dict[str, int] = {}
    for left, right in (
        ("training", "validation"),
        ("training", "test"),
        ("validation", "test"),
    ):
        key = f"{left}_{right}"
        diagnostics[f"global_{key}_protein_ids"] = len(
            values[left][0] & values[right][0]
        )
        diagnostics[f"global_{key}_exact_sequences"] = len(
            values[left][1] & values[right][1]
        )
    return diagnostics


def normalize_legacy_headers(
    source: Path, destination: Path, allow_singular: bool
) -> tuple[Path, list[Dict[str, str]]]:
    aliases: list[Dict[str, str]] = []
    singular_files: set[str] = set()
    for filename in required_csv_names():
        path = source / filename
        with path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle, strict=True)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValueError(f"Required CSV is empty: {path}") from exc
        if header and header[0] == "protein":
            singular_files.add(filename)
        elif not header or header[0] != "proteins":
            raise ValueError(
                f"{filename} must begin with proteins,sequences; found {header[:2]}"
            )
    if not singular_files:
        return source, aliases
    if not allow_singular:
        raise ValueError(
            "Legacy singular 'protein' header found but this run policy does not permit it: "
            + ", ".join(sorted(singular_files))
        )
    if destination.exists():
        raise ValueError(f"Normalized benchmark directory already exists: {destination}")
    destination.mkdir(parents=True)
    for filename in required_csv_names():
        source_path = source / filename
        destination_path = destination / filename
        if filename not in singular_files:
            shutil.copy2(source_path, destination_path)
            continue
        with source_path.open("r", newline="", encoding="utf-8-sig") as source_handle:
            rows = csv.reader(source_handle, strict=True)
            with destination_path.open("w", newline="", encoding="utf-8") as destination_handle:
                writer = csv.writer(destination_handle, lineterminator="\n")
                header = next(rows)
                header[0] = "proteins"
                writer.writerow(header)
                writer.writerows(rows)
        aliases.append(
            {
                "file": filename,
                "source_header": "protein",
                "materialized_header": "proteins",
                "reason": "PFP compatibility alias; source file was not modified",
            }
        )
    return destination, aliases


def validate_ontology_terms(ordered_terms: Dict[str, list[str]], obo_file: Path) -> Dict[str, Any]:
    namespaces, obsolete = read_obo(obo_file)
    result: Dict[str, Any] = {}
    for aspect, terms in ordered_terms.items():
        wanted = ASPECT_TO_NAMESPACE[aspect]
        missing = [term for term in terms if term not in namespaces]
        wrong = [term for term in terms if term in namespaces and namespaces[term] != wanted]
        old = [term for term in terms if term in obsolete]
        if missing or wrong or old:
            examples = {
                "missing": missing[:5],
                "wrong_namespace": wrong[:5],
                "obsolete": old[:5],
            }
            raise ValueError(f"{aspect} GO/OBO contract failed: {examples}")
        result[aspect] = {
            "namespace": wanted,
            "terms": len(terms),
            "ordered_term_sha256": sha256_lines(terms),
        }
    return result


def run_logged(command: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if completed.returncode:
        raise RuntimeError(
            f"Command failed with status {completed.returncode}; see {log_path}: "
            + " ".join(command)
        )


def iter_source_rows(path: Path) -> Iterable[Tuple[str, str, list[int]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, strict=True)
        next(reader)
        for row in reader:
            yield row[0], row[1], [index for index, value in enumerate(row[2:]) if value == "1"]


def verify_materialization(
    benchmark_dir: Path, data_dir: Path, source_terms: Dict[str, list[str]]
) -> Dict[str, Any]:
    summary: Dict[str, Any] = {}
    for aspect in ASPECTS:
        prefix = ASPECT_TO_CSV[aspect]
        terms = json.loads((data_dir / f"{aspect}_go_terms.json").read_text(encoding="utf-8"))
        if terms != source_terms[aspect]:
            raise ValueError(f"Materialized GO term identity/order mismatch for {aspect}")
        aspect_summary: Dict[str, Any] = {
            "terms": len(terms),
            "ordered_term_sha256": sha256_lines(terms),
            "splits": {},
        }
        for csv_split, pfp_split in PFP_SPLITS.items():
            names = [
                str(value)
                for value in np.load(
                    data_dir / f"{aspect}_{pfp_split}_names.npy", allow_pickle=True
                )
            ]
            sequences = json.loads(
                (data_dir / f"{aspect}_{pfp_split}_sequences.json").read_text(
                    encoding="utf-8"
                )
            )
            labels = ssp.load_npz(data_dir / f"{aspect}_{pfp_split}_labels.npz").tocsr()
            if labels.shape != (len(names), len(terms)):
                raise ValueError(f"Materialized shape mismatch for {aspect} {pfp_split}")
            source_count = 0
            for row_index, (protein_id, sequence, positive_indices) in enumerate(
                iter_source_rows(benchmark_dir / f"{prefix}-{csv_split}.csv")
            ):
                source_count += 1
                if row_index >= len(names):
                    raise ValueError(f"Materialized row count mismatch for {aspect} {pfp_split}")
                if names[row_index] != protein_id:
                    raise ValueError(f"Materialized ID/order mismatch for {aspect} {pfp_split}")
                if sequences.get(protein_id) != sequence:
                    raise ValueError(f"Materialized sequence mismatch for {protein_id}")
                observed = labels.indices[labels.indptr[row_index] : labels.indptr[row_index + 1]].tolist()
                if observed != positive_indices:
                    raise ValueError(f"Materialized labels mismatch for {aspect}/{protein_id}")
            if source_count != len(names):
                raise ValueError(f"Materialized row count mismatch for {aspect} {pfp_split}")
            aspect_summary["splits"][pfp_split] = {"proteins": len(names)}
        summary[aspect] = aspect_summary
    return summary


def prepared_reference_names() -> list[str]:
    names = [f"{aspect}_go_terms.json" for aspect in ASPECTS]
    for aspect in ASPECTS:
        for split in PFP_SPLITS.values():
            stem = f"{aspect}_{split}"
            names.extend(
                f"{stem}_{suffix}"
                for suffix in ("names.npy", "labels.npz", "sequences.json")
            )
    return sorted(names)


def file_digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bind_reference_archive(
    reference: Path,
    archive_path: Path,
    policy: Dict[str, Any],
) -> Dict[str, Any]:
    if not archive_path.is_file():
        raise FileNotFoundError(f"Prepared reference archive is missing: {archive_path}")
    expected_name = str(policy["name"])
    expected_size = int(policy["size_bytes"])
    checksum_algorithm = str(policy["checksum_algorithm"])
    expected_checksum = str(policy["checksum"])
    if archive_path.name != expected_name:
        raise ValueError(
            f"Prepared reference archive must be named {expected_name}: {archive_path}"
        )
    if archive_path.stat().st_size != expected_size:
        raise ValueError("Prepared reference archive size differs from the published artifact")
    observed_checksum = file_digest(archive_path, checksum_algorithm)
    if observed_checksum != expected_checksum:
        raise ValueError(
            "Prepared reference archive checksum differs from the published artifact"
        )

    expected_members = {
        f"data/{name}": name for name in prepared_reference_names()
    }
    observed_members: Dict[str, str] = {}
    with tarfile.open(archive_path, mode="r:gz") as archive:
        for member in archive:
            normalized = member.name.removeprefix("./")
            filename = expected_members.get(normalized)
            if filename is None:
                continue
            if filename in observed_members:
                raise ValueError(f"Prepared reference archive repeats data/{filename}")
            if not member.isfile():
                raise ValueError(
                    f"Prepared reference archive member is not a file: {normalized}"
                )
            source = archive.extractfile(member)
            if source is None:
                raise ValueError(
                    f"Cannot read prepared reference archive member: {normalized}"
                )
            digest = hashlib.sha256()
            with source:
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    digest.update(chunk)
            observed_members[filename] = digest.hexdigest()

    missing = sorted(set(prepared_reference_names()) - set(observed_members))
    if missing:
        raise ValueError(
            "Prepared reference archive lacks required data members: "
            + ", ".join(missing[:5])
        )
    disk_members = {
        name: sha256_file(reference / name) for name in prepared_reference_names()
    }
    if disk_members != observed_members:
        differing = sorted(
            name
            for name in disk_members
            if disk_members[name] != observed_members.get(name)
        )
        raise ValueError(
            "Prepared reference directory was not extracted from the supplied published "
            f"archive; differing={differing[:5]}"
        )
    return {
        "path": str(archive_path.resolve()),
        "name": archive_path.name,
        "size_bytes": archive_path.stat().st_size,
        "checksum_algorithm": checksum_algorithm,
        "checksum": observed_checksum,
        "member_count": len(observed_members),
        "member_sha256": observed_members,
    }


def compare_prepared(
    reference: Path,
    observed: Path,
    source_archive: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    compared = 0
    reference_artifacts: Dict[str, str] = {}
    for aspect in ASPECTS:
        for filename in (f"{aspect}_go_terms.json",):
            if json.loads((reference / filename).read_text()) != json.loads((observed / filename).read_text()):
                raise ValueError(f"Prepared reference differs: {filename}")
            reference_artifacts[filename] = sha256_file(reference / filename)
            compared += 1
        for split in PFP_SPLITS.values():
            stem = f"{aspect}_{split}"
            left_names = np.load(reference / f"{stem}_names.npy", allow_pickle=True)
            right_names = np.load(observed / f"{stem}_names.npy", allow_pickle=True)
            if not np.array_equal(left_names, right_names):
                raise ValueError(f"Prepared reference differs: {stem}_names.npy")
            left_labels = ssp.load_npz(reference / f"{stem}_labels.npz")
            right_labels = ssp.load_npz(observed / f"{stem}_labels.npz")
            if left_labels.shape != right_labels.shape or (left_labels != right_labels).nnz:
                raise ValueError(f"Prepared reference differs: {stem}_labels.npz")
            if json.loads((reference / f"{stem}_sequences.json").read_text()) != json.loads(
                (observed / f"{stem}_sequences.json").read_text()
            ):
                raise ValueError(f"Prepared reference differs: {stem}_sequences.json")
            for suffix in ("names.npy", "labels.npz", "sequences.json"):
                filename = f"{stem}_{suffix}"
                reference_artifacts[filename] = sha256_file(reference / filename)
            compared += 3
    reference_fingerprint = hashlib.sha256(
        json.dumps(
            reference_artifacts,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()
    result = {
        "passed": True,
        "artifacts_compared": compared,
        "reference": str(reference),
        "reference_artifact_sha256": reference_artifacts,
        "reference_fingerprint": reference_fingerprint,
    }
    if source_archive is not None:
        result["published_source_archive"] = source_archive
    return result


def bind_publication_evidence(
    evidence_by_name: Dict[str, Dict[str, str]],
    benchmark_dir: Path,
) -> Dict[str, Any] | None:
    """Bind a homology publication receipt to the exact nine selected CSVs."""
    if "output_manifest.json" not in evidence_by_name:
        return None
    manifest_path = Path(evidence_by_name["output_manifest.json"]["path"])
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    entries = manifest.get("files")
    if manifest.get("schema_version") != 1 or not isinstance(entries, list):
        raise ValueError("Benchmark output_manifest.json has an unsupported schema")
    if manifest.get("payload_file_count") != len(entries):
        raise ValueError("Benchmark output_manifest.json has an inconsistent payload count")
    selected: Dict[str, Dict[str, Any]] = {}
    required = set(required_csv_names())
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Benchmark output_manifest.json contains a non-object entry")
        name = Path(str(entry.get("path", ""))).name
        if name not in required:
            continue
        if name in selected:
            raise ValueError(f"Benchmark output manifest lists {name} more than once")
        selected[name] = entry
    missing = sorted(required - set(selected))
    if missing:
        raise ValueError(
            "Benchmark output manifest does not contain all selected CSVs: "
            + ", ".join(missing)
        )
    bindings = {}
    for name in required_csv_names():
        source = benchmark_dir / name
        entry = selected[name]
        expected_size = entry.get("size_bytes", entry.get("bytes"))
        if entry.get("sha256") != sha256_file(source) or expected_size != source.stat().st_size:
            raise ValueError(f"Benchmark output manifest does not bind selected CSV {name}")
        bindings[name] = {
            "manifest_path": str(entry["path"]),
            "sha256": entry["sha256"],
            "size_bytes": source.stat().st_size,
        }
    receipt_bindings = {}
    for name in ("validation_report.json",):
        if name not in evidence_by_name:
            continue
        matches = [entry for entry in entries if Path(str(entry.get("path", ""))).name == name]
        if len(matches) != 1:
            raise ValueError(f"Benchmark output manifest must list supplied {name} exactly once")
        entry = matches[0]
        source = Path(evidence_by_name[name]["path"])
        expected_size = entry.get("size_bytes", entry.get("bytes"))
        if entry.get("sha256") != sha256_file(source) or expected_size != source.stat().st_size:
            raise ValueError(f"Benchmark output manifest does not bind supplied {name}")
        receipt_bindings[name] = entry["sha256"]
    return {"passed": True, "csvs": bindings, "receipts": receipt_bindings}


def require_production_homology_publication(publication: Dict[str, Any]) -> None:
    """Reject valid fixtures and pilots at the production PFP boundary."""
    if (
        publication.get("benchmark_scope") != "dissertation-production"
        or publication.get("production_eligible") is not True
        or publication.get("fixture_mode") is not False
    ):
        raise ValueError(
            "Homology model execution requires a dissertation-production, "
            "production-eligible, non-fixture publication"
        )


def run_domain_validation(kind: str | None, benchmark_dir: Path) -> Dict[str, Any] | None:
    if kind is None:
        return None
    if kind != "homology-publication":
        raise ValueError(f"Unknown benchmark domain validator: {kind}")
    sys.path.insert(0, str(HOMOLOGY_SRC))
    from homology_cluster_benchmark.pipeline import validate_publication

    validate_publication(benchmark_dir)
    publication = json.loads(
        (benchmark_dir / "publication_metadata.json").read_text(encoding="utf-8")
    )
    require_production_homology_publication(publication)
    return {
        "validator": kind,
        "passed": True,
        "benchmark_scope": publication["benchmark_scope"],
        "production_eligible": publication["production_eligible"],
        "fixture_mode": publication["fixture_mode"],
        "identity_percent": publication.get("identity_percent"),
        "scientific_fingerprint": publication.get("scientific_fingerprint"),
        "publication_metadata_sha256": sha256_file(
            benchmark_dir / "publication_metadata.json"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--obo-file", type=Path, required=True)
    parser.add_argument("--pfp-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--reference-data-dir", type=Path)
    parser.add_argument("--reference-source-archive", type=Path)
    parser.add_argument("--benchmark-evidence", type=Path, action="append", default=[])
    args = parser.parse_args()

    benchmark_dir = args.benchmark_dir.resolve()
    obo_file = args.obo_file.resolve()
    pfp_root = args.pfp_root.resolve()
    config = load_run_config(args.config)
    contract_value = config["benchmark_contract"]
    contract = BenchmarkContract(
        id_overlap=contract_value["id_overlap"],
        sequence_overlap=contract_value["sequence_overlap"],
        protein_id_pattern=contract_value["protein_id_pattern"],
        sequence_pattern=contract_value["sequence_pattern"],
    )
    if not obo_file.is_file():
        raise FileNotFoundError(f"GO OBO file does not exist: {obo_file}")
    upstream = pfp_root / "scripts" / "prepare_cafa3_data.py"
    if not upstream.is_file():
        raise FileNotFoundError(f"PFP preparer does not exist: {upstream}")
    required_csvs = set(required_csv_names())
    extra_csvs = sorted(path.name for path in benchmark_dir.glob("*.csv") if path.name not in required_csvs)
    if extra_csvs:
        raise ValueError(
            "Benchmark directory contains unexpected CSV files: " + ", ".join(extra_csvs)
        )

    materialization_source, header_aliases = normalize_legacy_headers(
        benchmark_dir,
        args.data_dir.parent / "normalized_benchmark",
        bool(contract_value.get("allow_legacy_singular_protein_header", False)),
    )
    benchmark = parse_benchmark(materialization_source, contract)
    if benchmark.duplicate_rows and not contract_value.get("allow_duplicate_rows", False):
        raise BenchmarkError(
            f"Benchmark contains {benchmark.duplicate_rows} identical duplicate rows"
        )
    profiles, ordered_terms = csv_profile(
        materialization_source, bool(contract_value.get("allow_all_zero_rows", False))
    )
    ontology = validate_ontology_terms(ordered_terms, obo_file)
    evidence = []
    for path in args.benchmark_evidence:
        resolved = path.resolve()
        if not resolved.is_file():
            raise FileNotFoundError(f"Benchmark evidence is missing: {resolved}")
        evidence.append(
            {"path": str(resolved), "name": resolved.name, "sha256": sha256_file(resolved)}
        )
    if contract_value.get("require_benchmark_evidence", False) and not evidence:
        raise ValueError(
            "This benchmark policy requires domain-owned validation evidence; "
            "pass one or more --benchmark-evidence files"
        )
    evidence_by_name: Dict[str, Dict[str, str]] = {}
    for item in evidence:
        if item["name"] in evidence_by_name:
            raise ValueError(f"Repeated benchmark evidence filename: {item['name']}")
        evidence_by_name[item["name"]] = item
    required_evidence_names = contract_value.get("required_benchmark_evidence_names", [])
    missing_evidence = [name for name in required_evidence_names if name not in evidence_by_name]
    if missing_evidence:
        raise ValueError(
            "Required benchmark evidence is missing: " + ", ".join(missing_evidence)
        )
    if "validation_report.json" in evidence_by_name:
        validation = json.loads(
            Path(evidence_by_name["validation_report.json"]["path"]).read_text(encoding="utf-8")
        )
        if validation.get("valid") is not True:
            raise ValueError("Benchmark validation_report.json does not declare valid=true")
    if "RUN_COMPLETE.json" in evidence_by_name:
        marker = json.loads(
            Path(evidence_by_name["RUN_COMPLETE.json"]["path"]).read_text(encoding="utf-8")
        )
        if marker.get("complete") is not True:
            raise ValueError("Benchmark RUN_COMPLETE.json does not declare complete=true")
        if "output_manifest.json" in evidence_by_name:
            expected_manifest_hash = marker.get("manifest_sha256")
            observed_manifest_hash = evidence_by_name["output_manifest.json"]["sha256"]
            if expected_manifest_hash != observed_manifest_hash:
                raise ValueError(
                    "Benchmark completion marker does not bind the supplied output_manifest.json"
                )
        if "output_manifest.json" in required_evidence_names:
            fingerprint = marker.get("scientific_fingerprint")
            if not isinstance(fingerprint, str) or not re.fullmatch(r"[0-9a-f]{64}", fingerprint):
                raise ValueError(
                    "Benchmark completion marker lacks a valid scientific_fingerprint"
                )
    publication_binding = bind_publication_evidence(evidence_by_name, benchmark_dir)
    domain_validation = run_domain_validation(
        contract_value.get("domain_validator"), benchmark_dir
    )
    require_empty_output(args.data_dir)

    run_logged(
        [
            sys.executable,
            str(upstream),
            "--cafa3-dir",
            str(materialization_source),
            "--output-dir",
            str(args.data_dir),
        ],
        args.log_dir / "prepare_upstream.log",
    )
    run_logged(
        [
            sys.executable,
            str(FRAMEWORK_ROOT / "scripts" / "verification" / "verify_splits.py"),
            "--data-dir",
            str(args.data_dir),
            "--strict",
        ],
        args.log_dir / "verify_splits.log",
    )
    run_logged(
        [
            sys.executable,
            str(FRAMEWORK_ROOT / "scripts" / "embeddings" / "generate_embeddings_fasta.py"),
            "--data-dir",
            str(args.data_dir),
        ],
        args.log_dir / "generate_fasta.log",
    )
    materialized = verify_materialization(
        materialization_source, args.data_dir, ordered_terms
    )
    reference = None
    if args.reference_data_dir:
        reference_dir = args.reference_data_dir.resolve()
        reference_policy = config.get("reference_preparation")
        archive_binding = None
        if reference_policy:
            if not args.reference_source_archive:
                raise ValueError(
                    "This run policy requires --reference-source-archive to authenticate "
                    "the prepared reference"
                )
            archive_binding = bind_reference_archive(
                reference_dir,
                args.reference_source_archive.resolve(),
                reference_policy,
            )
        elif args.reference_source_archive:
            raise ValueError(
                "--reference-source-archive was supplied but the run config has no "
                "reference_preparation policy"
            )
        reference = compare_prepared(
            reference_dir, args.data_dir.resolve(), archive_binding
        )

    source_csv_hashes = {
        name: sha256_file(benchmark_dir / name) for name in required_csv_names()
    }
    scientific_fingerprint_payload = {
        "csv_sha256": source_csv_hashes,
        "obo_sha256": sha256_file(obo_file),
        "ordered_term_sha256": {
            aspect: sha256_lines(terms) for aspect, terms in ordered_terms.items()
        },
    }
    benchmark_fingerprint = hashlib.sha256(
        json.dumps(
            scientific_fingerprint_payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()
    report = {
        "schema_version": 1,
        "status": "passed",
        "benchmark_dir": str(benchmark_dir),
        "materialization_source": str(materialization_source),
        "benchmark_fingerprint": benchmark_fingerprint,
        "population_fingerprint": benchmark.fingerprint,
        "benchmark_fingerprint_payload": scientific_fingerprint_payload,
        "protein_count": len(benchmark.proteins),
        "required_csvs": required_csv_names(),
        "duplicate_rows": benchmark.duplicate_rows,
        "contract": contract_value,
        "global_cross_split_overlap_diagnostics": (
            global_cross_split_diagnostics(benchmark)
        ),
        "csvs": profiles,
        "source_csv_sha256": source_csv_hashes,
        "source_csv_bytes": {
            name: (benchmark_dir / name).stat().st_size for name in required_csv_names()
        },
        "header_compatibility_aliases": header_aliases,
        "benchmark_evidence": evidence,
        "benchmark_publication_binding": publication_binding,
        "domain_validation": domain_validation,
        "obo": {"path": str(obo_file), "sha256": sha256_file(obo_file)},
        "ontology_contract": ontology,
        "materialization": materialized,
        "reference_equivalence": reference,
        "zero_training_positive_terms_are_retained": True,
    }
    atomic_write_json(args.report, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
