from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

from . import __version__
from .benchmark import verify_input_identities
from .models import (
    BenchmarkData,
    EmbeddedProtein,
    PlanRecord,
    REGENERATE_MODALITIES,
    ReusePlan,
    ReusePlannerError,
)
from .planner import validate_plan


ACTION_TSV_COLUMNS: Tuple[str, ...] = (
    "protein_id",
    "sequence",
    "sequence_sha256",
    "action",
    "reason",
    "matching_embedded_benchmarks",
    "embedded_benchmark_memberships",
    "target_memberships",
    "regenerate_modalities",
)
KNOWN_TSV_COLUMNS: Tuple[str, ...] = (
    "protein_id",
    "sequence",
    "sequence_sha256",
    "embedded_benchmarks",
    "embedded_benchmark_memberships",
)
PAYLOAD_PATHS: Tuple[str, ...] = (
    "known_embedded_proteins.tsv",
    "regenerate_proteins.fasta",
    "regenerate_proteins.tsv",
    "regenerate_proteins.txt",
    "reuse_proteins.tsv",
    "reuse_proteins.txt",
    "run_manifest.json",
    "summary.json",
    "summary.md",
)


class ReportError(ReusePlannerError):
    pass


def write_reports(plan: ReusePlan, output_dir: Path) -> Path:
    validate_plan(plan)
    output = _absolute_output_path(output_dir)
    if os.path.lexists(output):
        raise ReportError("Refusing to overwrite existing output path: %s" % output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if not output.parent.is_dir():
        raise ReportError("Output parent is not a directory: %s" % output.parent)

    stage = Path(
        tempfile.mkdtemp(prefix=".%s.staging-" % output.name, dir=str(output.parent))
    )
    try:
        _write_payloads(plan, stage, output)
        _validate_staged_payloads(plan, stage)
        verify_input_identities((*plan.embedded_benchmarks, plan.target_benchmark))

        output_manifest = _build_output_manifest(stage)
        _write_json(stage / "output_manifest.json", output_manifest)
        _validate_output_manifest(stage)

        manifest_identity = _output_identity(stage / "output_manifest.json", stage)
        completion = {
            "complete": True,
            "counts": _counts(plan),
            "output_manifest": manifest_identity,
            "schema_version": 1,
        }
        _write_json(stage / "RUN_COMPLETE.json", completion)
        _validate_completion_marker(stage)
        verify_input_identities((*plan.embedded_benchmarks, plan.target_benchmark))

        if os.path.lexists(output):
            raise ReportError("Output path appeared during planning: %s" % output)
        os.rename(stage, output)
    finally:
        if os.path.lexists(stage):
            shutil.rmtree(stage, ignore_errors=True)
    return output


def canonical_command_arguments(plan: ReusePlan, output_dir: Path) -> List[str]:
    arguments = ["plan"]
    for benchmark in plan.embedded_benchmarks:
        arguments.extend(
            ["--embedded-benchmark", "%s=%s" % (benchmark.name, benchmark.directory)]
        )
    arguments.extend(
        [
            "--target-benchmark",
            "%s=%s" % (plan.target_benchmark.name, plan.target_benchmark.directory),
            "--output-dir",
            str(output_dir),
        ]
    )
    return arguments


def _write_payloads(plan: ReusePlan, stage: Path, output: Path) -> None:
    _write_action_tsv(stage / "reuse_proteins.tsv", plan.reuse_records)
    _write_action_tsv(stage / "regenerate_proteins.tsv", plan.regenerate_records)
    _write_id_text(stage / "reuse_proteins.txt", plan.reuse_records)
    _write_id_text(stage / "regenerate_proteins.txt", plan.regenerate_records)
    _write_fasta(stage / "regenerate_proteins.fasta", plan.regenerate_records)
    _write_known_tsv(stage / "known_embedded_proteins.tsv", plan.known_embedded_proteins)

    summary = _summary(plan)
    _write_json(stage / "summary.json", summary)
    _write_summary_markdown(stage / "summary.md", plan)
    _write_json(stage / "run_manifest.json", _run_manifest(plan, output))


def _write_action_tsv(path: Path, records: Sequence[PlanRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=ACTION_TSV_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(_plan_row(record))


def _write_known_tsv(path: Path, records: Sequence[EmbeddedProtein]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=KNOWN_TSV_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(_known_row(record))


def _write_id_text(path: Path, records: Sequence[PlanRecord]) -> None:
    path.write_text(
        "".join("%s\n" % record.protein_id for record in records),
        encoding="utf-8",
    )


def _write_fasta(path: Path, records: Sequence[PlanRecord]) -> None:
    path.write_text(
        "".join(">%s\n%s\n" % (record.protein_id, record.sequence) for record in records),
        encoding="utf-8",
    )


def _write_summary_markdown(path: Path, plan: ReusePlan) -> None:
    counts = _counts(plan)
    references = ", ".join(benchmark.name for benchmark in plan.embedded_benchmarks)
    lines = [
        "# Benchmark embedding reuse plan",
        "",
        "Comparison key: exact, case-sensitive `(protein ID, complete sequence)`.",
        "GO labels and ontology/split memberships are recorded but do not affect the action.",
        "",
        "## Benchmarks",
        "",
        "- Embedded benchmarks: %s" % references,
        "- Target benchmark: %s" % plan.target_benchmark.name,
        "",
        "## Counts",
        "",
        "| Measure | Count |",
        "|---|---:|",
        "| Known embedded proteins | %d |" % counts["known_embedded_proteins"],
        "| Target proteins | %d |" % counts["target_proteins"],
        "| Reuse | %d |" % counts["reuse"],
        "| Regenerate | %d |" % counts["regenerate"],
        "",
        "## Action semantics",
        "",
        "- `reuse`: do not regenerate the protein; preserve prior files and prior missing/masked behaviour.",
        "- `regenerate`: run all four modalities: prott5, text, structure, and ppi.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def _summary(plan: ReusePlan) -> Dict[str, Any]:
    return {
        "action_semantics": {
            "regenerate": (
                "Run all four PFP embedding-generation modalities for this protein."
            ),
            "reuse": (
                "Do not regenerate this protein; preserve prior embedding files and "
                "missing or masked modality behaviour."
            ),
        },
        "comparison_policy": {
            "case_sensitive": True,
            "comparison_key": ["protein_id", "complete_sequence"],
            "go_labels_affect_action": False,
            "memberships_affect_action": False,
            "sequence_only_matching": False,
        },
        "counts": _counts(plan),
        "embedded_benchmarks": [
            {"name": benchmark.name, "path": str(benchmark.directory)}
            for benchmark in plan.embedded_benchmarks
        ],
        "regenerate_modalities": list(REGENERATE_MODALITIES),
        "schema_version": 1,
        "target_benchmark": {
            "name": plan.target_benchmark.name,
            "path": str(plan.target_benchmark.directory),
        },
    }


def _run_manifest(plan: ReusePlan, output: Path) -> Dict[str, Any]:
    return {
        "benchmarks": {
            "embedded": [
                _benchmark_manifest(benchmark) for benchmark in plan.embedded_benchmarks
            ],
            "target": _benchmark_manifest(plan.target_benchmark),
        },
        "command": [
            "python",
            "-m",
            "pfp_benchmark_reuse",
            *canonical_command_arguments(plan, output),
        ],
        "command_arguments": canonical_command_arguments(plan, output),
        "counts": _counts(plan),
        "output_dir": str(output),
        "schema_version": 1,
        "tool": {"name": "pfp-benchmark-reuse", "version": __version__},
    }


def _benchmark_manifest(benchmark: BenchmarkData) -> Dict[str, Any]:
    return {
        "duplicate_occurrences": benchmark.duplicate_occurrences,
        "input_csvs": [
            {
                "relative_path": identity.relative_path,
                "resolved_path": str(identity.resolved_path),
                "sha256": identity.sha256,
                "size_bytes": identity.size_bytes,
            }
            for identity in benchmark.input_files
        ],
        "name": benchmark.name,
        "path": str(benchmark.directory),
        "protein_count": len(benchmark.proteins),
    }


def _counts(plan: ReusePlan) -> Dict[str, int]:
    return {
        "embedded_benchmarks": len(plan.embedded_benchmarks),
        "known_embedded_proteins": len(plan.known_embedded_proteins),
        "regenerate": len(plan.regenerate_records),
        "reuse": len(plan.reuse_records),
        "target_proteins": len(plan.target_benchmark.proteins),
    }


def _plan_row(record: PlanRecord) -> Dict[str, str]:
    return {
        "action": record.action,
        "embedded_benchmark_memberships": _json_list(
            record.embedded_benchmark_memberships
        ),
        "matching_embedded_benchmarks": _json_list(
            record.matching_embedded_benchmarks
        ),
        "protein_id": record.protein_id,
        "reason": record.reason,
        "regenerate_modalities": _json_list(record.regenerate_modalities),
        "sequence": record.sequence,
        "sequence_sha256": record.sequence_sha256,
        "target_memberships": _json_list(record.target_memberships),
    }


def _known_row(record: EmbeddedProtein) -> Dict[str, str]:
    return {
        "embedded_benchmark_memberships": _json_list(
            record.embedded_benchmark_memberships
        ),
        "embedded_benchmarks": _json_list(record.embedded_benchmarks),
        "protein_id": record.protein_id,
        "sequence": record.sequence,
        "sequence_sha256": record.sequence_sha256,
    }


def _json_list(values: Iterable[str]) -> str:
    return json.dumps(list(values), ensure_ascii=False, separators=(",", ":"))


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _build_output_manifest(stage: Path) -> Dict[str, Any]:
    files = [_output_identity(stage / relative_path, stage) for relative_path in PAYLOAD_PATHS]
    return {"file_count": len(files), "files": files, "schema_version": 1}


def _output_identity(path: Path, root: Path) -> Dict[str, Any]:
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _validate_staged_payloads(plan: ReusePlan, stage: Path) -> None:
    expected_by_action = {
        "reuse": [_plan_row(record) for record in plan.reuse_records],
        "regenerate": [_plan_row(record) for record in plan.regenerate_records],
    }
    observed_by_action = {
        "reuse": _read_tsv(stage / "reuse_proteins.tsv", ACTION_TSV_COLUMNS),
        "regenerate": _read_tsv(
            stage / "regenerate_proteins.tsv", ACTION_TSV_COLUMNS
        ),
    }
    if observed_by_action != expected_by_action:
        raise ReportError("Written action TSVs do not match the validated plan")

    reuse_ids = [row["protein_id"] for row in observed_by_action["reuse"]]
    regenerate_ids = [row["protein_id"] for row in observed_by_action["regenerate"]]
    target_ids = sorted(plan.target_benchmark.proteins)
    if set(reuse_ids) & set(regenerate_ids) or sorted(reuse_ids + regenerate_ids) != target_ids:
        raise ReportError("Written action TSVs are not an exact binary target partition")
    if any(row["action"] != "reuse" for row in observed_by_action["reuse"]):
        raise ReportError("reuse_proteins.tsv contains a non-reuse action")
    if any(row["action"] != "regenerate" for row in observed_by_action["regenerate"]):
        raise ReportError("regenerate_proteins.tsv contains a non-regenerate action")

    if (stage / "reuse_proteins.txt").read_text(encoding="utf-8").splitlines() != reuse_ids:
        raise ReportError("reuse_proteins.txt does not match reuse_proteins.tsv")
    if (
        stage / "regenerate_proteins.txt"
    ).read_text(encoding="utf-8").splitlines() != regenerate_ids:
        raise ReportError("regenerate_proteins.txt does not match regenerate_proteins.tsv")

    fasta = _read_fasta(stage / "regenerate_proteins.fasta")
    expected_fasta = [
        (
            record.protein_id,
            plan.target_benchmark.proteins[record.protein_id].sequence,
        )
        for record in plan.regenerate_records
    ]
    if fasta != expected_fasta:
        raise ReportError("Regenerate FASTA does not match target IDs and sequences")

    expected_known = [_known_row(record) for record in plan.known_embedded_proteins]
    if _read_tsv(stage / "known_embedded_proteins.tsv", KNOWN_TSV_COLUMNS) != expected_known:
        raise ReportError("known_embedded_proteins.tsv does not match the embedded union")

    summary = _read_json(stage / "summary.json")
    run_manifest = _read_json(stage / "run_manifest.json")
    if summary.get("counts") != _counts(plan) or run_manifest.get("counts") != _counts(plan):
        raise ReportError("JSON report counts do not match the validated plan")


def _validate_output_manifest(stage: Path) -> None:
    manifest = _read_json(stage / "output_manifest.json")
    files = manifest.get("files")
    if not isinstance(files, list):
        raise ReportError("output_manifest.json has no file list")
    if manifest.get("file_count") != len(PAYLOAD_PATHS) or len(files) != len(PAYLOAD_PATHS):
        raise ReportError("output_manifest.json has an incorrect file count")
    observed_paths = [item.get("path") for item in files if isinstance(item, dict)]
    if observed_paths != list(PAYLOAD_PATHS):
        raise ReportError("output_manifest.json has an incorrect or unsorted path set")
    for item in files:
        path = stage / item["path"]
        if not path.is_file():
            raise ReportError("Manifested output is missing: %s" % item["path"])
        if _output_identity(path, stage) != item:
            raise ReportError("Manifest identity does not match: %s" % item["path"])


def _validate_completion_marker(stage: Path) -> None:
    completion = _read_json(stage / "RUN_COMPLETE.json")
    if completion.get("complete") is not True:
        raise ReportError("RUN_COMPLETE.json does not mark the run complete")
    expected = _output_identity(stage / "output_manifest.json", stage)
    if completion.get("output_manifest") != expected:
        raise ReportError("RUN_COMPLETE.json does not authenticate output_manifest.json")


def _read_tsv(path: Path, expected_columns: Sequence[str]) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if reader.fieldnames != list(expected_columns):
            raise ReportError("Unexpected TSV columns in %s" % path.name)
        return [dict(row) for row in reader]


def _read_fasta(path: Path) -> List[Tuple[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) % 2:
        raise ReportError("Regenerate FASTA has an incomplete record")
    records: List[Tuple[str, str]] = []
    for index in range(0, len(lines), 2):
        header, sequence = lines[index : index + 2]
        if not header.startswith(">") or not header[1:] or not sequence:
            raise ReportError("Regenerate FASTA contains a malformed record")
        records.append((header[1:], sequence))
    if len(records) != len({protein_id for protein_id, _ in records}):
        raise ReportError("Regenerate FASTA contains duplicate protein IDs")
    return records


def _read_json(path: Path) -> Dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ReportError("Invalid JSON output %s: %s" % (path.name, exc)) from exc
    if not isinstance(payload, dict):
        raise ReportError("JSON output must contain an object: %s" % path.name)
    return payload


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _absolute_output_path(output_dir: Path) -> Path:
    expanded = output_dir.expanduser()
    if not expanded.is_absolute():
        expanded = Path.cwd() / expanded
    absolute = Path(os.path.abspath(str(expanded)))
    if os.path.lexists(absolute):
        return absolute
    return absolute.parent.resolve() / absolute.name
