#!/usr/bin/env python3
"""Audit how published homology benchmarks change across identity thresholds."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys


ASPECTS = ("bp", "cc", "mf")
SPLITS = ("training", "validation", "test")
CSV_NAMES = tuple(f"{aspect}-{split_name}.csv" for aspect in ASPECTS for split_name in SPLITS)
ASSIGNMENT_NAME = "protein_cluster_assignments.tsv"
CLUSTER_SPLIT_NAME = "cluster_split_assignments.tsv"
SCHEMA_NAME = "homology-threshold-progression-audit"
SCHEMA_VERSION = 1
LABEL_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(frozen=True)
class CsvRow:
    split: str
    sequence_sha256: str
    row_sha256: str


@dataclass
class Benchmark:
    label: str
    root: Path
    file_evidence: dict[str, dict[str, object]]
    proteins: dict[str, tuple[str, str]]
    csv_rows: dict[str, dict[str, CsvRow]]
    unassigned_proteins: int
    full_assignment_evidence: dict[str, object] | None


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _log(message: str) -> None:
    print(f"[{_utc_now()}] {message}", file=sys.stderr, flush=True)


def _sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_spec(raw: str) -> tuple[str, Path]:
    if "=" not in raw:
        raise argparse.ArgumentTypeError("benchmark must be LABEL=PATH")
    label, path = raw.split("=", 1)
    if not LABEL_RE.fullmatch(label):
        raise argparse.ArgumentTypeError(f"invalid benchmark label: {label!r}")
    return label, Path(path).expanduser().resolve()


def _read_assignments(path: Path) -> tuple[dict[str, tuple[str, str]], int]:
    proteins: dict[str, tuple[str, str]] = {}
    unassigned = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        required = {"uniprot_accession", "mmseqs_cluster_id", "split"}
        if not required.issubset(header):
            raise ValueError(f"{path} lacks required columns: {sorted(required - set(header))}")
        indices = {name: header.index(name) for name in required}
        for line_number, raw_line in enumerate(handle, start=2):
            columns = raw_line.rstrip("\r\n").split("\t")
            if len(columns) != len(header):
                raise ValueError(f"Malformed assignment row at {path}:{line_number}")
            protein = columns[indices["uniprot_accession"]]
            cluster = columns[indices["mmseqs_cluster_id"]]
            split_name = columns[indices["split"]]
            if not protein or protein in proteins:
                raise ValueError(f"Missing or duplicate protein at {path}:{line_number}: {protein!r}")
            if split_name and split_name not in SPLITS:
                raise ValueError(f"Unknown split at {path}:{line_number}: {split_name!r}")
            if bool(cluster) != bool(split_name):
                raise ValueError(
                    f"Cluster/split assignment disagree at {path}:{line_number}: "
                    f"cluster={cluster!r}, split={split_name!r}"
                )
            unassigned += int(not split_name)
            proteins[protein] = (cluster, split_name)
    return proteins, unassigned


def _read_csv(path: Path, expected_split: str) -> tuple[dict[str, CsvRow], dict[str, object]]:
    rows: dict[str, CsvRow] = {}
    digest = hashlib.sha256()
    size_bytes = 0
    with path.open("rb") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"Empty CSV: {path}")
        digest.update(header)
        size_bytes += len(header)
        first, separator, remainder = header.partition(b",")
        second, separator_two, _ = remainder.partition(b",")
        if first != b"proteins" or not separator or second != b"sequences" or not separator_two:
            raise ValueError(f"Unexpected first CSV columns in {path}")
        header_sha256 = _sha256_bytes(header.rstrip(b"\r\n"))
        for line_number, raw_line in enumerate(handle, start=2):
            digest.update(raw_line)
            size_bytes += len(raw_line)
            protein_raw, separator, remainder = raw_line.partition(b",")
            sequence_raw, separator_two, _ = remainder.partition(b",")
            if not separator or not separator_two:
                raise ValueError(f"Malformed CSV row at {path}:{line_number}")
            try:
                protein = protein_raw.decode("ascii")
                sequence = sequence_raw.decode("ascii")
            except UnicodeDecodeError as error:
                raise ValueError(f"Non-ASCII protein or sequence at {path}:{line_number}") from error
            if not protein or not sequence or protein in rows:
                raise ValueError(f"Missing or duplicate protein at {path}:{line_number}: {protein!r}")
            rows[protein] = CsvRow(
                split=expected_split,
                sequence_sha256=_sha256_bytes(sequence_raw),
                row_sha256=_sha256_bytes(raw_line.rstrip(b"\r\n")),
            )
    return rows, {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size_bytes": size_bytes,
        "data_rows": len(rows),
        "header_sha256": header_sha256,
    }


def _load_benchmark(
    label: str, root: Path, full_assignment: Path | None = None
) -> Benchmark:
    if not root.is_dir():
        raise FileNotFoundError(root)
    required = (*CSV_NAMES, ASSIGNMENT_NAME, CLUSTER_SPLIT_NAME)
    missing = [name for name in required if not (root / name).is_file()]
    if missing:
        raise FileNotFoundError(f"{root} is missing: {', '.join(missing)}")

    _log(f"{label}: reading protein-to-cluster assignments")
    proteins, unassigned = _read_assignments(root / ASSIGNMENT_NAME)
    file_evidence: dict[str, dict[str, object]] = {}
    for name in (ASSIGNMENT_NAME, CLUSTER_SPLIT_NAME):
        path = root / name
        file_evidence[name] = {
            "path": str(path),
            "sha256": _sha256(path),
            "size_bytes": path.stat().st_size,
        }

    csv_rows: dict[str, dict[str, CsvRow]] = {aspect: {} for aspect in ASPECTS}
    for name in CSV_NAMES:
        aspect, split_with_suffix = name.split("-", 1)
        split_name = split_with_suffix.removesuffix(".csv")
        _log(f"{label}: hashing and validating {name}")
        rows, evidence = _read_csv(root / name, split_name)
        overlap = set(rows).intersection(csv_rows[aspect])
        if overlap:
            raise ValueError(
                f"{label}/{aspect} has {len(overlap)} proteins in multiple splits; "
                f"sample={sorted(overlap)[:5]}"
            )
        for protein in rows:
            assignment = proteins.get(protein)
            if assignment is None:
                raise ValueError(f"{label}/{name} protein missing from assignment table: {protein}")
            if assignment[1] != split_name:
                raise ValueError(
                    f"{label}/{name} split disagrees with assignment table for {protein}: "
                    f"{assignment[1]!r}"
                )
        csv_rows[aspect].update(rows)
        file_evidence[name] = evidence

    assigned = sum(bool(split_name) for _, split_name in proteins.values())
    csv_union = set().union(*(set(rows) for rows in csv_rows.values()))
    if len(csv_union) != assigned:
        missing_from_csv = sorted(
            protein for protein, (_, split_name) in proteins.items()
            if split_name and protein not in csv_union
        )
        raise ValueError(
            f"{label} assigned-protein population differs from nine-CSV union: "
            f"assigned={assigned}, csv_union={len(csv_union)}, sample={missing_from_csv[:5]}"
        )
    full_assignment_evidence = None
    if full_assignment is not None:
        if not full_assignment.is_file():
            raise FileNotFoundError(full_assignment)
        _log(f"{label}: hashing complete UniRef50 assignment file")
        full_assignment_evidence = {
            "path": str(full_assignment),
            "sha256": _sha256(full_assignment),
            "size_bytes": full_assignment.stat().st_size,
        }
    return Benchmark(
        label, root, file_evidence, proteins, csv_rows, unassigned,
        full_assignment_evidence,
    )


def _ratio(numerator: int | float, denominator: int | float) -> float | None:
    return numerator / denominator if denominator else None


def _partition_metrics(
    left: dict[str, tuple[str, str]], right: dict[str, tuple[str, str]]
) -> dict[str, object]:
    if set(left) != set(right):
        raise ValueError("Protein assignment universes differ")
    pair_counts: Counter[tuple[str, str]] = Counter()
    left_counts: Counter[str] = Counter()
    right_counts: Counter[str] = Counter()
    left_blocks: defaultdict[str, list[str]] = defaultdict(list)
    right_blocks: defaultdict[str, list[str]] = defaultdict(list)
    for protein in left:
        left_cluster = left[protein][0] or f"__unassigned_left__:{protein}"
        right_cluster = right[protein][0] or f"__unassigned_right__:{protein}"
        pair_counts[left_cluster, right_cluster] += 1
        left_counts[left_cluster] += 1
        right_counts[right_cluster] += 1
        left_blocks[left_cluster].append(protein)
        right_blocks[right_cluster].append(protein)

    together_both = sum(count * (count - 1) // 2 for count in pair_counts.values())
    together_left = sum(count * (count - 1) // 2 for count in left_counts.values())
    together_right = sum(count * (count - 1) // 2 for count in right_counts.values())
    n = len(left)
    total_pairs = n * (n - 1) // 2
    expected = together_left * together_right / total_pairs if total_pairs else 0.0
    maximum = (together_left + together_right) / 2
    exact_left = Counter(hashlib.sha256("\0".join(sorted(members)).encode()).hexdigest() for members in left_blocks.values())
    exact_right = Counter(hashlib.sha256("\0".join(sorted(members)).encode()).hexdigest() for members in right_blocks.values())
    exact_blocks = exact_left & exact_right
    left_sizes = {
        hashlib.sha256("\0".join(sorted(members)).encode()).hexdigest(): len(members)
        for members in left_blocks.values()
    }
    exact_members = sum(left_sizes[key] * count for key, count in exact_blocks.items())
    return {
        "proteins": n,
        "left_clusters": len(left_counts),
        "right_clusters": len(right_counts),
        "exact_cluster_blocks": sum(exact_blocks.values()),
        "proteins_in_exact_cluster_blocks": exact_members,
        "partitions_exactly_identical": exact_members == n and len(left_counts) == len(right_counts),
        "pairs_together_in_left": together_left,
        "pairs_together_in_right": together_right,
        "pairs_together_in_both": together_both,
        "left_pair_recall": _ratio(together_both, together_left),
        "right_pair_recall": _ratio(together_both, together_right),
        "pair_jaccard": _ratio(together_both, together_left + together_right - together_both),
        "fowlkes_mallows": _ratio(together_both, math.sqrt(together_left * together_right)),
        "adjusted_rand_index": _ratio(together_both - expected, maximum - expected),
    }


def _transition_counts(
    left: dict[str, str], right: dict[str, str]
) -> tuple[dict[str, int], dict[str, object]]:
    matrix: Counter[str] = Counter()
    for protein in set(left).union(right):
        matrix[f"{left.get(protein, 'absent')}->{right.get(protein, 'absent')}"] += 1
    same = sum(count for transition, count in matrix.items() if transition.split("->")[0] == transition.split("->")[1])
    return dict(sorted(matrix.items())), {
        "union_proteins": len(set(left).union(right)),
        "common_proteins": len(set(left).intersection(right)),
        "left_only_proteins": len(set(left) - set(right)),
        "right_only_proteins": len(set(right) - set(left)),
        "same_state_proteins": same,
        "changed_state_proteins": len(set(left).union(right)) - same,
    }


def _compare(left: Benchmark, right: Benchmark) -> dict[str, object]:
    if set(left.proteins) != set(right.proteins):
        raise ValueError(f"Protein universes differ between {left.label} and {right.label}")
    assigned_left = {
        protein: assignment for protein, assignment in left.proteins.items() if assignment[1]
    }
    assigned_right = {
        protein: assignment for protein, assignment in right.proteins.items() if assignment[1]
    }
    global_left = {protein: split_name for protein, (_, split_name) in assigned_left.items()}
    global_right = {protein: split_name for protein, (_, split_name) in assigned_right.items()}
    global_matrix, global_summary = _transition_counts(global_left, global_right)
    aspects: dict[str, object] = {}
    for aspect in ASPECTS:
        left_rows = left.csv_rows[aspect]
        right_rows = right.csv_rows[aspect]
        matrix, summary = _transition_counts(
            {protein: row.split for protein, row in left_rows.items()},
            {protein: row.split for protein, row in right_rows.items()},
        )
        common = set(left_rows).intersection(right_rows)
        sequence_disagreements = sum(
            left_rows[protein].sequence_sha256 != right_rows[protein].sequence_sha256
            for protein in common
        )
        row_disagreements = sum(
            left_rows[protein].row_sha256 != right_rows[protein].row_sha256
            for protein in common
        )
        if sequence_disagreements:
            raise ValueError(
                f"Sequence content differs for {sequence_disagreements} {aspect} proteins "
                f"between {left.label} and {right.label}"
            )
        summary.update({
            "transition_matrix": matrix,
            "common_sequence_disagreements": sequence_disagreements,
            "common_complete_row_disagreements": row_disagreements,
        })
        aspects[aspect] = summary

    full_assignment_byte_identical = None
    if left.full_assignment_evidence is not None and right.full_assignment_evidence is not None:
        full_assignment_byte_identical = (
            left.full_assignment_evidence["sha256"]
            == right.full_assignment_evidence["sha256"]
        )
    core_files_identical = all(
        left.file_evidence[name]["sha256"] == right.file_evidence[name]["sha256"]
        for name in (*CSV_NAMES, ASSIGNMENT_NAME, CLUSTER_SPLIT_NAME)
    )
    return {
        "left": left.label,
        "right": right.label,
        "all_required_files_byte_identical": (
            core_files_identical and full_assignment_byte_identical is not False
        ),
        "full_assignment_byte_identical": full_assignment_byte_identical,
        "byte_identical_files": [
            name for name in (*CSV_NAMES, ASSIGNMENT_NAME, CLUSTER_SPLIT_NAME)
            if left.file_evidence[name]["sha256"] == right.file_evidence[name]["sha256"]
        ],
        "global_split": {**global_summary, "transition_matrix": global_matrix},
        "retained_partition": _partition_metrics(assigned_left, assigned_right),
        "aspects": aspects,
    }


def _benchmark_summary(benchmark: Benchmark) -> dict[str, object]:
    split_counts = Counter(split_name or "unassigned" for _, split_name in benchmark.proteins.values())
    return {
        "label": benchmark.label,
        "path": str(benchmark.root),
        "protein_universe": len(benchmark.proteins),
        "assigned_proteins": len(benchmark.proteins) - benchmark.unassigned_proteins,
        "unassigned_proteins": benchmark.unassigned_proteins,
        "split_counts": dict(sorted(split_counts.items())),
        "aspect_populations": {aspect: len(benchmark.csv_rows[aspect]) for aspect in ASPECTS},
        "files": benchmark.file_evidence,
        "full_assignment": benchmark.full_assignment_evidence,
    }


def _pct(value: int, total: int) -> str:
    return f"{100 * value / total:.3f}%" if total else "n/a"


def _metric(value: object) -> str:
    return f"{value:.6f}" if isinstance(value, float) else str(value)


def _write_outputs(output: Path, payload: dict[str, object]) -> None:
    output.mkdir(parents=True, exist_ok=False)
    (output / "audit.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    with (output / "file_hashes.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("threshold\tfile\tsha256\tsize_bytes\tdata_rows\n")
        for benchmark in payload["benchmarks"]:
            for name, evidence in benchmark["files"].items():
                handle.write(
                    f"{benchmark['label']}\t{name}\t{evidence['sha256']}\t"
                    f"{evidence['size_bytes']}\t{evidence.get('data_rows', '')}\n"
                )
            evidence = benchmark["full_assignment"]
            if evidence is not None:
                handle.write(
                    f"{benchmark['label']}\tfull_cluster_assignments.tsv.gz\t"
                    f"{evidence['sha256']}\t{evidence['size_bytes']}\t\n"
                )

    with (output / "split_transitions.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("left\tright\taspect\ttransition\tproteins\n")
        for comparison in payload["adjacent_comparisons"]:
            for transition, count in comparison["global_split"]["transition_matrix"].items():
                handle.write(
                    f"{comparison['left']}\t{comparison['right']}\tall\t{transition}\t{count}\n"
                )
            for aspect, summary in comparison["aspects"].items():
                for transition, count in summary["transition_matrix"].items():
                    handle.write(
                        f"{comparison['left']}\t{comparison['right']}\t{aspect}\t"
                        f"{transition}\t{count}\n"
                    )

    lines = [
        "# Homology threshold progression audit",
        "",
        "## Scope",
        "",
        "This report compares the published framework-stream benchmarks in the supplied order. It hashes every model CSV and retained-assignment artifact, checks the exact protein/sequence contracts, and measures adjacent split and retained-cluster changes. The full UniRef50 partition comparison is a separate, more expensive audit.",
        "",
        "## Threshold inputs",
        "",
        "| Threshold | Protein universe | Assigned | Unassigned | BP | CC | MF |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for benchmark in payload["benchmarks"]:
        lines.append(
            f"| {benchmark['label']} | {benchmark['protein_universe']:,} | "
            f"{benchmark['assigned_proteins']:,} | {benchmark['unassigned_proteins']:,} | "
            f"{benchmark['aspect_populations']['bp']:,} | "
            f"{benchmark['aspect_populations']['cc']:,} | "
            f"{benchmark['aspect_populations']['mf']:,} |"
        )
    lines.extend([
        "",
        "## Adjacent comparisons",
        "",
        "| Pair | All artifacts identical | Global split moves | Retained partition ARI | Exact-cluster protein coverage | BP moves | CC moves | MF moves |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for comparison in payload["adjacent_comparisons"]:
        partition = comparison["retained_partition"]
        proteins = int(partition["proteins"])
        exact_members = int(partition["proteins_in_exact_cluster_blocks"])
        lines.append(
            f"| {comparison['left']} to {comparison['right']} | "
            f"{'yes' if comparison['all_required_files_byte_identical'] else 'no'} | "
            f"{comparison['global_split']['changed_state_proteins']:,} | "
            f"{_metric(partition['adjusted_rand_index'])} | "
            f"{exact_members:,}/{proteins:,} ({_pct(exact_members, proteins)}) | "
            f"{comparison['aspects']['bp']['changed_state_proteins']:,} | "
            f"{comparison['aspects']['cc']['changed_state_proteins']:,} | "
            f"{comparison['aspects']['mf']['changed_state_proteins']:,} |"
        )
    lines.extend([
        "",
        "## Interpretation boundary",
        "",
        "- A high ARI means the retained annotated proteins are grouped similarly; it does not prove every pair of sequences satisfies the nominal identity threshold.",
        "- Split moves show how much the actual model inputs change between thresholds.",
        "- Byte-identical 10% and 5% artifacts would prove those two published benchmarks are operationally identical, not that MMseqs found every remote-homology edge.",
        "- Root-excluded performance diagnostics are needed alongside this construction audit before interpreting flat model scores.",
    ])
    (output / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    manifest = []
    for path in sorted(output.iterdir()):
        if path.is_file():
            manifest.append({
                "path": path.name,
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            })
    (output / "output_manifest.json").write_text(
        json.dumps({"files": manifest}, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    complete = {
        "complete": True,
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "audit_sha256": _sha256(output / "audit.json"),
        "summary_sha256": _sha256(output / "summary.md"),
    }
    (output / "RUN_COMPLETE.json").write_text(
        json.dumps(complete, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def run(args: argparse.Namespace) -> dict[str, object]:
    specs = [_parse_spec(raw) for raw in args.benchmark]
    if len(specs) < 2:
        raise ValueError("At least two --benchmark LABEL=PATH values are required")
    labels = [label for label, _ in specs]
    if len(set(labels)) != len(labels):
        raise ValueError("Benchmark labels must be unique")
    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise FileExistsError(output)

    summaries: list[dict[str, object]] = []
    comparisons: list[dict[str, object]] = []
    previous: Benchmark | None = None
    full_specs = {
        full_label: full_path
        for full_label, full_path in (
            _parse_spec(raw) for raw in getattr(args, "full_assignment", [])
        )
    }
    unknown_full_labels = set(full_specs) - set(labels)
    if unknown_full_labels:
        raise ValueError(
            f"Full-assignment labels lack benchmarks: {sorted(unknown_full_labels)}"
        )
    if full_specs and set(full_specs) != set(labels):
        raise ValueError("Supply --full-assignment for every benchmark label or none")

    for label, root in specs:
        current = _load_benchmark(label, root, full_specs.get(label))
        summaries.append(_benchmark_summary(current))
        if previous is not None:
            _log(f"comparing {previous.label} to {current.label}")
            comparisons.append(_compare(previous, current))
        previous = current
    payload: dict[str, object] = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": _utc_now(),
        "benchmarks": summaries,
        "adjacent_comparisons": comparisons,
    }
    _write_outputs(output, payload)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark",
        action="append",
        required=True,
        metavar="LABEL=PATH",
        help="Published benchmark directory; repeat in comparison order.",
    )
    parser.add_argument(
        "--full-assignment",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help="Complete UniRef50 representative/member assignment; repeat for every label.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = run(args)
    print(json.dumps({
        "benchmarks": len(payload["benchmarks"]),
        "comparisons": len(payload["adjacent_comparisons"]),
        "output_dir": str(args.output_dir.expanduser().resolve()),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
