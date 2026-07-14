from __future__ import annotations

import csv
from itertools import zip_longest
from pathlib import Path
import statistics
from typing import Iterable

import pandas as pd

from .config import PREFIX_TO_NAMESPACE, SPLITS, BuildConfig
from .inputs import open_text
from .labels import _protein_clusters
from .models import (
    ClusterInfo,
    LabelBuildResult,
    MappingDecision,
    SplitAssignment,
    ValidationReport,
)
from .mmseqs import ClusterIndex, CommandSpec
from .ontology import Ontology
from .splitting import assign_splits
from .uniref import SEQUENCE_RE, UniRefIndex


class _BoundedErrors(list[str]):
    def __init__(self, limit: int):
        super().__init__()
        self.limit = limit
        self.total = 0

    def append(self, item: str) -> None:
        self.total += 1
        if len(self) < self.limit:
            super().append(item)


def _cluster_manifest_errors(
    output_dir: Path,
    assignments: dict[str, SplitAssignment],
    cluster_index: ClusterIndex,
    uniref: UniRefIndex,
) -> list[str]:
    errors: list[str] = []
    with (output_dir / "cluster_split_assignments.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_assignments = (
            assignments[cluster_id] for cluster_id in sorted(assignments)
        )
        missing = object()
        for row_number, pair in enumerate(
            zip_longest(reader, expected_assignments, fillvalue=missing), start=2
        ):
            row, expected = pair
            if row is missing or expected is missing:
                errors.append(f"cluster-assignment-row-count:{row_number}")
                break
            assert isinstance(row, dict)
            expected_values = {
                "mmseqs_cluster_id": expected.cluster_id,
                "split": expected.split,
                "uniref90_member_count": str(expected.member_count),
                "qualifying_uniprot_count": str(expected.labelled_protein_count),
                "assignment_stage": expected.stage,
            }
            if any(row.get(key) != value for key, value in expected_values.items()):
                errors.append(f"cluster-assignment-content:{row_number}:{expected.cluster_id}")
                if len(errors) >= 50:
                    break

    with open_text(output_dir / "retained_cluster_members.tsv.gz") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_rows = (
            (cluster_id, member_id, digest, length)
            for cluster_id, member_id, digest, length
            in cluster_index.iter_assignments_with_metadata_for_clusters(
                uniref, assignments
            )
        )
        missing = object()
        for row_number, pair in enumerate(
            zip_longest(reader, expected_rows, fillvalue=missing), start=2
        ):
            row, expected_row = pair
            if row is missing or expected_row is missing:
                errors.append(f"retained-member-row-count:{row_number}")
                break
            assert isinstance(row, dict)
            cluster_id, member_id, digest, length = expected_row
            expected_assignment = assignments[cluster_id]
            expected_values = {
                "mmseqs_cluster_id": cluster_id,
                "split": expected_assignment.split,
                "uniref90_id": member_id,
                "sequence_sha256": digest,
                "sequence_length": str(length),
            }
            if any(row.get(key) != value for key, value in expected_values.items()):
                errors.append(f"retained-member-content:{row_number}:{member_id}")
                if len(errors) >= 50:
                    break
    with open_text(output_dir / "mmseqs_cluster_membership.tsv.gz") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_rows = cluster_index.iter_assignments()
        missing = object()
        for row_number, pair in enumerate(
            zip_longest(reader, expected_rows, fillvalue=missing), start=2
        ):
            row, expected_row = pair
            if row is missing or expected_row is missing:
                errors.append(f"canonical-membership-row-count:{row_number}")
                break
            assert isinstance(row, dict)
            cluster_id, member_id = expected_row
            if row != {"mmseqs_cluster_id": cluster_id, "uniref90_id": member_id}:
                errors.append(f"canonical-membership-content:{row_number}:{member_id}")
                if len(errors) >= 50:
                    break
    return errors


def validate_term_universe_artifacts(
    output_dir: Path, min_count: int
) -> tuple[list[str], tuple[str, ...]]:
    """Recount terms only from serialized development artifacts, independently of build state."""
    errors: list[str] = []
    support: dict[str, set[str]] = {}
    for name in ("train_data_train.pkl", "train_data_valid.pkl"):
        frame = pd.read_pickle(output_dir / name)
        if list(frame.columns) != ["proteins", "sequences", "annotations"]:
            errors.append(f"schema:{name}")
            continue
        for row in frame.itertuples(index=False):
            protein = str(row.proteins)
            for term in set(row.annotations):
                support.setdefault(str(term), set()).add(protein)
    expected = tuple(sorted(term for term, proteins in support.items() if len(proteins) >= min_count))
    terms = pd.read_pickle(output_dir / "terms.pkl")
    if list(terms.columns) != ["terms"]:
        errors.append("schema:terms.pkl")
        observed: list[str] = []
    else:
        observed = [str(term) for term in terms["terms"].tolist()]
    if observed != sorted(observed):
        errors.append("order:terms.pkl")
    if len(observed) != len(set(observed)):
        errors.append("duplicate:terms.pkl")
    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))
    if missing:
        errors.append("missing-eligible:" + ",".join(missing[:20]))
    if unexpected:
        errors.append("unexpected-or-test-only:" + ",".join(unexpected[:20]))
    return errors, expected


def split_balance_metrics(
    assignments: dict[str, SplitAssignment], config: BuildConfig
) -> dict[str, float]:
    member_counts = {
        split: sum(item.member_count for item in assignments.values() if item.split == split)
        for split in SPLITS
    }
    total = sum(member_counts.values())
    development = member_counts["training"] + member_counts["validation"]
    development_ratio = development / total if total else 0.0
    training_ratio = member_counts["training"] / development if development else 0.0
    return {
        "development_ratio": development_ratio,
        "training_ratio_within_development": training_ratio,
        "development_deviation": abs(development_ratio - config.development_fraction),
        "training_deviation": abs(
            training_ratio - config.training_fraction_within_development
        ),
    }


def production_balance_within_tolerance(
    assignments: dict[str, SplitAssignment], config: BuildConfig
) -> bool:
    metrics = split_balance_metrics(assignments, config)
    return (
        config.fixture_mode
        or config.split_policy != "sequence-balanced"
        or (
            metrics["development_deviation"] <= 0.05
            and metrics["training_deviation"] <= 0.05
        )
    )


def _check(report: ValidationReport, name: str, condition: bool, detail: str, **metrics: object) -> None:
    report.add_check(name, condition, detail, **metrics)


def _annotation_audit_errors(output_dir: Path, goa, config: BuildConfig) -> list[str]:
    errors: list[str] = []
    with open_text(output_dir / "qualifying_annotations.tsv.gz") as handle:
        accepted_rows = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
    if accepted_rows != goa.counters["kept_rows"]:
        errors.append(f"accepted-row-count:{accepted_rows}!={goa.counters['kept_rows']}")

    observed_counts: dict[tuple[str, ...], int] = {}
    with (output_dir / "annotation_decision_counts.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (
                row["disposition"], row["reason"], row["evidence_code"],
                row["database"], row["object_type"], row["aspect"],
            )
            observed_counts[key] = int(row["rows"])
    if observed_counts != dict(goa.annotation_decision_counts):
        errors.append("complete-decision-counts")

    sample_counts: dict[str, int] = {}
    with open_text(output_dir / "excluded_annotations_sample.tsv.gz") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            reason = row["rejection_reason"]
            sample_counts[reason] = sample_counts.get(reason, 0) + 1
    if sample_counts != dict(goa.excluded_sample_counts):
        errors.append("excluded-sample-counts")
    if any(count > config.excluded_sample_per_reason for count in sample_counts.values()):
        errors.append("excluded-sample-bound")
    for forbidden in (
        "uniref90_clusters.tsv", "uniref90_to_mmseqs_cluster.tsv",
        "retained_cluster_members.tsv", "qualifying_annotations.tsv", "excluded_annotations.tsv",
    ):
        if (output_dir / forbidden).exists():
            errors.append(f"forbidden-uncompressed-or-duplicate:{forbidden}")
    return errors


def validate_outputs(
    output_dir: Path,
    config: BuildConfig,
    ontology: Ontology,
    goa,
    labels: LabelBuildResult,
    mappings: list[MappingDecision],
    assignments: dict[str, SplitAssignment],
    commands: Iterable[CommandSpec],
    cluster_index: ClusterIndex,
    uniref: UniRefIndex,
    uniref_count: int,
    mmseqs_member_count: int,
    global_sequence_conflicts: list[dict[str, str]],
) -> ValidationReport:
    report = ValidationReport()
    _check(
        report, "mmseqs_assignment_completeness", uniref_count == mmseqs_member_count,
        "Every frozen UniRef90 entry must receive exactly one MMseqs2 assignment",
        uniref90_entries=uniref_count, assigned_members=mmseqs_member_count,
    )
    manifest_errors = _cluster_manifest_errors(
        output_dir, assignments, cluster_index, uniref
    )
    annotation_audit_errors = _annotation_audit_errors(output_dir, goa, config)
    _check(
        report, "bounded_compressed_annotation_audits", not annotation_audit_errors,
        "Accepted annotations and complete rejection decisions reconcile; detailed exclusions "
        "are deterministic bounded samples and large manifests use canonical gzip artifacts",
        errors=annotation_audit_errors,
        excluded_sample_per_reason=config.excluded_sample_per_reason,
    )
    uniref_sequence_conflicts = [
        item for item in global_sequence_conflicts
        if "," in item["uniref90_splits"]
    ]
    annotated_sequence_conflicts = [
        item for item in global_sequence_conflicts
        if "," in item["uniprot_splits"]
    ]
    _check(
        report, "retained_cluster_manifest_consistency", not manifest_errors,
        "Every retained cluster occurs once and every retained member agrees with its global split",
        clusters=len(assignments), errors=manifest_errors[:50],
    )
    observed_splits = {assignment.split for assignment in assignments.values()}
    _check(
        report, "all_three_splits_present", observed_splits == set(SPLITS),
        "Train, validation, and test must all receive retained clusters",
        observed=sorted(observed_splits),
    )

    retained_proxy = {
        key: ClusterInfo(key, value.member_count, value.labelled_protein_count)
        for key, value in assignments.items()
    }
    replay = assign_splits(
        retained_proxy,
        config.split_policy,
        config.seed,
        config.development_fraction,
        config.training_fraction_within_development,
    )
    _check(
        report, "deterministic_split_replay",
        {key: value.split for key, value in replay.items()}
        == {key: value.split for key, value in assignments.items()},
        "Replaying the split algorithm with the same seed gives identical assignments",
    )

    mapped_cluster = _protein_clusters(mappings)
    global_split = {
        protein: assignments[cluster].split
        for protein, cluster in mapped_cluster.items()
        if cluster in assignments
    }
    protein_seen: dict[str, str] = {}
    leakage: list[str] = []
    for split in SPLITS:
        frame = labels.frames[split]
        for row in frame.itertuples(index=False):
            protein = str(row.proteins)
            previous = protein_seen.setdefault(protein, split)
            if previous != split:
                leakage.append(protein)
            if global_split.get(protein) != split:
                leakage.append(protein)
    _check(
        report, "global_protein_disjointness", not leakage,
        "No UniProtKB accession may cross train/validation/test and every label row follows its cluster",
        conflicts=sorted(set(leakage))[:20],
    )
    _check(
        report, "global_exact_sequence_disjointness", not annotated_sequence_conflicts,
        "No exact available sequence for any mapped qualifying UniProtKB protein may cross splits, "
        "including proteins later lacking an evaluable PFP term",
        conflicts=annotated_sequence_conflicts,
    )
    _check(
        report, "retained_uniref90_exact_sequence_disjointness", not uniref_sequence_conflicts,
        "No exact frozen UniRef90 scaffold sequence may cross train/validation/test",
        conflicts=uniref_sequence_conflicts,
    )
    _check(
        report, "global_retained_exact_sequence_disjointness", not global_sequence_conflicts,
        "No exact sequence crosses splits across the combined retained UniRef90 scaffold and "
        "mapped qualifying UniProtKB population",
        conflicts=global_sequence_conflicts,
    )

    universe = set(labels.term_universe)
    outside = {split: [] for split in ("validation", "test")}
    restriction_errors: list[str] = []
    for split in outside:
        proteins = set(labels.frames[split]["proteins"])
        outside[split] = sorted({
            term
            for protein in proteins
            for term in labels.restricted_annotations.get(str(protein), ())
            if term not in universe
        })
        for row in labels.frames[split].itertuples(index=False):
            expected = tuple(sorted(set(row.annotations) & universe))
            if labels.restricted_annotations.get(str(row.proteins), ()) != expected:
                restriction_errors.append(str(row.proteins))
    _check(
        report, "evaluation_terms_development_defined",
        not outside["validation"] and not outside["test"] and not restriction_errors,
        "All split labels are restricted to the complete development-defined term universe",
        outside=outside, restriction_errors=restriction_errors[:20],
    )

    cluster_command = next((item for item in commands if item.stage == "cluster"), None)
    command_text = cluster_command.argv if cluster_command else ()
    required_fragments = {
        "--min-seq-id": f"{config.identity:.2f}",
        "-c": "0.8",
        "--cov-mode": "0",
        "--cluster-mode": "0",
        "--alignment-mode": "3",
        "--seq-id-mode": "0",
        "--cluster-reassign": "1",
        "-s": "7.5",
        "-e": "1e-4",
    }
    command_valid = cluster_command is not None
    for flag, expected in required_fragments.items():
        if flag not in command_text:
            command_valid = False
            continue
        index = command_text.index(flag)
        command_valid = command_valid and index + 1 < len(command_text) and command_text[index + 1] == expected
    _check(
        report, "locked_mmseqs_parameters", command_valid,
        "The recorded MMseqs2 command contains the locked identity, coverage, modes, sensitivity, "
        "and E-value policy",
        required=required_fragments,
    )
    createdb_command = next((item for item in commands if item.stage == "createdb"), None)
    createdb_text = createdb_command.argv if createdb_command else ()
    deterministic_createdb = (
        "--shuffle" in createdb_text
        and createdb_text[createdb_text.index("--shuffle") + 1] == "0"
    )
    _check(
        report, "deterministic_mmseqs_input_order", deterministic_createdb,
        "MMseqs2 createdb has input shuffling explicitly disabled",
    )

    csv_errors = _BoundedErrors(50)
    csv_rows: dict[str, int] = {}
    for prefix, namespace in PREFIX_TO_NAMESPACE.items():
        expected_terms = [
            term for term in labels.term_universe if ontology.namespace(term) == namespace
        ]
        exact_header = ["proteins", "sequences", *expected_terms]
        for split in SPLITS:
            path = output_dir / f"{prefix}-{split}.csv"
            if not path.is_file():
                csv_errors.append(f"missing:{path.name}")
                continue
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.reader(handle)
                header = next(reader, None)
                if header != exact_header:
                    csv_errors.append(f"exact-header:{path.name}")
                    continue
                if len(header) < 3 and not config.allow_empty_fixture_outputs:
                    csv_errors.append(f"no-go-columns:{path.name}")
                def expected_rows():
                    for label_row in labels.frames[split].sort_values(
                        "proteins", kind="stable"
                    ).itertuples(index=False):
                        annotation_set = set(label_row.annotations)
                        if not annotation_set.intersection(expected_terms):
                            continue
                        yield [
                            str(label_row.proteins),
                            str(label_row.sequences),
                            *("1" if term in annotation_set else "0" for term in expected_terms),
                        ]

                missing = object()
                rows = 0
                for row_number, pair in enumerate(
                    zip_longest(reader, expected_rows(), fillvalue=missing), start=2
                ):
                    row, expected_row = pair
                    if row is missing or expected_row is missing:
                        csv_errors.append(f"row-count:{path.name}")
                        if row is missing:
                            continue
                    assert isinstance(row, list)
                    if expected_row is not missing and row != expected_row:
                        csv_errors.append(f"exact-content:{path.name}:{row_number}")
                    if len(row) != len(header):
                        csv_errors.append(f"width:{path.name}:{row_number}")
                        continue
                    if not row[1] or SEQUENCE_RE.fullmatch(row[1]) is None:
                        csv_errors.append(f"sequence:{path.name}:{row_number}")
                    if any(value not in {"0", "1"} for value in row[2:]):
                        csv_errors.append(f"binary:{path.name}:{row_number}")
                    if not any(value == "1" for value in row[2:]):
                        csv_errors.append(f"all-zero:{path.name}:{row_number}")
                    if global_split.get(row[0]) != split:
                        csv_errors.append(f"global-split:{path.name}:{row[0]}")
                    rows += 1
                if rows == 0 and not config.allow_empty_fixture_outputs:
                    csv_errors.append(f"empty:{path.name}")
                csv_rows[path.name] = rows
    _check(
        report, "pfp_csv_contract", not csv_errors,
        "All nine CSVs have the immutable PFP schema, binary labels, stable training-defined GO columns, "
        "correct namespaces, nonempty valid sequences, and global split consistency",
        error_count=csv_errors.total, errors=list(csv_errors), rows=csv_rows,
    )

    pickle_errors: list[str] = []
    columns = ["proteins", "sequences", "annotations"]
    expected_pickles = {
        "train_data_train.pkl": labels.frames["training"].loc[:, columns].reset_index(drop=True),
        "train_data_valid.pkl": labels.frames["validation"].loc[:, columns].reset_index(drop=True),
        "test_data.pkl": labels.frames["test"].loc[:, columns].reset_index(drop=True),
    }
    expected_pickles["train_data.pkl"] = pd.concat(
        [expected_pickles["train_data_train.pkl"], expected_pickles["train_data_valid.pkl"]],
        ignore_index=True,
    ).sort_values("proteins", kind="stable").reset_index(drop=True)
    for name, expected_frame in expected_pickles.items():
        path = output_dir / name
        try:
            frame = pd.read_pickle(path)
            pd.testing.assert_frame_equal(frame, expected_frame, check_exact=True)
        except AssertionError as exc:
            pickle_errors.append(f"exact-content:{name}:{str(exc).splitlines()[0]}")
        except Exception as exc:  # pragma: no cover - exercised through failure paths
            pickle_errors.append(f"read:{name}:{exc}")
    try:
        terms = pd.read_pickle(output_dir / "terms.pkl")
        if list(terms.columns) != ["terms"] or terms["terms"].tolist() != list(labels.term_universe):
            pickle_errors.append("schema-or-order:terms.pkl")
    except Exception as exc:  # pragma: no cover
        pickle_errors.append(f"read:terms.pkl:{exc}")
    _check(
        report, "deepgoplus_pickle_contract", not pickle_errors,
        "The five DeepGOPlus-shaped pickle files exactly match the expected split rows, "
        "sequences, annotations, and training-defined terms",
        errors=pickle_errors,
    )
    independent_term_errors, independent_terms = validate_term_universe_artifacts(
        output_dir, config.min_count
    )
    _check(
        report, "independent_development_term_recount", not independent_term_errors,
        "terms.pkl is independently recounted from serialized training plus validation proteins; "
        "test-only and below-threshold terms are excluded",
        errors=independent_term_errors,
        recounted_term_count=len(independent_terms),
    )
    _check(
        report, "attrition_conservation",
        sum(labels.row_attrition_counts.values()) == labels.intended_annotation_rows
        and sum(labels.protein_attrition_counts.values()) == labels.intended_accessions,
        "Row-level and protein-candidate attrition buckets are mutually exclusive and reconcile "
        "to their declared denominators",
        intended_annotation_rows=labels.intended_annotation_rows,
        row_buckets=dict(labels.row_attrition_counts),
        intended_accessions=labels.intended_accessions,
        protein_buckets=dict(labels.protein_attrition_counts),
    )

    member_counts = {split: 0 for split in SPLITS}
    cluster_counts = {split: 0 for split in SPLITS}
    labelled_counts = {split: 0 for split in SPLITS}
    for assignment in assignments.values():
        member_counts[assignment.split] += assignment.member_count
        cluster_counts[assignment.split] += 1
        labelled_counts[assignment.split] += assignment.labelled_protein_count
    total_members = sum(member_counts.values()) or 1
    total_clusters = sum(cluster_counts.values()) or 1
    balance = split_balance_metrics(assignments, config)
    development_ratio = balance["development_ratio"]
    conditional_training_ratio = balance["training_ratio_within_development"]
    development_deviation = balance["development_deviation"]
    training_deviation = balance["training_deviation"]
    if config.split_policy == "sequence-balanced" and (
        development_deviation > 0.02 or training_deviation > 0.02
    ):
        report.add_warning(
            "sequence_balance",
            "Whole-cluster indivisibility leaves a sequence balance more than two percentage "
            "points from a requested stage target",
            requested_development=config.development_fraction,
            achieved_development=development_ratio,
            development_percentage_point_deviation=development_deviation * 100,
            requested_training_within_development=config.training_fraction_within_development,
            achieved_training_within_development=conditional_training_ratio,
            training_percentage_point_deviation=training_deviation * 100,
        )
    balance_within_production_tolerance = production_balance_within_tolerance(
        assignments, config
    )
    _check(
        report, "production_sequence_balance_tolerance",
        balance_within_production_tolerance,
        "Production sequence-balanced stages must each remain within five percentage points; "
        "fixture runs report but do not fail impossible indivisible cases",
        development_percentage_point_deviation=development_deviation * 100,
        training_percentage_point_deviation=training_deviation * 100,
        enforced=not config.fixture_mode and config.split_policy == "sequence-balanced",
    )
    giant = [
        assignment.cluster_id for assignment in assignments.values()
        if assignment.member_count >= config.giant_cluster_threshold
    ]
    if giant:
        report.add_warning(
            "giant_clusters", "Retained giant clusters can dominate balance and runtime",
            threshold=config.giant_cluster_threshold, clusters=giant[:20], count=len(giant),
        )
    sizes = [assignment.member_count for assignment in assignments.values()]
    if sizes:
        median_size = float(statistics.median(sizes))
        largest_share = max(sizes) / sum(sizes)
        max_to_median = max(sizes) / median_size if median_size else float("inf")
        if largest_share >= 0.50 or max_to_median >= 10:
            report.add_warning(
                "cluster_size_imbalance",
                "The retained cluster-size distribution is highly imbalanced",
                maximum=max(sizes), median=median_size,
                max_to_median=max_to_median, largest_member_share=largest_share,
            )
    report.add_check(
        "split_ratios_reported", True,
        "Requested and achieved cluster, UniRef90-member, and labelled-UniProt ratios were calculated",
        requested_development=config.development_fraction,
        member_counts=member_counts,
        cluster_counts=cluster_counts,
        labelled_protein_counts=labelled_counts,
        member_ratios={key: value / total_members for key, value in member_counts.items()},
        cluster_ratios={key: value / total_clusters for key, value in cluster_counts.items()},
        requested_training_within_development=config.training_fraction_within_development,
        achieved_development_member_ratio=development_ratio,
        achieved_training_member_ratio_within_development=conditional_training_ratio,
        largest_indivisible_cluster=(
            max(sizes) if sizes else 0
        ),
        largest_indivisible_cluster_share=(
            max(sizes) / sum(sizes) if sizes else 0
        ),
    )
    return report


def write_validation_report(report: ValidationReport, output_dir: Path) -> tuple[Path, Path]:
    import json

    json_path = output_dir / "validation_report.json"
    markdown_path = output_dir / "validation_report.md"
    payload = {"valid": report.valid, "checks": report.checks, "warnings": report.warnings}
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# Validation report", "", f"- Overall contract status: **{'PASS' if report.valid else 'FAIL'}**",
        f"- Checks: {len(report.checks)}", f"- Warnings: {len(report.warnings)}", "",
        "## Checks", "",
    ]
    for item in report.checks:
        lines.append(f"- [{'x' if item['passed'] else ' '}] `{item['name']}` — {item['detail']}")
    lines.extend(["", "## Warnings", ""])
    if report.warnings:
        for item in report.warnings:
            lines.append(f"- `{item['name']}` — {item['detail']}")
    else:
        lines.append("- None.")
    lines.extend([
        "", "## Boundary", "",
        "These checks establish software, schema, provenance, split, and leakage properties. "
        "They do not establish biological optimality or low-identity homology validity.", "",
    ])
    markdown_path.write_text("\n".join(lines), encoding="utf-8")
    return json_path, markdown_path
