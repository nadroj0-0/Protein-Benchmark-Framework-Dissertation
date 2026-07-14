from __future__ import annotations

import argparse
import csv
import json
import logging
import os
from pathlib import Path

import pandas as pd

from .attrition import load_attrition_policy
from .authorization import validate_pilot_approval
from .config import (
    SUPPORTED_IDENTITIES,
    UNIPROT_SOURCE_SCOPES,
    BuildConfig,
    parse_identity,
)
from .models import BuildResult, InputSpec
from .inputs import sha256_file
from .mmseqs import build_mmseqs_commands
from .pipeline import build_benchmark, validate_publication
from .provenance import (
    git_state,
    publish,
    staging_output,
    verify_output_manifest,
    write_completion_marker,
    write_output_manifest,
)


LOGGER = logging.getLogger(__name__)


def _add_input(parser: argparse.ArgumentParser, option: str, label: str) -> None:
    parser.add_argument(f"--{option}", type=Path, help=f"Explicit local {label} path")
    parser.add_argument(f"--{option}-url", help=f"Explicit frozen source URL for {label}")
    parser.add_argument(f"--{option}-sha256", help=f"Expected SHA-256 for {label}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="homology-cluster-benchmark",
        description="Build Daniel's frozen UniRef90/MMseqs2 whole-cluster PFP benchmark",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("build", help="Build one or all locked identity thresholds")
    build.add_argument(
        "--identity", action="append", default=[], metavar="{30,25,20,15,10,5,all}",
        help="Locked identity percentage; repeat or use 'all'",
    )
    build.add_argument("--split-policy", choices=("cluster-count-random", "sequence-balanced"), default="sequence-balanced")
    build.add_argument(
        "--training-population", choices=("annotated-only", "all-cluster-members"),
        default="annotated-only",
        help="annotated-only is implemented; all-cluster-members is explicitly unsupported",
    )
    _add_input(build, "uniref90-fasta", "UniRef90 FASTA")
    _add_input(build, "idmapping", "idmapping_selected.tab")
    build.add_argument(
        "--uniprot-source-scope",
        choices=UNIPROT_SOURCE_SCOPES,
        help="Required production supervised-population scope; UniRef90 remains the clustering scaffold",
    )
    _add_input(build, "uniprot-sprot-sequences", "frozen Swiss-Prot DAT")
    _add_input(build, "uniprot-trembl-sequences", "frozen TrEMBL DAT")
    _add_input(build, "goa", "frozen GOA GAF")
    _add_input(build, "go-obo", "frozen GO OBO")
    build.add_argument("--cluster-assignments", type=Path, help="Precomputed MMseqs2 createtsv output; intended for validated fixture tests")
    build.add_argument(
        "--fixture-mode", action="store_true",
        help=(
            "Permit synthetic/smoke inputs, precomputed assignments, and fixture thresholds; "
            "output is marked non-production"
        ),
    )
    build.add_argument("--mmseqs-bin", default=os.environ.get("MMSEQS_BIN", "mmseqs"))
    build.add_argument("--expected-mmseqs-version")
    build.add_argument("--frozen-input-manifest", type=Path)
    build.add_argument("--attrition-policy", type=Path)
    build.add_argument("--attrition-override", type=Path)
    build.add_argument("--framework-revision")
    build.add_argument(
        "--diagnostic-pilot", action="store_true",
        help="Publish measured 30% pilot evidence as non-production without self-approval",
    )
    build.add_argument("--output-dir", type=Path, required=True, help="Root beneath which identity/policy/population directories are published")
    build.add_argument(
        "--temp-dir", type=Path,
        default=Path(os.environ.get("TMPDIR", "/tmp")) / "homology-cluster-benchmark",
    )
    build.add_argument("--threads", type=int, default=1)
    build.add_argument("--requested-slots", type=int)
    build.add_argument("--allocated-slots", type=int)
    build.add_argument("--run-id", default=os.environ.get("RUN_ID", "local"))
    build.add_argument("--seed", type=int, default=0)
    build.add_argument("--min-count", type=int, default=50)
    build.add_argument(
        "--development-fraction", type=float, default=0.80,
        help="Locked at 0.80; exposed only so incompatible invocations fail explicitly",
    )
    build.add_argument(
        "--training-fraction-within-development", type=float, default=0.90,
        help="Locked at 0.90; exposed only so incompatible invocations fail explicitly",
    )
    build.add_argument(
        "--sensitivity", type=float, default=7.5,
        help="Locked at 7.5 by the reviewed MMseqs2 contract",
    )
    build.add_argument(
        "--evalue", type=float, default=1e-4,
        help="Locked at 1e-4 by the reviewed MMseqs2 contract",
    )
    build.add_argument("--scratch-safety-multiplier", type=float, default=8.0)
    build.add_argument("--minimum-free-disk-gb", type=float, default=0.0)
    build.add_argument("--persistent-results-root", type=Path)
    build.add_argument("--mmseqs-work-multiplier", type=float, default=8.0)
    build.add_argument("--publication-safety-multiplier", type=float, default=2.0)
    build.add_argument("--excluded-sample-per-reason", type=int, default=1000)
    build.add_argument("--uniprot-release", default="2026_02")
    build.add_argument("--goa-release", default="234")
    build.add_argument("--ontology-release", default="releases/2026-06-15")
    build.add_argument("--no-downloads", action="store_true")
    build.add_argument("--no-strict-qc", action="store_true")
    build.add_argument("--no-relationships", action="store_true")
    build.add_argument("--keep-temp", action="store_true")
    build.add_argument("--dry-run", action="store_true", help="Print command previews without resolving inputs or writing outputs")
    build.add_argument("--allow-empty-fixture-outputs", action="store_true", help=argparse.SUPPRESS)
    build.add_argument("--verbosity", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")

    validate = subparsers.add_parser("validate", help="Verify a published run manifest, marker, and saved validation status")
    validate.add_argument("--run-dir", type=Path, required=True)
    summarize = subparsers.add_parser(
        "summarize", help="Validate and aggregate six independently published threshold runs"
    )
    summarize.add_argument("--output-dir", type=Path, required=True)
    summarize.add_argument(
        "--run-dir", type=Path, action="append", required=True,
        help="Published run leaf or parent containing exactly one publication; repeat exactly six times",
    )
    authorize = subparsers.add_parser(
        "authorize-array",
        help="Validate reviewed attrition policy and 30% pilot approval before qsub",
    )
    authorize.add_argument("--attrition-policy", type=Path, required=True)
    authorize.add_argument("--pilot-approval", type=Path, required=True)
    authorize.add_argument("--pilot-completion-marker", type=Path, required=True)
    authorize.add_argument("--pilot-attrition-report", type=Path, required=True)
    authorize.add_argument("--pilot-run-dir", type=Path, required=True)
    authorize.add_argument("--pilot-task-context", type=Path, required=True)
    authorize.add_argument("--pilot-measurement-evidence", type=Path, required=True)
    authorize.add_argument("--frozen-input-manifest", type=Path, required=True)
    authorize.add_argument("--framework-revision", required=True)
    authorize.add_argument("--uniprot-source-scope", choices=UNIPROT_SOURCE_SCOPES, required=True)
    authorize.add_argument("--split-policy", choices=("cluster-count-random", "sequence-balanced"), required=True)
    authorize.add_argument("--training-population", required=True)
    authorize.add_argument("--expected-mmseqs-version", required=True)
    authorize.add_argument("--uniprot-release", default="2026_02")
    authorize.add_argument("--goa-release", default="234")
    authorize.add_argument("--ontology-release", default="releases/2026-06-15")
    return parser


def _input_spec(
    args: argparse.Namespace, name: str, release: str, source_population: str = "shared"
) -> InputSpec:
    attribute = name.replace("-", "_")
    return InputSpec(
        name=attribute,
        path=getattr(args, attribute),
        url=getattr(args, attribute + "_url"),
        expected_sha256=getattr(args, attribute + "_sha256"),
        release=release,
        source_population=source_population,
    )


def _identities(values: list[str]) -> tuple[float, ...]:
    values = values or ["30"]
    if any(value.strip().lower() == "all" for value in values):
        if len(values) != 1:
            raise ValueError("Use --identity all by itself")
        return SUPPORTED_IDENTITIES
    parsed = tuple(parse_identity(value) for value in values)
    if len(set(parsed)) != len(parsed):
        raise ValueError("Duplicate identity thresholds were requested")
    return parsed


def _config(args: argparse.Namespace, identity: float) -> BuildConfig:
    source_scope = args.uniprot_source_scope or (
        "sprot-only" if args.fixture_mode else ""
    )
    return BuildConfig(
        identity=identity,
        output_dir=args.output_dir,
        temp_dir=args.temp_dir,
        uniref90_fasta=_input_spec(
            args, "uniref90-fasta", args.uniprot_release,
            "uniref90-clustering-scaffold",
        ),
        idmapping=_input_spec(
            args, "idmapping", args.uniprot_release, "uniprotkb-shared-mapping"
        ),
        uniprot_source_scope=source_scope,
        uniprot_sprot_sequences=_input_spec(
            args, "uniprot-sprot-sequences", args.uniprot_release, "sprot"
        ),
        uniprot_trembl_sequences=_input_spec(
            args, "uniprot-trembl-sequences", args.uniprot_release, "trembl"
        ),
        goa=_input_spec(args, "goa", args.goa_release, "uniprotkb-goa"),
        go_obo=_input_spec(args, "go-obo", args.ontology_release, "gene-ontology"),
        split_policy=args.split_policy,
        training_population=args.training_population,
        mmseqs_bin=args.mmseqs_bin,
        expected_mmseqs_version=args.expected_mmseqs_version,
        cluster_assignments=args.cluster_assignments,
        frozen_input_manifest=args.frozen_input_manifest,
        attrition_policy=args.attrition_policy,
        attrition_override=args.attrition_override,
        framework_revision=args.framework_revision,
        fixture_mode=args.fixture_mode,
        diagnostic_pilot=args.diagnostic_pilot,
        threads=args.threads,
        requested_slots=args.requested_slots,
        allocated_slots=args.allocated_slots,
        run_id=args.run_id,
        seed=args.seed,
        min_count=args.min_count,
        development_fraction=args.development_fraction,
        training_fraction_within_development=args.training_fraction_within_development,
        sensitivity=args.sensitivity,
        evalue=args.evalue,
        include_relationships=not args.no_relationships,
        allow_downloads=not args.no_downloads,
        strict_qc=not args.no_strict_qc,
        allow_empty_fixture_outputs=args.allow_empty_fixture_outputs,
        keep_temp=args.keep_temp,
        release_uniprot=args.uniprot_release,
        release_goa=args.goa_release,
        release_ontology=args.ontology_release,
        scratch_safety_multiplier=args.scratch_safety_multiplier,
        minimum_free_disk_bytes=int(args.minimum_free_disk_gb * 1024 ** 3),
        persistent_results_root=args.persistent_results_root,
        mmseqs_work_multiplier=args.mmseqs_work_multiplier,
        publication_safety_multiplier=args.publication_safety_multiplier,
        excluded_sample_per_reason=args.excluded_sample_per_reason,
    )


def _preview(config: BuildConfig) -> dict[str, object]:
    config.validate(require_pinned_inputs=False)
    path = config.uniref90_fasta.path
    if path is None:
        path = Path(f"<download:{config.uniref90_fasta.url or 'UNRESOLVED_UNIREF90_URL'}>")
    commands = build_mmseqs_commands(config, path, config.temp_dir / config.identity_directory / "mmseqs")
    return {
        "identity": config.identity,
        "benchmark_scope": config.benchmark_scope,
        "uniprot_source_scope": config.uniprot_source_scope,
        "framework_revision": config.framework_revision,
        "attrition_policy": str(config.attrition_policy) if config.attrition_policy else None,
        "requested_slots": config.requested_slots,
        "allocated_slots": config.allocated_slots,
        "mmseqs_threads": config.threads,
        "final_output": str(
            config.output_dir / config.publication_relative_path
        ),
        "mmseqs_commands": [command.display for command in commands],
        "downloads_enabled": config.allow_downloads,
        "note": "Preview only; no path, release, hash, executable, or biological validation was performed.",
    }


def _resolve_run_dir(candidate: Path) -> Path:
    resolved = candidate.expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(resolved)
    if (resolved / "publication_metadata.json").is_file() and (
        resolved / "RUN_COMPLETE.json"
    ).is_file():
        return resolved
    matches = sorted({
        path.parent.resolve()
        for path in resolved.rglob("publication_metadata.json")
        if (path.parent / "RUN_COMPLETE.json").is_file()
        and json.loads(path.read_text(encoding="utf-8")).get("benchmark_scope")
        != "all-thresholds-summary"
    })
    if len(matches) != 1:
        raise ValueError(
            f"Run root must contain exactly one benchmark publication; {resolved} has {len(matches)}"
        )
    return matches[0]


def _cross_threshold_reports(
    runs,
    root: Path,
    split_policy: str | None = None,
    training_population: str | None = None,
    seed: int | None = None,
    min_count: int | None = None,
) -> Path:
    expected_percentages = {int(identity * 100) for identity in SUPPORTED_IDENTITIES}
    candidates = [item.output_dir if isinstance(item, BuildResult) else Path(item) for item in runs]
    if len(candidates) != 6:
        raise ValueError("Cross-threshold reporting requires exactly six --run-dir values")
    run_dirs = [_resolve_run_dir(path) for path in candidates]
    if len(set(run_dirs)) != 6:
        raise ValueError("Cross-threshold run directories must be six distinct publications")

    child_metadata: dict[int, dict[str, object]] = {}
    child_dirs: dict[int, Path] = {}
    for run_dir in run_dirs:
        validate_publication(run_dir)
        metadata = json.loads(
            (run_dir / "publication_metadata.json").read_text(encoding="utf-8")
        )
        if metadata.get("benchmark_scope") not in {"fixture-only", "dissertation-production"}:
            raise ValueError(f"Child run has invalid benchmark scope: {run_dir}")
        percent = int(metadata.get("identity_percent", -1))
        if percent in child_metadata:
            raise ValueError(f"Duplicate identity threshold in child publications: {percent}")
        child_metadata[percent] = metadata
        child_dirs[percent] = run_dir
    observed_percentages = set(child_metadata)
    if observed_percentages != expected_percentages:
        raise ValueError(
            "Cross-threshold reporting requires each locked identity exactly once; "
            f"observed={sorted(observed_percentages, reverse=True)}"
        )
    fingerprints = {str(item.get("scientific_fingerprint")) for item in child_metadata.values()}
    frozen_manifest_hashes = {
        str(item.get("frozen_input_manifest_sha256")) for item in child_metadata.values()
    }
    source_scopes = {str(item.get("uniprot_source_scope")) for item in child_metadata.values()}
    framework_revisions = {str(item.get("framework_revision")) for item in child_metadata.values()}
    attrition_policy_hashes = {
        str(item.get("attrition_policy_sha256")) for item in child_metadata.values()
    }
    if (
        len(fingerprints) != 1
        or len(frozen_manifest_hashes) != 1
        or len(source_scopes) != 1
        or len(framework_revisions) != 1
        or len(attrition_policy_hashes) != 1
    ):
        raise ValueError(
            "The six threshold runs do not share one source scope, framework revision, "
            "frozen manifest, attrition policy, and scientific fingerprint"
        )
    first = child_metadata[max(child_metadata)]
    fingerprint_payload = first.get("scientific_fingerprint_payload")
    if not isinstance(fingerprint_payload, dict):
        raise ValueError("Child publication lacks a validated scientific fingerprint payload")
    actual_split_policy = str(first["split_policy"])
    actual_population = str(first["training_population"])
    actual_seed = int(first["seed"])
    actual_min_count = int(first["min_count"])
    legacy_expectations = (
        ("split_policy", split_policy, actual_split_policy),
        ("training_population", training_population, actual_population),
        ("seed", seed, actual_seed),
        ("min_count", min_count, actual_min_count),
    )
    for label, expected, observed in legacy_expectations:
        if expected is not None and expected != observed:
            raise ValueError(f"Aggregate {label} selector does not match child metadata")
    production_values = {bool(item["production_eligible"]) for item in child_metadata.values()}
    fixture_values = {bool(item["fixture_mode"]) for item in child_metadata.values()}
    if len(production_values) != 1 or len(fixture_values) != 1:
        raise ValueError("Production/fixture eligibility differs across child publications")
    production_eligible = next(iter(production_values))
    fixture_mode = next(iter(fixture_values))
    source_scope = next(iter(source_scopes))
    framework_revision = next(iter(framework_revisions))
    attrition_policy_sha256 = next(iter(attrition_policy_hashes))
    aggregate_repository_state = git_state(Path(__file__).resolve().parents[4])
    if production_eligible and (
        aggregate_repository_state.get("commit") != first.get("repository_commit")
        or aggregate_repository_state.get("dirty") is not False
    ):
        raise ValueError(
            "Production aggregation requires the same clean framework commit as every child; "
            f"observed={aggregate_repository_state}"
        )
    final = (
        root.resolve() / f"source_{source_scope}" / f"framework_{framework_revision[:12]}"
        / "all_thresholds_summary" / actual_split_policy / actual_population
        / f"seed_{actual_seed}" / f"min_count_{actual_min_count}"
    )
    if any(final.is_relative_to(child) for child in child_dirs.values()):
        raise ValueError("Aggregate output must not be nested inside a child publication")
    with staging_output(final) as stage:
        universes: dict[int, dict[str, set[str]]] = {}
        for percent, run_dir in child_dirs.items():
            universes[percent] = {"all": set(pd.read_pickle(run_dir / "terms.pkl")["terms"])}
            by_namespace = {"biological_process": set(), "cellular_component": set(), "molecular_function": set()}
            with (run_dir / "go_term_summary.tsv").open("r", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row["in_development_universe"] == "1":
                        by_namespace[row["namespace"]].add(row["go_id"])
            universes[percent].update(by_namespace)

        overlap_rows = []
        percentages = sorted(universes, reverse=True)
        for left_index, left in enumerate(percentages):
            for right in percentages[left_index + 1:]:
                for namespace in sorted(universes[left]):
                    a, b = universes[left][namespace], universes[right][namespace]
                    union = a | b
                    overlap_rows.append({
                        "identity_a": left, "identity_b": right, "namespace": namespace,
                        "terms_a": len(a), "terms_b": len(b), "intersection": len(a & b),
                        "union": len(union), "jaccard": len(a & b) / len(union) if union else 1.0,
                    })
        with (stage / "term_universe_overlap.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(overlap_rows[0]) if overlap_rows else [
                "identity_a", "identity_b", "namespace", "terms_a", "terms_b", "intersection", "union", "jaccard"
            ], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(overlap_rows)

        change_rows = []
        for previous, current in zip(percentages, percentages[1:]):
            for namespace in sorted(universes[previous]):
                gained = universes[current][namespace] - universes[previous][namespace]
                lost = universes[previous][namespace] - universes[current][namespace]
                for term in sorted(gained):
                    change_rows.append({"from_identity": previous, "to_identity": current, "namespace": namespace, "change": "gained", "go_id": term})
                for term in sorted(lost):
                    change_rows.append({"from_identity": previous, "to_identity": current, "namespace": namespace, "change": "lost", "go_id": term})
        with (stage / "term_universe_changes.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["from_identity", "to_identity", "namespace", "change", "go_id"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(change_rows)
        assignment_hashes = {
            percent: sha256_file(run_dir / "mmseqs_cluster_membership.tsv.gz")
            for percent, run_dir in child_dirs.items()
        }
        identical_rows = []
        for left_index, left in enumerate(percentages):
            for right in percentages[left_index + 1:]:
                if assignment_hashes[left] == assignment_hashes[right]:
                    identical_rows.append({
                        "identity_a": left,
                        "identity_b": right,
                        "assignment_sha256": assignment_hashes[left],
                        "warning": "identical MMseqs2 cluster assignment manifests",
                    })
        with (stage / "identical_cluster_assignments.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["identity_a", "identity_b", "assignment_sha256", "warning"],
                delimiter="\t", lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(identical_rows)
        benchmark_rows = []
        split_rows = []
        attrition_rows = []
        for percent in percentages:
            run_dir = child_dirs[percent]
            benchmark = json.loads(
                (run_dir / "benchmark_summary.json").read_text(encoding="utf-8")
            )
            benchmark_rows.append({
                "identity_percent": percent,
                "uniprot_source_scope": source_scope,
                **benchmark["counts"],
            })
            with (run_dir / "split_summary.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    split_rows.append({"identity_percent": percent, **row})
            attrition = json.loads(
                (run_dir / "attrition_report.json").read_text(encoding="utf-8")
            )
            for metric in attrition["metrics"]:
                attrition_rows.append({
                    "identity_percent": percent,
                    "metric": metric["name"],
                    "numerator": metric["numerator"],
                    "denominator": metric["denominator"],
                    "ratio": metric["ratio"],
                    "allowed_ratio": metric["allowed_ratio"],
                    "bound_type": metric["bound_type"],
                    "passed": int(metric["passed"]),
                })
        for filename, rows in (
            ("threshold_metrics.tsv", benchmark_rows),
            ("threshold_split_metrics.tsv", split_rows),
            ("threshold_attrition_metrics.tsv", attrition_rows),
        ):
            with (stage / filename).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(rows[0]) if rows else [],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)
        (stage / "all_thresholds_summary.json").write_text(json.dumps({
            "identities": percentages,
            "uniprot_source_scope": source_scope,
            "framework_revision": framework_revision,
            "frozen_input_manifest_sha256": next(iter(frozen_manifest_hashes)),
            "attrition_policy_sha256": attrition_policy_sha256,
            "split_policy": actual_split_policy,
            "training_population": actual_population,
            "seed": actual_seed,
            "min_count": actual_min_count,
            "production_eligible": production_eligible,
            "scientific_fingerprint": next(iter(fingerprints)),
            "runs": {str(percent): str(path) for percent, path in sorted(child_dirs.items())},
            "term_counts": {
                str(identity): {namespace: len(terms) for namespace, terms in namespaces.items()}
                for identity, namespaces in universes.items()
            },
            "identical_cluster_assignment_pairs": identical_rows,
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        checks = [
            {
                "name": "all_six_identities_exactly_once",
                "passed": percentages == sorted(expected_percentages, reverse=True),
                "detail": "The aggregate contains exactly Daniel's six locked thresholds",
                "observed": percentages,
            },
            {
                "name": "child_publications_valid",
                "passed": len(child_metadata) == 6,
                "detail": "Every child run passed its manifest, marker, and semantic validation",
            },
            {
                "name": "cross_threshold_reports_present",
                "passed": all((stage / name).is_file() for name in (
                    "term_universe_overlap.tsv", "term_universe_changes.tsv",
                    "identical_cluster_assignments.tsv", "all_thresholds_summary.json",
                    "threshold_metrics.tsv", "threshold_split_metrics.tsv",
                    "threshold_attrition_metrics.tsv",
                )),
                "detail": "All required cross-threshold reports were written",
            },
        ]
        validation_payload = {
            "valid": all(item["passed"] for item in checks),
            "checks": checks,
            "warnings": ([{
                "name": "identical_cluster_assignments",
                "detail": "One or more thresholds produced identical normalized assignments",
                "pairs": identical_rows,
            }] if identical_rows else []),
        }
        (stage / "validation_report.json").write_text(
            json.dumps(validation_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (stage / "validation_report.md").write_text(
            "# Cross-threshold validation report\n\n"
            f"- Overall contract status: **{'PASS' if validation_payload['valid'] else 'FAIL'}**\n"
            f"- Identities: {', '.join(str(value) for value in percentages)}\n"
            f"- Production eligible: {str(production_eligible).lower()}\n",
            encoding="utf-8",
        )
        if not validation_payload["valid"]:
            raise ValueError("Cross-threshold validation failed before publication")
        parameter_keys = (
            "split_policy", "development_fraction", "training_fraction_within_development",
            "training_population", "seed", "min_count", "evidence_codes",
            "include_relationships", "root_policy", "coverage", "cov_mode", "cluster_mode",
            "alignment_mode", "seq_id_mode", "createdb_shuffle", "cluster_reassign",
            "sensitivity", "evalue", "uniprot_release", "goa_release", "ontology_release",
            "uniprot_source_scope",
        )
        parameters = {key: fingerprint_payload[key] for key in parameter_keys}
        (stage / "parameters.json").write_text(
            json.dumps(parameters, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        first_dir = child_dirs[max(child_dirs)]
        import shutil

        shutil.copyfile(
            first_dir / "frozen_input_manifest.json", stage / "frozen_input_manifest.json"
        )
        child_references = {
            str(percent): {
                "run_dir": str(run_dir),
                "output_manifest_sha256": sha256_file(run_dir / "output_manifest.json"),
                "publication_metadata_sha256": sha256_file(
                    run_dir / "publication_metadata.json"
                ),
            }
            for percent, run_dir in sorted(child_dirs.items())
        }
        input_manifest = {
            "schema_version": 1,
            "kind": "all-thresholds-child-references",
            "uniprot_source_scope": source_scope,
            "framework_revision": framework_revision,
            "frozen_input_manifest_sha256": next(iter(frozen_manifest_hashes)),
            "attrition_policy_sha256": attrition_policy_sha256,
            "children": child_references,
        }
        (stage / "input_manifest.json").write_text(
            json.dumps(input_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        runtime = {
            "repository": aggregate_repository_state,
            "runtime": {"mmseqs2": {
                "expected_version": first["expected_mmseqs_version"],
                "observed_version_token": first["observed_mmseqs_version"],
                "resolved_executable": first["mmseqs_resolved_executable"],
                "executable_sha256": first["mmseqs_executable_sha256"],
            }},
        }
        (stage / "run_provenance.json").write_text(
            json.dumps(runtime, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        publication_metadata = {
            "schema_version": 1,
            "fixture_mode": fixture_mode,
            "production_eligible": production_eligible,
            "benchmark_scope": "all-thresholds-summary",
            "uniprot_source_scope": source_scope,
            "framework_revision": framework_revision,
            "run_id": "aggregate",
            "identity_percent": None,
            "identities": percentages,
            "split_policy": actual_split_policy,
            "training_population": actual_population,
            "seed": actual_seed,
            "min_count": actual_min_count,
            "run_input_manifest_sha256": sha256_file(stage / "input_manifest.json"),
            "frozen_input_manifest_sha256": sha256_file(stage / "frozen_input_manifest.json"),
            "attrition_policy_sha256": attrition_policy_sha256,
            "attrition_report_sha256": None,
            "attrition_override_sha256": None,
            "attrition_policy_passed": all(
                bool(item.get("attrition_policy_passed")) for item in child_metadata.values()
            ),
            "attrition_override_valid": any(
                bool(item.get("attrition_override_valid")) for item in child_metadata.values()
            ),
            "requested_slots": first.get("requested_slots"),
            "allocated_slots": first.get("allocated_slots"),
            "mmseqs_threads": first.get("mmseqs_threads"),
            "expected_mmseqs_version": first["expected_mmseqs_version"],
            "observed_mmseqs_version": first["observed_mmseqs_version"],
            "mmseqs_resolved_executable": first["mmseqs_resolved_executable"],
            "mmseqs_executable_sha256": first["mmseqs_executable_sha256"],
            "repository_commit": first["repository_commit"],
            "scientific_fingerprint_payload": fingerprint_payload,
            "scientific_fingerprint": next(iter(fingerprints)),
        }
        from .provenance import write_publication_metadata

        write_publication_metadata(stage, publication_metadata)
        write_output_manifest(stage)
        publish(stage, final)
    try:
        verify_output_manifest(final)
        write_completion_marker(final)
        validate_publication(final)
    except BaseException:
        import shutil

        shutil.rmtree(final, ignore_errors=True)
        raise
    return final


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if args.command == "validate":
        validate_publication(args.run_dir)
        print(f"Valid published homology benchmark: {args.run_dir.resolve()}")
        return 0
    if args.command == "summarize":
        try:
            aggregate = _cross_threshold_reports(args.run_dir, args.output_dir)
            print(f"Published cross-threshold term analysis: {aggregate}")
            return 0
        except (OSError, ValueError, RuntimeError) as exc:
            parser.error(str(exc))
    if args.command == "authorize-array":
        try:
            pilot_run_dir = args.pilot_run_dir.resolve()
            if args.pilot_completion_marker.resolve() != pilot_run_dir / "RUN_COMPLETE.json":
                raise ValueError(
                    "Pilot completion marker must be RUN_COMPLETE.json in --pilot-run-dir"
                )
            if args.pilot_attrition_report.resolve() != pilot_run_dir / "attrition_report.json":
                raise ValueError(
                    "Pilot attrition report must be attrition_report.json in --pilot-run-dir"
                )
            validate_publication(pilot_run_dir)
            manifest_hash = sha256_file(args.frozen_input_manifest)
            reviewed_policy, reviewed_policy_sha256 = load_attrition_policy(
                args.attrition_policy,
                source_scope=args.uniprot_source_scope,
                expected_releases={
                    "uniprot_uniref": args.uniprot_release,
                    "goa": args.goa_release,
                    "ontology": args.ontology_release,
                },
                framework_commit=args.framework_revision,
                frozen_input_manifest_sha256=manifest_hash,
            )
            validate_pilot_approval(
                args.pilot_approval,
                completion_marker_path=args.pilot_completion_marker,
                attrition_report_path=args.pilot_attrition_report,
                task_context_path=args.pilot_task_context,
                measurement_evidence_path=args.pilot_measurement_evidence,
                framework_commit=args.framework_revision,
                frozen_input_manifest_sha256=manifest_hash,
                source_scope=args.uniprot_source_scope,
                split_policy=args.split_policy,
                training_population=args.training_population,
                mmseqs_version=args.expected_mmseqs_version,
                reviewed_attrition_policy=reviewed_policy,
                reviewed_attrition_policy_sha256=reviewed_policy_sha256,
            )
            print("Reviewed pilot approval and attrition policy authorize array preview/submission")
            return 0
        except (OSError, ValueError, RuntimeError) as exc:
            parser.error(str(exc))

    logging.basicConfig(level=getattr(logging, args.verbosity), format="%(levelname)s %(message)s")
    try:
        identities = _identities(args.identity)
        configs = [_config(args, identity) for identity in identities]
        if args.dry_run:
            print(json.dumps([_preview(config) for config in configs], indent=2, sort_keys=True))
            return 0
        for config in configs:
            config.validate()
        if len(configs) > 1:
            missing_local = [
                spec.name
                for spec in (
                    configs[0].uniref90_fasta,
                    configs[0].idmapping,
                    *(getattr(configs[0], name) for name in configs[0].selected_uniprot_input_names),
                    configs[0].goa,
                    configs[0].go_obo,
                )
                if spec is not None
                if spec.path is None
            ]
            if missing_local:
                raise ValueError(
                    "--identity all requires all frozen inputs to be staged as local paths so "
                    "multi-gigabyte sources are not downloaded once per threshold; missing local="
                    + ",".join(missing_local)
                )
        results = []
        for config in configs:
            LOGGER.info("Building identity %d%%", int(config.identity * 100))
            result = build_benchmark(config)
            results.append(result)
            print(f"Published validated benchmark: {result.output_dir}")
        if len(results) == len(SUPPORTED_IDENTITIES):
            aggregate = _cross_threshold_reports(results, args.output_dir)
            print(f"Published cross-threshold term analysis: {aggregate}")
        return 0
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
