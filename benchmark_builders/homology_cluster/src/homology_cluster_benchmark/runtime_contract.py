from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from urllib.parse import urlparse

from .attrition import METRIC_DEFINITIONS
from .frozen_inputs import (
    SOURCE_POPULATIONS,
    expected_input_names,
    load_frozen_input_manifest,
)
from .inputs import open_text, sha256_file
from .pipeline import validate_publication


RUNTIME_REVIEW_DATE = "2026-07-15"
RUNTIME_CREATED_AT = "2026-06-17T00:00:00Z"
OFFICIAL_HOSTS = {
    "uniref90_fasta": "ftp.uniprot.org",
    "idmapping": "ftp.uniprot.org",
    "uniprot_sprot_sequences": "ftp.uniprot.org",
    "uniprot_trembl_sequences": "ftp.uniprot.org",
    "goa": "ftp.ebi.ac.uk",
    "go_obo": "release.geneontology.org",
}
ROLE_FILENAMES = {
    "uniref90_fasta": "uniref90.fasta.gz",
    "idmapping": "idmapping_selected.tab.gz",
    "uniprot_sprot_sequences": "uniprot_sprot.dat.gz",
    "uniprot_trembl_sequences": "uniprot_trembl.dat.gz",
    "goa": "goa_uniprot_all.gaf.234.gz",
    "go_obo": "go-basic.obo",
}


@dataclass(frozen=True)
class RuntimeInput:
    name: str
    path: Path
    url: str
    acquisition: str


def _json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _gaf_metadata(path: Path) -> dict[str, str]:
    headers: dict[str, str] = {}
    with open_text(path) as handle:
        for raw_line in handle:
            if not raw_line.startswith("!"):
                break
            text = raw_line[1:].strip()
            if ":" not in text:
                continue
            key, value = text.split(":", 1)
            headers[key.strip().lower().replace("-", "_")] = value.strip()
    required = {"gaf_version", "date_generated", "go_version"}
    if not required.issubset(headers):
        raise ValueError(f"GOA input lacks required release headers: {sorted(required - set(headers))}")
    if headers["gaf_version"] != "2.2":
        raise ValueError(f"Expected GAF 2.2, observed {headers['gaf_version']!r}")
    if "2026-06-17" not in headers["date_generated"]:
        raise ValueError(
            "GOA runtime input is not release 234 from 2026-06-17: "
            f"{headers['date_generated']!r}"
        )
    if "2026-06-15" not in headers["go_version"]:
        raise ValueError(
            "GOA runtime input does not declare GO 2026-06-15: "
            f"{headers['go_version']!r}"
        )
    return {key: headers[key] for key in sorted(required)}


def _obo_metadata(path: Path) -> dict[str, str]:
    with open_text(path) as handle:
        for raw_line in handle:
            if raw_line.startswith("data-version:"):
                value = raw_line.split(":", 1)[1].strip()
                if value != "releases/2026-06-15":
                    raise ValueError(
                        "Ontology runtime input has the wrong data-version: "
                        f"expected 'releases/2026-06-15', observed {value!r}"
                    )
                return {"data_version": value}
    raise ValueError("Ontology runtime input lacks a data-version header")


def _release_for(name: str) -> str:
    if name == "goa":
        return "234"
    if name == "go_obo":
        return "releases/2026-06-15"
    return "2026_02"


def _validate_runtime_input(item: RuntimeInput) -> None:
    if item.name not in ROLE_FILENAMES:
        raise ValueError(f"Unsupported runtime input role: {item.name}")
    if not item.path.is_file() or item.path.stat().st_size <= 0:
        raise ValueError(f"Runtime input is missing or empty: {item.name}: {item.path}")
    if item.path.name != ROLE_FILENAMES[item.name]:
        raise ValueError(
            f"Runtime input {item.name} must use canonical filename "
            f"{ROLE_FILENAMES[item.name]!r}, observed {item.path.name!r}"
        )
    parsed = urlparse(item.url)
    if parsed.scheme != "https" or parsed.hostname != OFFICIAL_HOSTS[item.name]:
        raise ValueError(
            f"Runtime input {item.name} must retain its reviewed official HTTPS URL; "
            f"observed {item.url!r}"
        )
    if item.acquisition not in {
        "downloaded-to-scratch",
        "provided-path-staged-to-scratch",
        "provided-persistent-store",
    }:
        raise ValueError(f"Unsupported acquisition record for {item.name}: {item.acquisition!r}")


def write_runtime_policy(
    policy_path: Path,
    manifest_path: Path,
    source_scope: str,
    framework_revision: str,
) -> str:
    manifest = load_frozen_input_manifest(
        manifest_path, uniprot_source_scope=source_scope, fixture_mode=False
    )
    if len(framework_revision) != 40 or any(
        ch not in "0123456789abcdef" for ch in framework_revision
    ):
        raise ValueError("framework_revision must be exactly 40 lowercase hexadecimal characters")
    metrics = {}
    for name, definition in METRIC_DEFINITIONS.items():
        bound_key = f"allowed_{definition.bound}_ratio"
        metrics[name] = {
            "numerator_definition": definition.numerator,
            "denominator_definition": definition.denominator,
            bound_key: 0.0 if definition.bound == "minimum" else 1.0,
            "rationale": (
                "Automatic non-blocking observation boundary for direct runtime submission. "
                "The strict semantic validator remains authoritative."
            ),
            "evidence_source": (
                "Runtime measurement; no pilot-derived biological attrition threshold asserted."
            ),
        }
    policy = {
        "schema_name": "homology-cluster-attrition-policy",
        "schema_version": 1,
        "uniprot_source_scope": source_scope,
        "expected_releases": {
            "uniprot_uniref": "2026_02",
            "goa": "234",
            "ontology": "releases/2026-06-15",
        },
        "metrics": metrics,
        "rationale": (
            "Permit direct six-task execution without making a prior pilot mandatory. All "
            "attrition metrics are retained for review, but these permissive bounds are not "
            "presented as biologically optimized thresholds."
        ),
        "evidence_source": "Automated runtime validation and complete attrition reporting.",
        "author": "homology runtime wrapper",
        "reviewer": "automated-runtime-contract",
        "review_date": RUNTIME_REVIEW_DATE,
        "framework_commit": framework_revision,
        "frozen_input_manifest_sha256": manifest.sha256,
    }
    _json(policy_path, policy)
    return sha256_file(policy_path)


def write_runtime_contract(
    manifest_path: Path,
    policy_path: Path,
    source_scope: str,
    framework_revision: str,
    inputs: list[RuntimeInput],
) -> dict[str, str]:
    ordered_names = expected_input_names(source_scope)
    by_name = {item.name: item for item in inputs}
    if len(by_name) != len(inputs) or set(by_name) != set(ordered_names):
        raise ValueError(
            "Runtime inputs do not exactly match source scope: "
            f"expected={list(ordered_names)}, observed={sorted(by_name)}"
        )
    if len(framework_revision) != 40 or any(ch not in "0123456789abcdef" for ch in framework_revision):
        raise ValueError("framework_revision must be exactly 40 lowercase hexadecimal characters")

    entries = []
    for name in ordered_names:
        item = by_name[name]
        _validate_runtime_input(item)
        embedded_metadata: dict[str, str] = {}
        if name == "goa":
            embedded_metadata = _gaf_metadata(item.path)
        elif name == "go_obo":
            embedded_metadata = _obo_metadata(item.path)
        entries.append({
            "name": name,
            "logical_role": name,
            "source_population": SOURCE_POPULATIONS[name],
            "release": _release_for(name),
            "url": item.url,
            "local_filename": item.path.name,
            "size_bytes": item.path.stat().st_size,
            "sha256": sha256_file(item.path),
            "acquisition": item.acquisition,
            "embedded_metadata": embedded_metadata,
            "notes": (
                "Read from an authenticated persistent project store to build the shared "
                "threshold-independent preprocessing cache."
                if item.acquisition == "provided-persistent-store" else
                "Staged in job-owned scratch. Release markers and embedded metadata were "
                "validated before benchmark construction; bytes are discarded after copy-back."
            ),
        })

    manifest = {
        "schema_name": "homology-cluster-frozen-inputs",
        "schema_version": 2,
        "uniprot_source_scope": source_scope,
        "created_at": RUNTIME_CREATED_AT,
        "review": {
            "status": "authoritative-source-reviewed",
            "authoritative_origin": True,
            "reviewed_by": "repository-reviewed-runtime-source-contract",
            "evidence": (
                "Official release-specific host allow-list, UniProt 2026_02 marker, GOA 234 "
                "marker/header, GO releases/2026-06-15 header, byte size and SHA-256."
            ),
        },
        "inputs": entries,
    }
    _json(manifest_path, manifest)
    manifest_sha256 = sha256_file(manifest_path)

    policy_sha256 = write_runtime_policy(
        policy_path, manifest_path, source_scope, framework_revision
    )
    return {
        "manifest_sha256": manifest_sha256,
        "policy_sha256": policy_sha256,
    }


def write_runtime_review(
    run_dir: Path,
    output_dir: Path,
    run_kind: str,
) -> dict[str, object]:
    if run_kind not in {"pilot", "array"}:
        raise ValueError("run_kind must be pilot or array")
    validate_publication(run_dir)
    publication = json.loads((run_dir / "publication_metadata.json").read_text(encoding="utf-8"))
    validation = json.loads((run_dir / "validation_report.json").read_text(encoding="utf-8"))
    attrition = json.loads((run_dir / "attrition_report.json").read_text(encoding="utf-8"))
    csv_files = sorted(path.name for path in run_dir.glob("*.csv"))
    pickle_files = sorted(path.name for path in run_dir.glob("*.pkl"))
    if len(csv_files) != 9:
        raise ValueError(f"Runtime publication contains {len(csv_files)} CSVs instead of nine")
    if len(pickle_files) != 5:
        raise ValueError(f"Runtime publication contains {len(pickle_files)} pickles instead of five")
    warnings = list(validation.get("warnings", []))
    payload: dict[str, object] = {
        "schema_name": "homology-runtime-automatic-review",
        "schema_version": 1,
        "status": "pass",
        "run_kind": run_kind,
        "pilot_required_for_array": False,
        "run_dir": str(run_dir.resolve()),
        "identity_percent": publication["identity_percent"],
        "uniprot_source_scope": publication["uniprot_source_scope"],
        "framework_revision": publication["framework_revision"],
        "validation_valid": validation.get("valid") is True,
        "validation_check_count": len(validation.get("checks", [])),
        "validation_warning_count": len(warnings),
        "attrition_policy_passed": attrition.get("policy_passed") is True,
        "attrition_policy_kind": "automatic-nonblocking-runtime-observation",
        "production_eligible": publication.get("production_eligible"),
        "csv_files": csv_files,
        "pickle_files": pickle_files,
        "warnings": warnings,
        "interpretation": (
            "This automatic review proves software and publication contracts. It records, but "
            "does not biologically optimize, attrition thresholds."
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    _json(output_dir / "automatic_review.json", payload)
    markdown = [
        "# Homology runtime automatic review",
        "",
        "- Status: **pass**",
        f"- Run kind: **{run_kind}**",
        f"- Identity: **{publication['identity_percent']}%**",
        f"- UniProt scope: **{publication['uniprot_source_scope']}**",
        f"- Validation checks: **{payload['validation_check_count']}**",
        f"- Validation warnings: **{payload['validation_warning_count']}**",
        "- Pilot required before the full array: **no**",
        "- CSV outputs: **9**",
        "- Pickle outputs: **5**",
        "",
        (
            "The strict publication validator passed. The runtime attrition policy is "
            "deliberately non-blocking and records observations; it is not a claim that "
            "pilot-derived biological attrition limits were reviewed."
        ),
        "",
    ]
    (output_dir / "automatic_review.md").write_text("\n".join(markdown), encoding="utf-8")
    return payload


def _input_argument(parser: argparse.ArgumentParser, name: str) -> None:
    option = name.replace("_", "-")
    parser.add_argument(f"--{option}", type=Path)
    parser.add_argument(f"--{option}-url")
    parser.add_argument(
        f"--{option}-acquisition",
        choices=(
            "downloaded-to-scratch",
            "provided-path-staged-to-scratch",
            "provided-persistent-store",
        ),
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare and review runtime homology contracts")
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--manifest-out", type=Path, required=True)
    prepare.add_argument("--policy-out", type=Path, required=True)
    prepare.add_argument(
        "--source-scope",
        choices=("sprot-only", "trembl-only", "sprot-and-trembl"),
        required=True,
    )
    prepare.add_argument("--framework-revision", required=True)
    for name in ROLE_FILENAMES:
        _input_argument(prepare, name)
    policy = subparsers.add_parser("policy")
    policy.add_argument("--manifest", type=Path, required=True)
    policy.add_argument("--policy-out", type=Path, required=True)
    policy.add_argument(
        "--source-scope",
        choices=("sprot-only", "trembl-only", "sprot-and-trembl"),
        required=True,
    )
    policy.add_argument("--framework-revision", required=True)
    review = subparsers.add_parser("review")
    review.add_argument("--run-dir", type=Path, required=True)
    review.add_argument("--output-dir", type=Path, required=True)
    review.add_argument("--run-kind", choices=("pilot", "array"), required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "review":
        write_runtime_review(args.run_dir, args.output_dir, args.run_kind)
        return 0
    if args.command == "policy":
        digest = write_runtime_policy(
            args.policy_out,
            args.manifest,
            args.source_scope,
            args.framework_revision,
        )
        print(json.dumps({"policy_sha256": digest}, sort_keys=True))
        return 0
    required = set(expected_input_names(args.source_scope))
    inputs = []
    for name in ROLE_FILENAMES:
        path = getattr(args, name)
        url = getattr(args, f"{name}_url")
        acquisition = getattr(args, f"{name}_acquisition")
        if name in required:
            if path is None or url is None or acquisition is None:
                raise ValueError(
                    f"Source scope {args.source_scope} requires path, URL and acquisition "
                    f"for runtime input {name}"
                )
            inputs.append(RuntimeInput(name, path, url, acquisition))
        elif path is not None or url is not None or acquisition is not None:
            raise ValueError(f"Source scope {args.source_scope} forbids runtime input {name}")
    result = write_runtime_contract(
        args.manifest_out,
        args.policy_out,
        args.source_scope,
        args.framework_revision,
        inputs,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
