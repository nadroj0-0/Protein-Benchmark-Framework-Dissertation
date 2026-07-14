from __future__ import annotations

import csv
from dataclasses import replace
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.cli import _cross_threshold_reports
from homology_cluster_benchmark.config import SUPPORTED_IDENTITIES
from homology_cluster_benchmark.models import InputSpec
from homology_cluster_benchmark.pipeline import build_benchmark, validate_publication
from homology_cluster_benchmark.inputs import sha256_file
from homology_cluster_benchmark.provenance import PUBLICATION_MARKER_KEYS

from tests.helpers import FIXTURES, fixture_config


REQUIRED_MANIFESTS = {
    "input_manifest.json", "parameters.json", "run_provenance.json", "mmseqs_commands.tsv",
    "frozen_input_manifest.json", "publication_metadata.json", "uniprot_to_uniref90.tsv",
    "mmseqs_cluster_membership.tsv.gz",
    "protein_cluster_assignments.tsv", "cluster_split_assignments.tsv", "retained_clusters.tsv",
    "retained_cluster_members.tsv.gz", "qualifying_annotations.tsv.gz",
    "excluded_annotations_sample.tsv.gz", "annotation_decision_counts.tsv",
    "attrition_summary.json", "attrition_summary.tsv", "split_balance_summary.json",
    "mapping_summary.json", "evidence_summary.tsv", "go_term_summary.tsv", "taxonomy_summary.tsv",
    "split_summary.tsv", "cluster_size_summary.tsv", "benchmark_summary.json", "benchmark_summary.md",
    "disk_preflight.json", "validation_report.json", "validation_report.md",
    "output_manifest.json", "RUN_COMPLETE.json",
}
REQUIRED_CSVS = {
    f"{ontology}-{split}.csv"
    for ontology in ("bp", "cc", "mf")
    for split in ("training", "validation", "test")
}
REQUIRED_PICKLES = {
    "train_data.pkl", "train_data_train.pkl", "train_data_valid.pkl", "test_data.pkl", "terms.pkl",
}


def _rewrite_hashed_publication_metadata(run_dir: Path, publication: dict) -> None:
    publication_path = run_dir / "publication_metadata.json"
    publication_path.write_text(
        json.dumps(publication, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest_path = run_dir / "output_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    entry = next(
        item for item in manifest["files"] if item["path"] == "publication_metadata.json"
    )
    entry["size_bytes"] = publication_path.stat().st_size
    entry["sha256"] = sha256_file(publication_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    marker_path = run_dir / "RUN_COMPLETE.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    marker["publication_metadata_sha256"] = sha256_file(publication_path)
    marker["manifest_sha256"] = sha256_file(manifest_path)
    for key in PUBLICATION_MARKER_KEYS:
        marker[key] = publication.get(key)
    marker_path.write_text(
        json.dumps(marker, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _rewrite_fixture_repository_commit(run_dir: Path, commit: str) -> None:
    provenance_path = run_dir / "run_provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["repository"]["commit"] = commit
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    publication_path = run_dir / "publication_metadata.json"
    publication = json.loads(publication_path.read_text(encoding="utf-8"))
    publication["repository_commit"] = commit
    publication["framework_revision"] = commit
    fingerprint_payload = publication["scientific_fingerprint_payload"]
    fingerprint_payload["repository_commit"] = commit
    fingerprint_payload["framework_revision"] = commit
    publication["scientific_fingerprint"] = hashlib.sha256(
        json.dumps(
            fingerprint_payload, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()
    publication_path.write_text(
        json.dumps(publication, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest_path = run_dir / "output_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for relative, path in (
        ("publication_metadata.json", publication_path),
        ("run_provenance.json", provenance_path),
    ):
        entry = next(item for item in manifest["files"] if item["path"] == relative)
        entry["size_bytes"] = path.stat().st_size
        entry["sha256"] = sha256_file(path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    marker_path = run_dir / "RUN_COMPLETE.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    marker["publication_metadata_sha256"] = sha256_file(publication_path)
    marker["manifest_sha256"] = sha256_file(manifest_path)
    for key in PUBLICATION_MARKER_KEYS:
        marker[key] = publication.get(key)
    marker_path.write_text(
        json.dumps(marker, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


class PipelineTests(unittest.TestCase):
    def test_end_to_end_fixture_publishes_complete_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = build_benchmark(fixture_config(root / "outputs", root / "temp"))
            names = {path.name for path in result.files}
            self.assertTrue(REQUIRED_MANIFESTS | REQUIRED_CSVS | REQUIRED_PICKLES <= names)
            validate_publication(result.output_dir)
            marker = json.loads((result.output_dir / "RUN_COMPLETE.json").read_text())
            self.assertTrue(marker["complete"])
            self.assertTrue(marker["post_publication_hash_verification"])
            self.assertFalse(marker["production_eligible"])
            self.assertEqual(marker["benchmark_scope"], "fixture-only")
            self.assertTrue((result.output_dir / "logs" / "mmseqs" / "NOT_EXECUTED.json").is_file())
            validation = json.loads((result.output_dir / "validation_report.json").read_text())
            self.assertTrue(validation["valid"])

            terms = pd.read_pickle(result.output_dir / "terms.pkl")["terms"].tolist()
            self.assertTrue({"GO:0008150", "GO:0005575", "GO:0003674"} <= set(terms))
            term_index = {term: index for index, term in enumerate(terms)}
            for root_term in ("GO:0008150", "GO:0005575", "GO:0003674"):
                self.assertIsInstance(term_index[root_term], int)
            for prefix, root_term in (
                ("bp", "GO:0008150"), ("cc", "GO:0005575"),
                ("mf", "GO:0003674"),
            ):
                for split in ("training", "validation", "test"):
                    self.assertIn(
                        root_term,
                        pd.read_csv(result.output_dir / f"{prefix}-{split}.csv").columns,
                    )
            for name in REQUIRED_CSVS:
                frame = pd.read_csv(result.output_dir / name)
                self.assertEqual(frame.columns[:2].tolist(), ["proteins", "sequences"])
                self.assertGreater(len(frame.columns), 2)
                self.assertGreater(len(frame), 0)
            for name in REQUIRED_PICKLES - {"terms.pkl"}:
                self.assertEqual(
                    list(pd.read_pickle(result.output_dir / name).columns),
                    ["proteins", "sequences", "annotations"],
                )

            mapping = json.loads((result.output_dir / "mapping_summary.json").read_text())
            self.assertEqual(mapping["mapping_status_counts"]["ambiguous"], 1)
            self.assertEqual(mapping["mapping_status_counts"]["unmapped-absent"], 1)
            summary = json.loads((result.output_dir / "benchmark_summary.json").read_text())
            self.assertEqual(summary["counts"]["retained_mmseqs_clusters"], 5)
            retained = pd.read_csv(result.output_dir / "retained_cluster_members.tsv.gz", sep="\t")
            self.assertNotIn("UniRef90_U6", set(retained["uniref90_id"]))

    def test_atomic_failure_cleans_stage_and_never_writes_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bad = root / "bad.tsv"
            cluster_lines = [
                line for line in (FIXTURES / "clusters.tsv").read_text().splitlines() if line
            ]
            bad.write_text("\n".join(cluster_lines[:-1]) + "\n")
            config = fixture_config(root / "outputs", root / "temp", cluster_assignments=bad)
            with self.assertRaisesRegex(ValueError, "missing_members"):
                build_benchmark(config)
            final = (
                root / "outputs" / "identity_30" / "sequence-balanced" / "annotated-only"
                / "seed_0" / "min_count_1"
            )
            self.assertFalse(final.exists())
            self.assertFalse(any(root.glob("outputs/**/RUN_COMPLETE.json")))
            self.assertFalse(list(root.glob("outputs/**/.*.staging-*")))

    def test_deterministic_payload_hashes_match_across_reruns(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = build_benchmark(fixture_config(root / "out-a", root / "temp"))
            second = build_benchmark(fixture_config(root / "out-b", root / "temp"))
            first_manifest = json.loads((first.output_dir / "output_manifest.json").read_text())
            second_manifest = json.loads((second.output_dir / "output_manifest.json").read_text())
            first_hashes = {
                item["path"]: item["sha256"] for item in first_manifest["files"]
                if item["deterministic_payload"]
            }
            second_hashes = {
                item["path"]: item["sha256"] for item in second_manifest["files"]
                if item["deterministic_payload"]
            }
            self.assertEqual(first_hashes, second_hashes)

    def test_exact_sequence_cross_split_leak_fails_before_publication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = build_benchmark(fixture_config(root / "baseline", root / "temp"))
            with (baseline.output_dir / "protein_cluster_assignments.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            train_id = next(row["uniprot_accession"] for row in rows if row["split"] == "training")
            test_id = next(row["uniprot_accession"] for row in rows if row["split"] == "test")
            sequences = {}
            headers = {}
            for line in (FIXTURES / "uniprot.fasta").read_text().splitlines():
                if not line:
                    continue
                if line.startswith(">"):
                    protein = line.split("|")[1]
                    headers[protein] = line
                else:
                    sequences[protein] = line
            sequences[test_id] = sequences[train_id]
            modified = root / "duplicate.fasta"
            modified.write_text("".join(
                f"{headers[protein]}\n{sequences[protein]}\n" for protein in headers
            ))
            config = fixture_config(
                root / "leak", root / "temp",
                uniprot_sprot_sequences=InputSpec(
                    "uniprot_sprot_sequences", modified, release="2026_02",
                    source_population="sprot",
                ),
            )
            with self.assertRaisesRegex(ValueError, "global_exact_sequence_disjointness"):
                build_benchmark(config)
            self.assertFalse(
                (
                    root / "leak" / "identity_30" / "sequence-balanced" / "annotated-only"
                    / "seed_0" / "min_count_1"
                ).exists()
            )

    def test_uniprot_to_retained_uniref_cross_population_sequence_leak_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = build_benchmark(fixture_config(root / "baseline", root / "temp"))
            retained = pd.read_csv(
                baseline.output_dir / "retained_cluster_members.tsv.gz", sep="\t"
            )
            scaffold_id = retained.loc[
                retained["split"] == "training", "uniref90_id"
            ].iloc[0]
            mappings = pd.read_csv(
                baseline.output_dir / "protein_cluster_assignments.tsv", sep="\t"
            )
            protein_id = mappings.loc[
                (mappings["split"] == "test") & (mappings["mapping_status"] == "mapped"),
                "uniprot_accession",
            ].iloc[0]

            uniref_sequences = {}
            current = None
            for line in (FIXTURES / "uniref90.fasta").read_text().splitlines():
                if line.startswith(">"):
                    current = line[1:].split()[0]
                elif line:
                    uniref_sequences[current] = line
            headers, sequences = {}, {}
            for line in (FIXTURES / "uniprot.fasta").read_text().splitlines():
                if line.startswith(">"):
                    current = line.split("|")[1]
                    headers[current] = line
                elif line:
                    sequences[current] = line
            sequences[protein_id] = uniref_sequences[scaffold_id]
            modified = root / "cross-population.fasta"
            modified.write_text("".join(
                f"{headers[protein]}\n{sequences[protein]}\n" for protein in headers
            ))
            config = fixture_config(
                root / "cross-population", root / "temp",
                uniprot_sprot_sequences=InputSpec(
                    "uniprot_sprot_sequences", modified, release="2026_02",
                    source_population="sprot",
                ),
            )
            with self.assertRaisesRegex(ValueError, "global_retained_exact_sequence_disjointness"):
                build_benchmark(config)

    def test_manifest_extra_file_and_marker_tamper_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = build_benchmark(fixture_config(root / "outputs", root / "temp"))
            extra = result.output_dir / "unexpected.txt"
            extra.write_text("not listed\n")
            with self.assertRaisesRegex(ValueError, "does not reconcile"):
                validate_publication(result.output_dir)
            extra.unlink()
            marker_path = result.output_dir / "RUN_COMPLETE.json"
            marker = json.loads(marker_path.read_text())
            marker["manifest_sha256"] = "0" * 64
            marker_path.write_text(json.dumps(marker))
            with self.assertRaisesRegex(ValueError, "does not agree"):
                validate_publication(result.output_dir)

    def test_publication_policy_marker_forgery_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = build_benchmark(fixture_config(root / "outputs", root / "temp"))
            marker_path = result.output_dir / "RUN_COMPLETE.json"
            original = json.loads(marker_path.read_text())
            mutations = {
                "fixture_mode": False,
                "production_eligible": True,
                "benchmark_scope": "dissertation-production",
            }
            for key, value in mutations.items():
                with self.subTest(key=key):
                    changed = dict(original)
                    changed[key] = value
                    marker_path.write_text(json.dumps(changed))
                    with self.assertRaisesRegex(ValueError, "metadata mismatch"):
                        validate_publication(result.output_dir)
                    marker_path.write_text(json.dumps(original))
            changed = {**original, **mutations}
            marker_path.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "metadata mismatch"):
                validate_publication(result.output_dir)

    def test_attrition_report_is_bound_to_the_published_run_input_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = build_benchmark(fixture_config(root / "outputs", root / "temp"))
            attrition_path = result.output_dir / "attrition_report.json"
            attrition = json.loads(attrition_path.read_text(encoding="utf-8"))
            attrition["input_manifest_sha256"] = "0" * 64
            attrition_path.write_text(
                json.dumps(attrition, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            manifest_path = result.output_dir / "output_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            entry = next(
                item for item in manifest["files"]
                if item["path"] == "attrition_report.json"
            )
            entry["size_bytes"] = attrition_path.stat().st_size
            entry["sha256"] = sha256_file(attrition_path)
            manifest_path.write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            publication_path = result.output_dir / "publication_metadata.json"
            publication = json.loads(publication_path.read_text(encoding="utf-8"))
            publication["attrition_report_sha256"] = sha256_file(attrition_path)
            _rewrite_hashed_publication_metadata(result.output_dir, publication)
            with self.assertRaisesRegex(ValueError, "run-input-manifest binding"):
                validate_publication(result.output_dir)

    def test_malformed_review_policy_fails_before_large_input_hashing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = replace(
                fixture_config(root / "outputs", root / "temp"),
                frozen_input_manifest=root / "frozen.json",
                attrition_policy=root / "policy.json",
            )
            with (
                mock.patch(
                    "homology_cluster_benchmark.pipeline.load_frozen_input_manifest",
                    return_value=mock.Mock(sha256="a" * 64),
                ),
                mock.patch(
                    "homology_cluster_benchmark.pipeline.load_attrition_policy",
                    side_effect=ValueError("template placeholder"),
                ),
                mock.patch(
                    "homology_cluster_benchmark.pipeline._resolved_inputs"
                ) as resolve_inputs,
            ):
                with self.assertRaisesRegex(ValueError, "template placeholder"):
                    build_benchmark(config)
            resolve_inputs.assert_not_called()

    def test_rehashed_scope_and_missing_fingerprint_metadata_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            scoped = build_benchmark(fixture_config(root / "scope", root / "temp"))
            publication = json.loads(
                (scoped.output_dir / "publication_metadata.json").read_text()
            )
            publication["benchmark_scope"] = "dissertation-production"
            _rewrite_hashed_publication_metadata(scoped.output_dir, publication)
            with self.assertRaisesRegex(ValueError, "benchmark_scope"):
                validate_publication(scoped.output_dir)

            missing = build_benchmark(fixture_config(root / "fingerprint", root / "temp"))
            publication = json.loads(
                (missing.output_dir / "publication_metadata.json").read_text()
            )
            publication.pop("scientific_fingerprint")
            _rewrite_hashed_publication_metadata(missing.output_dir, publication)
            with self.assertRaisesRegex(ValueError, "missing required keys"):
                validate_publication(missing.output_dir)

    def test_independent_term_recount_rejects_reauthored_terms_tamper(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for mode in ("missing", "unexpected"):
                result = build_benchmark(
                    fixture_config(root / mode, root / "temp")
                )
                terms_path = result.output_dir / "terms.pkl"
                terms = pd.read_pickle(terms_path)
                if mode == "missing":
                    terms = terms.iloc[1:].reset_index(drop=True)
                else:
                    terms = pd.concat(
                        [terms, pd.DataFrame({"terms": ["GO:9999999"]})], ignore_index=True
                    ).sort_values("terms", kind="stable").reset_index(drop=True)
                terms.to_pickle(terms_path)
                manifest_path = result.output_dir / "output_manifest.json"
                manifest = json.loads(manifest_path.read_text())
                entry = next(item for item in manifest["files"] if item["path"] == "terms.pkl")
                entry["size_bytes"] = terms_path.stat().st_size
                entry["sha256"] = sha256_file(terms_path)
                manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
                marker_path = result.output_dir / "RUN_COMPLETE.json"
                marker = json.loads(marker_path.read_text())
                marker["manifest_sha256"] = sha256_file(manifest_path)
                marker_path.write_text(json.dumps(marker, indent=2, sort_keys=True) + "\n")
                with self.assertRaisesRegex(ValueError, "independent development recount"):
                    validate_publication(result.output_dir)

    def test_bounded_exclusion_sample_keeps_complete_counts_and_deterministic_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            goa = root / "many-rejections.gaf"
            extra = "".join(
                f"UniProtKB\tPX{index}\tPX{index}\tinvolved_in\tGO:0009987\tPMID:9\tIEA\t\tP\t\t\tprotein\ttaxon:9606\t20260617\tUniProt\t\t\n"
                for index in range(25)
            )
            goa.write_text((FIXTURES / "goa.gaf").read_text() + extra)
            outputs = []
            for name in ("a", "b"):
                result = build_benchmark(fixture_config(
                    root / name, root / "temp",
                    goa=InputSpec(
                        "goa", goa, release="234", source_population="uniprotkb-goa"
                    ),
                    excluded_sample_per_reason=3,
                ))
                outputs.append(result.output_dir)
                with gzip.open(
                    result.output_dir / "excluded_annotations_sample.tsv.gz", "rt"
                ) as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(
                    sum(row["rejection_reason"] == "evidence_code" for row in rows), 3
                )
                counts = pd.read_csv(
                    result.output_dir / "annotation_decision_counts.tsv", sep="\t"
                )
                rejected_iea = counts[
                    (counts["disposition"] == "rejected")
                    & (counts["reason"] == "evidence_code")
                    & (counts["evidence_code"] == "IEA")
                ]
                self.assertEqual(int(rejected_iea["rows"].sum()), 26)
                self.assertFalse((result.output_dir / "uniref90_clusters.tsv").exists())
                self.assertFalse((result.output_dir / "uniref90_to_mmseqs_cluster.tsv").exists())
            self.assertEqual(
                hashlib.sha256(
                    (outputs[0] / "excluded_annotations_sample.tsv.gz").read_bytes()
                ).hexdigest(),
                hashlib.sha256(
                    (outputs[1] / "excluded_annotations_sample.tsv.gz").read_bytes()
                ).hexdigest(),
            )

    def test_seed_and_min_count_scopes_coexist_but_exact_scope_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = build_benchmark(fixture_config(root / "outputs", root / "temp"))
            second = build_benchmark(fixture_config(
                root / "outputs", root / "temp", seed=1
            ))
            third = build_benchmark(fixture_config(
                root / "outputs", root / "temp", min_count=2
            ))
            self.assertEqual(len({first.output_dir, second.output_dir, third.output_dir}), 3)
            for result in (first, second, third):
                validate_publication(result.output_dir)
            with self.assertRaises(FileExistsError):
                build_benchmark(fixture_config(root / "outputs", root / "temp"))

    def test_unsupported_population_fails_before_repository_or_input_work(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = fixture_config(
                root / "outputs", root / "temp", training_population="all-cluster-members"
            )
            with mock.patch("homology_cluster_benchmark.pipeline.git_state") as git_probe:
                with self.assertRaisesRegex(ValueError, "intentionally unsupported"):
                    build_benchmark(config)
            git_probe.assert_not_called()

    def test_all_six_fixture_thresholds_publish_valid_scoped_aggregate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for identity in SUPPORTED_IDENTITIES:
                build_benchmark(fixture_config(
                    root / f"timestamped-{int(identity * 100)}", root / "temp", identity=identity
                ))
            aggregate = _cross_threshold_reports(
                [root / f"timestamped-{int(identity * 100)}" for identity in SUPPORTED_IDENTITIES],
                root / "aggregate",
            )
            validate_publication(aggregate)
            validation = json.loads((aggregate / "validation_report.json").read_text())
            marker = json.loads((aggregate / "RUN_COMPLETE.json").read_text())
            self.assertTrue(validation["valid"])
            self.assertEqual(
                json.loads((aggregate / "all_thresholds_summary.json").read_text())["identities"],
                [30, 25, 20, 15, 10, 5],
            )
            self.assertEqual(marker["benchmark_scope"], "all-thresholds-summary")
            self.assertFalse(marker["production_eligible"])
            self.assertFalse(list(aggregate.glob("*-training.csv")))
            roots = [
                root / f"timestamped-{int(identity * 100)}" for identity in SUPPORTED_IDENTITIES
            ]
            with self.assertRaisesRegex(ValueError, "exactly six"):
                _cross_threshold_reports(roots[:5], root / "missing")
            with self.assertRaisesRegex(ValueError, "distinct"):
                _cross_threshold_reports([*roots[:5], roots[0]], root / "duplicate")
            different = build_benchmark(fixture_config(
                root / "timestamped-different-5", root / "temp", identity=0.05, seed=1
            ))
            with self.assertRaisesRegex(ValueError, "scientific fingerprint"):
                _cross_threshold_reports([*roots[:5], different.output_dir], root / "mismatch")

            base = fixture_config(root / "unused", root / "temp", identity=0.05)
            trembl = build_benchmark(replace(
                base,
                output_dir=root / "timestamped-trembl-5",
                uniprot_source_scope="trembl-only",
                uniprot_sprot_sequences=None,
                uniprot_trembl_sequences=replace(
                    base.uniprot_sprot_sequences,
                    name="uniprot_trembl_sequences",
                    source_population="trembl",
                ),
            ))
            with self.assertRaisesRegex(ValueError, "source scope|do not share"):
                _cross_threshold_reports(
                    [*roots[:5], trembl.output_dir], root / "mixed-source"
                )

            changed_manifest = build_benchmark(replace(
                base,
                output_dir=root / "timestamped-manifest-5",
                uniref90_fasta=replace(
                    base.uniref90_fasta,
                    url="https://synthetic.invalid/alternate-uniref90.fasta",
                ),
            ))
            with self.assertRaisesRegex(ValueError, "frozen manifest|do not share"):
                _cross_threshold_reports(
                    [*roots[:5], changed_manifest.output_dir], root / "mixed-manifest"
                )

            changed_commit = build_benchmark(replace(
                base, output_dir=root / "timestamped-commit-5"
            ))
            _rewrite_fixture_repository_commit(changed_commit.output_dir, "b" * 40)
            validate_publication(changed_commit.output_dir)
            with self.assertRaisesRegex(ValueError, "framework revision|do not share"):
                _cross_threshold_reports(
                    [*roots[:5], changed_commit.output_dir], root / "mixed-commit"
                )

    def test_immutable_pfp_preparation_consumes_all_nine_csvs(self):
        pfp_root = Path(
            os.environ.get("PFP_REPOSITORY", str(Path.home() / "temporary" / "PFP"))
        )
        preparer = pfp_root / "scripts" / "prepare_cafa3_data.py"
        if not preparer.is_file():
            self.skipTest("immutable PFP checkout is unavailable for ingestion integration")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = build_benchmark(fixture_config(root / "outputs", root / "temp"))
            processed = root / "processed"
            completed = subprocess.run(
                [
                    sys.executable, str(preparer), "--cafa3-dir", str(result.output_dir),
                    "--output-dir", str(processed),
                ],
                cwd=pfp_root,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout)
            for ontology, root_term in (
                ("BPO", "GO:0008150"), ("CCO", "GO:0005575"),
                ("MFO", "GO:0003674"),
            ):
                terms = json.loads((processed / f"{ontology}_go_terms.json").read_text())
                self.assertIn(root_term, terms)
                for split in ("train", "valid", "test"):
                    self.assertTrue((processed / f"{ontology}_{split}_labels.npz").is_file())

    def test_all_cluster_members_is_guarded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = fixture_config(
                root / "outputs", root / "temp", training_population="all-cluster-members"
            )
            with self.assertRaisesRegex(ValueError, "intentionally unsupported"):
                build_benchmark(config)


if __name__ == "__main__":
    unittest.main()
