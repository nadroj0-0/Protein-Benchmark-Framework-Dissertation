from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from homology_cluster_benchmark.attrition import (
    METRIC_DEFINITIONS,
    evaluate_attrition,
    load_attrition_policy,
    observation,
)
from homology_cluster_benchmark.inputs import sha256_file


COMMIT = "a" * 40
MANIFEST_HASH = "b" * 64
INPUT_HASH = "c" * 64
RELEASES = {
    "uniprot_uniref": "2026_02",
    "goa": "234",
    "ontology": "releases/2026-06-15",
}


def _policy_payload(bound: float = 0.5) -> dict:
    metrics = {}
    for name, definition in METRIC_DEFINITIONS.items():
        metrics[name] = {
            "numerator_definition": definition.numerator,
            "denominator_definition": definition.denominator,
            f"allowed_{definition.bound}_ratio": bound,
            "rationale": "Reviewed fixture boundary.",
            "evidence_source": "reviewed synthetic 30 percent pilot",
        }
    return {
        "schema_name": "homology-cluster-attrition-policy",
        "schema_version": 1,
        "uniprot_source_scope": "sprot-only",
        "expected_releases": RELEASES,
        "metrics": metrics,
        "rationale": "Exercise exact policy boundaries.",
        "evidence_source": "reviewed synthetic 30 percent pilot",
        "author": "fixture-author",
        "reviewer": "fixture-reviewer",
        "review_date": "2026-07-14",
        "framework_commit": COMMIT,
        "frozen_input_manifest_sha256": MANIFEST_HASH,
    }


def _observations(value: float) -> dict:
    return {
        name: observation(name, value, 1.0)
        for name in METRIC_DEFINITIONS
    }


class AttritionPolicyTests(unittest.TestCase):
    def test_policy_schema_and_bindings_are_strict(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "policy.json"
            path.write_text(json.dumps(_policy_payload()), encoding="utf-8")
            payload, digest = load_attrition_policy(
                path,
                source_scope="sprot-only",
                expected_releases=RELEASES,
                framework_commit=COMMIT,
                frozen_input_manifest_sha256=MANIFEST_HASH,
            )
            self.assertEqual(digest, sha256_file(path))
            self.assertEqual(set(payload["metrics"]), set(METRIC_DEFINITIONS))
            for label, kwargs, message in (
                ("scope", {"source_scope": "trembl-only"}, "wrong UniProt"),
                ("commit", {"framework_commit": "d" * 40}, "wrong framework"),
                ("manifest", {"frozen_input_manifest_sha256": "e" * 64}, "wrong frozen"),
                ("release", {"expected_releases": {**RELEASES, "goa": "999"}}, "wrong frozen releases"),
            ):
                with self.subTest(label=label), self.assertRaisesRegex(ValueError, message):
                    load_attrition_policy(
                        path,
                        source_scope=kwargs.get("source_scope", "sprot-only"),
                        expected_releases=kwargs.get("expected_releases", RELEASES),
                        framework_commit=kwargs.get("framework_commit", COMMIT),
                        frozen_input_manifest_sha256=kwargs.get(
                            "frozen_input_manifest_sha256", MANIFEST_HASH
                        ),
                    )
            malformed = root / "malformed.json"
            malformed.write_text("[]", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "root must be an object"):
                load_attrition_policy(
                    malformed,
                    source_scope="sprot-only",
                    expected_releases=RELEASES,
                    framework_commit=COMMIT,
                    frozen_input_manifest_sha256=MANIFEST_HASH,
                )
            placeholder = _policy_payload()
            placeholder["reviewer"] = "REPLACE_ME"
            path.write_text(json.dumps(placeholder), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "template placeholder"):
                load_attrition_policy(
                    path,
                    source_scope="sprot-only",
                    expected_releases=RELEASES,
                    framework_commit=COMMIT,
                    frozen_input_manifest_sha256=MANIFEST_HASH,
                )

    def test_every_metric_accepts_exact_boundary_and_fails_just_beyond(self):
        policy = _policy_payload(0.5)
        exact = _observations(0.5)
        report = evaluate_attrition(
            policy, "f" * 64, exact,
            source_scope="sprot-only", framework_commit=COMMIT,
            input_manifest_sha256=INPUT_HASH,
        )
        self.assertTrue(report["policy_passed"])
        diagnostic = evaluate_attrition(
            policy, "f" * 64, exact,
            source_scope="sprot-only", framework_commit=COMMIT,
            input_manifest_sha256=INPUT_HASH, diagnostic=True,
        )
        self.assertTrue(diagnostic["policy_passed"])
        self.assertFalse(diagnostic["production_authorized"])
        for name, definition in METRIC_DEFINITIONS.items():
            with self.subTest(metric=name):
                changed = _observations(0.5)
                outside = 0.499999 if definition.bound == "minimum" else 0.500001
                changed[name] = observation(name, outside, 1.0)
                failed = evaluate_attrition(
                    policy, "f" * 64, changed,
                    source_scope="sprot-only", framework_commit=COMMIT,
                    input_manifest_sha256=INPUT_HASH,
                    diagnostic=True,
                )
                self.assertFalse(failed["policy_passed"])
                self.assertEqual([item["name"] for item in failed["failed_metrics"]], [name])
                self.assertTrue(failed["diagnostic"])

    def test_reviewed_override_must_exactly_bind_failures_and_input_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            policy = _policy_payload(0.5)
            observed = _observations(0.5)
            observed["goa_to_selected_uniprot_mapping_ratio"] = observation(
                "goa_to_selected_uniprot_mapping_ratio", 0.4, 1
            )
            override = {
                "schema_name": "homology-cluster-attrition-override",
                "schema_version": 1,
                "failed_metrics": [{
                    "name": "goa_to_selected_uniprot_mapping_ratio",
                    "observed_ratio": 0.4,
                }],
                "justification": "Reviewed fixture exception.",
                "reviewer": "fixture-reviewer",
                "review_date": "2026-07-14",
                "input_manifest_sha256": INPUT_HASH,
                "framework_commit": COMMIT,
                "uniprot_source_scope": "sprot-only",
                "pilot_or_run_identifier": "fixture-pilot-1",
            }
            path = root / "override.json"
            path.write_text(json.dumps(override), encoding="utf-8")
            report = evaluate_attrition(
                policy, "f" * 64, observed,
                source_scope="sprot-only", framework_commit=COMMIT,
                input_manifest_sha256=INPUT_HASH, override_path=path,
            )
            self.assertTrue(report["production_authorized"])
            self.assertTrue(report["override_valid"])
            override["justification"] = "REPLACE_WITH_REVIEWED_JUSTIFICATION"
            path.write_text(json.dumps(override), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "template placeholder"):
                evaluate_attrition(
                    policy, "f" * 64, observed,
                    source_scope="sprot-only", framework_commit=COMMIT,
                    input_manifest_sha256=INPUT_HASH, override_path=path,
                )
            override["justification"] = "Reviewed fixture exception."
            override["input_manifest_sha256"] = "0" * 64
            path.write_text(json.dumps(override), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "wrong input-manifest"):
                evaluate_attrition(
                    policy, "f" * 64, observed,
                    source_scope="sprot-only", framework_commit=COMMIT,
                    input_manifest_sha256=INPUT_HASH, override_path=path,
                )
            path.write_text("[]", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "root must be an object"):
                evaluate_attrition(
                    policy, "f" * 64, observed,
                    source_scope="sprot-only", framework_commit=COMMIT,
                    input_manifest_sha256=INPUT_HASH, override_path=path,
                )


if __name__ == "__main__":
    unittest.main()
