from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from homology_cluster_benchmark.attrition import METRIC_DEFINITIONS, observation
from homology_cluster_benchmark.authorization import validate_pilot_approval
from homology_cluster_benchmark.inputs import sha256_file


COMMIT = "a" * 40
MANIFEST_HASH = "b" * 64
REVIEWED_POLICY_HASH = "e" * 64


def _reviewed_policy(minimum: float = 0.4, maximum: float = 0.6) -> dict:
    return {
        "metrics": {
            name: {
                f"allowed_{definition.bound}_ratio": (
                    minimum if definition.bound == "minimum" else maximum
                )
            }
            for name, definition in METRIC_DEFINITIONS.items()
        }
    }


class PilotApprovalTests(unittest.TestCase):
    def _files(self, root: Path) -> tuple[Path, Path, Path, Path, Path]:
        marker = root / "RUN_COMPLETE.json"
        marker.write_text(json.dumps({
            "complete": True,
            "benchmark_scope": "diagnostic-pilot",
            "production_eligible": False,
            "identity_percent": 30,
            "framework_revision": COMMIT,
            "repository_commit": COMMIT,
            "frozen_input_manifest_sha256": MANIFEST_HASH,
            "run_input_manifest_sha256": "c" * 64,
            "uniprot_source_scope": "sprot-only",
            "split_policy": "sequence-balanced",
            "training_population": "annotated-only",
            "observed_mmseqs_version": "15-6f452",
            "run_id": "pilot-run",
        }), encoding="utf-8")
        attrition = root / "attrition_report.json"
        attrition.write_text(json.dumps({
            "schema_name": "homology-cluster-attrition-report",
            "schema_version": 1,
            "diagnostic": True,
            "production_authorized": False,
            "uniprot_source_scope": "sprot-only",
            "framework_commit": COMMIT,
            "input_manifest_sha256": "c" * 64,
            "policy_sha256": "d" * 64,
            "metrics": [
                observation(name, 0.5, 1.0) for name in METRIC_DEFINITIONS
            ],
        }), encoding="utf-8")
        task_context = root / "hpc_task_context.json"
        task_context.write_text(json.dumps({
            "job_id": "12345",
            "sge_task_id": 1,
            "identity_percent": 30,
            "uniprot_source_scope": "sprot-only",
            "run_id": "pilot-run",
            "framework_revision": COMMIT,
            "requested_smp_slots": 8,
            "nslots": 8,
            "mmseqs_threads": 8,
        }), encoding="utf-8")
        measurements = root / "pilot_measurements.json"
        measurements.write_text(json.dumps({
            "schema_name": "homology-cluster-pilot-measurement-evidence",
            "schema_version": 1,
            "pilot_job_id": "12345",
            "pilot_task_id": 1,
            "pilot_identity_percent": 30,
            "run_id": "pilot-run",
            "framework_commit": COMMIT,
            "uniprot_source_scope": "sprot-only",
            "successful_completion_marker_sha256": sha256_file(marker),
            "runtime_seconds": 10.0,
            "peak_memory_bytes": 100,
            "scratch_peak_bytes": 200,
            "output_size_bytes": 50,
            "measurement_sources": {
                "runtime": "reviewed Grid Engine accounting",
                "peak_memory": "reviewed Grid Engine accounting",
                "scratch_peak": "reviewed task scratch monitoring",
                "output_size": "reviewed du measurement",
            },
            "reviewer": "fixture-reviewer",
            "review_date": "2026-07-14",
            "evidence_notes": "Reviewed fixture measurement evidence.",
        }), encoding="utf-8")
        approval = root / "approval.json"
        approval.write_text(json.dumps({
            "schema_name": "homology-cluster-pilot-approval",
            "schema_version": 1,
            "approved": True,
            "pilot_job_id": "12345",
            "pilot_run_id": "pilot-run",
            "pilot_task_id": 1,
            "pilot_identity_percent": 30,
            "successful_completion_marker_sha256": sha256_file(marker),
            "framework_commit": COMMIT,
            "frozen_input_manifest_sha256": MANIFEST_HASH,
            "uniprot_source_scope": "sprot-only",
            "split_policy": "sequence-balanced",
            "training_population": "annotated-only",
            "mmseqs_version": "15-6f452",
            "attrition_report_sha256": sha256_file(attrition),
            "reviewed_attrition_policy_sha256": REVIEWED_POLICY_HASH,
            "pilot_task_context_sha256": sha256_file(task_context),
            "pilot_measurement_evidence_sha256": sha256_file(measurements),
            "runtime_seconds": 10.0,
            "peak_memory_bytes": 100,
            "scratch_peak_bytes": 200,
            "output_size_bytes": 50,
            "validation_outcome": "pass",
            "reviewer": "fixture-reviewer",
            "review_date": "2026-07-14",
            "evidence_notes": "Reviewed fixture output and QC.",
        }), encoding="utf-8")
        return approval, marker, attrition, task_context, measurements

    def _validate(
        self,
        approval: Path,
        marker: Path,
        attrition: Path,
        task_context: Path,
        measurements: Path,
        **changes,
    ):
        values = {
            "completion_marker_path": marker,
            "attrition_report_path": attrition,
            "task_context_path": task_context,
            "measurement_evidence_path": measurements,
            "framework_commit": COMMIT,
            "frozen_input_manifest_sha256": MANIFEST_HASH,
            "source_scope": "sprot-only",
            "split_policy": "sequence-balanced",
            "training_population": "annotated-only",
            "mmseqs_version": "15-6f452",
            "reviewed_attrition_policy": _reviewed_policy(),
            "reviewed_attrition_policy_sha256": REVIEWED_POLICY_HASH,
        }
        values.update(changes)
        return validate_pilot_approval(approval, **values)

    def test_valid_reviewed_pilot_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = self._files(Path(tmp))
            self.assertTrue(self._validate(*files)["approved"])

    def test_approval_rejects_all_common_contract_mismatches(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            approval, marker, attrition, task_context, measurements = self._files(root)
            cases = {
                "wrong framework": {"framework_commit": "c" * 40},
                "wrong frozen": {"frozen_input_manifest_sha256": "d" * 64},
                "wrong source": {"source_scope": "trembl-only"},
                "wrong split": {"split_policy": "cluster-count-random"},
                "wrong MMseqs2": {"mmseqs_version": "14-7e284"},
            }
            for label, changes in cases.items():
                with self.subTest(label=label), self.assertRaises(ValueError):
                    self._validate(
                        approval, marker, attrition, task_context, measurements, **changes
                    )

            payload = json.loads(approval.read_text(encoding="utf-8"))
            for key, value in (
                ("pilot_task_id", 2),
                ("pilot_identity_percent", 25),
                ("validation_outcome", "fail"),
                ("approved", False),
            ):
                with self.subTest(key=key):
                    changed = dict(payload)
                    changed[key] = value
                    approval.write_text(json.dumps(changed), encoding="utf-8")
                    with self.assertRaises(ValueError):
                        self._validate(
                            approval, marker, attrition, task_context, measurements
                        )
            approval.write_text(json.dumps(payload), encoding="utf-8")
            marker.write_text(marker.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "marker_sha256 mismatch"):
                self._validate(approval, marker, attrition, task_context, measurements)

    def test_reviewed_policy_and_run_id_are_bound_to_pilot_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            approval, marker, attrition, task_context, measurements = self._files(root)
            failing_policy = _reviewed_policy(minimum=0.6, maximum=0.4)
            with self.assertRaisesRegex(ValueError, "does not accept"):
                self._validate(
                    approval,
                    marker,
                    attrition,
                    task_context,
                    measurements,
                    reviewed_attrition_policy=failing_policy,
                )

            approval, marker, attrition, task_context, measurements = self._files(root)
            marker_payload = json.loads(marker.read_text(encoding="utf-8"))
            marker_payload["run_id"] = "different-run"
            marker.write_text(json.dumps(marker_payload), encoding="utf-8")
            measurement_payload = json.loads(measurements.read_text(encoding="utf-8"))
            measurement_payload["successful_completion_marker_sha256"] = sha256_file(marker)
            measurements.write_text(json.dumps(measurement_payload), encoding="utf-8")
            approval_payload = json.loads(approval.read_text(encoding="utf-8"))
            approval_payload["successful_completion_marker_sha256"] = sha256_file(marker)
            approval_payload["pilot_measurement_evidence_sha256"] = sha256_file(measurements)
            approval.write_text(json.dumps(approval_payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "required binding run_id"):
                self._validate(approval, marker, attrition, task_context, measurements)

    def test_review_placeholders_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            approval, marker, attrition, task_context, measurements = self._files(root)
            payload = json.loads(approval.read_text(encoding="utf-8"))
            payload["reviewer"] = "REPLACE_ME"
            approval.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "template placeholder"):
                self._validate(approval, marker, attrition, task_context, measurements)

    def test_context_measurements_and_attrition_are_evidence_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            approval, marker, attrition, task_context, measurements = self._files(root)
            context_payload = json.loads(task_context.read_text(encoding="utf-8"))
            context_payload["nslots"] = 7
            task_context.write_text(json.dumps(context_payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "task_context_sha256 mismatch"):
                self._validate(approval, marker, attrition, task_context, measurements)

            approval, marker, attrition, task_context, measurements = self._files(root)
            approval_payload = json.loads(approval.read_text(encoding="utf-8"))
            approval_payload["runtime_seconds"] = 0
            approval.write_text(json.dumps(approval_payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "finite positive"):
                self._validate(approval, marker, attrition, task_context, measurements)

            approval, marker, attrition, task_context, measurements = self._files(root)
            attrition_payload = json.loads(attrition.read_text(encoding="utf-8"))
            attrition_payload["uniprot_source_scope"] = "trembl-only"
            attrition.write_text(json.dumps(attrition_payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "attrition_report_sha256 mismatch"):
                self._validate(approval, marker, attrition, task_context, measurements)


if __name__ == "__main__":
    unittest.main()
