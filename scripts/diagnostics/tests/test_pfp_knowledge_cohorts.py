from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


FRAMEWORK = Path(__file__).parents[3]
DIAGNOSTICS = FRAMEWORK / "scripts" / "diagnostics"
MODEL_EXECUTION = FRAMEWORK / "scripts" / "model_execution"
sys.path.insert(0, str(MODEL_EXECUTION))
sys.path.insert(0, str(DIAGNOSTICS))

import test_pfp_label_sensitivity as sensitivity_tests  # noqa: E402
from test_temporal_annotation_ledger import write_tsv  # noqa: E402


BUILD_LEDGER = DIAGNOSTICS / "build_temporal_annotation_ledger.py"
EVALUATE = DIAGNOSTICS / "evaluate_pfp_knowledge_cohorts.py"


class PfpKnowledgeCohortTests(unittest.TestCase):
    def test_accepted_global_no_knowledge_has_no_known_comparator(self) -> None:
        helper = sensitivity_tests.PfpLabelSensitivityTests(methodName="runTest")
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            obo = root / "go.obo"
            helper.make_obo(obo)
            prediction = helper.make_prediction_artifact(root, obo)
            scope = root / "scope.tsv"
            presence = root / "presence.tsv"
            t0_direct = root / "t0_direct.tsv"
            t1_direct = root / "t1_direct.tsv"
            t0_closure = root / "t0_closure.tsv"
            t1_closure = root / "t1_closure.tsv"
            proteins = [["ROOT_ONLY"], ["DEEP"], ["ALL_ZERO"]]
            write_tsv(scope, ["protein_id"], proteins)
            write_tsv(presence, ["protein_id"], proteins)
            write_tsv(t0_direct, ["protein_id", "aspect", "go_term"], [])
            write_tsv(t0_closure, ["protein_id", "aspect", "go_term"], [])
            write_tsv(
                t1_direct,
                ["protein_id", "aspect", "go_term"],
                [
                    ["ROOT_ONLY", "BPO", "GO:0008150"],
                    ["DEEP", "BPO", "GO:0009987"],
                ],
            )
            write_tsv(
                t1_closure,
                ["protein_id", "aspect", "go_term"],
                [
                    ["ROOT_ONLY", "BPO", "GO:0008150"],
                    ["DEEP", "BPO", "GO:0008150"],
                    ["DEEP", "BPO", "GO:0009987"],
                ],
            )
            ledger = root / "ledger"
            build = subprocess.run(
                [
                    sys.executable,
                    str(BUILD_LEDGER),
                    "--t0-direct-annotations",
                    str(t0_direct),
                    "--t1-direct-annotations",
                    str(t1_direct),
                    "--t0-closure-annotations",
                    str(t0_closure),
                    "--t1-closure-annotations",
                    str(t1_closure),
                    "--t0-protein-presence",
                    str(presence),
                    "--t1-protein-presence",
                    str(presence),
                    "--protein-scope",
                    str(scope),
                    "--output-dir",
                    str(ledger),
                    "--t0-snapshot",
                    "t0",
                    "--t1-snapshot",
                    "t1",
                    "--evidence-policy-id",
                    "supervisor",
                    "--graph-policy-id",
                    "cafa_narrow_is_a_part_of",
                    "--relationship",
                    "is_a",
                    "--relationship",
                    "part_of",
                    "--expected-global-knowledge-state",
                    "no_qualifying",
                    "--benchmark-id",
                    "fixture",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(build.returncode, 0, build.stderr)

            output = root / "knowledge"
            result = subprocess.run(
                [
                    sys.executable,
                    str(EVALUATE),
                    "--prediction-manifest",
                    str(prediction),
                    "--temporal-ledger-dir",
                    str(ledger),
                    "--truth-graph-policy-id",
                    "cafa_narrow_is_a_part_of",
                    "--bootstrap-replicates",
                    "5",
                    "--output-dir",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(
                (output / "knowledge_cohort_analysis.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                report["scientific_label"],
                "accepted_global_no_knowledge_cohort_inventory",
            )
            self.assertEqual(
                report["aspects"]["BPO"]["cohort_counts"]["global_no_qualifying"],
                3,
            )
            with (output / "cohort_metrics.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            known = next(
                row
                for row in rows
                if row["aspect"] == "BPO"
                and row["cohort"] == "global_known_qualifying"
                and row["target_component"] == "combined_t1_truth"
            )
            self.assertEqual(known["status"], "not_evaluable_empty_cohort")


if __name__ == "__main__":
    unittest.main()
