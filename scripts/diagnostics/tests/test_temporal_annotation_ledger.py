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
sys.path.insert(0, str(DIAGNOSTICS))

from temporal_annotation_common import (  # noqa: E402
    build_temporal_state_rows,
    descriptive_cohort_masks,
)


BUILDER = DIAGNOSTICS / "build_temporal_annotation_ledger.py"


def write_tsv(path: Path, fields: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


class TemporalAnnotationLedgerTests(unittest.TestCase):
    def test_direct_states_and_closure_transitions_are_kept_separate(self) -> None:
        direct0 = {
            "P1": {
                "BPO": frozenset(),
                "CCO": frozenset({"GO:0000002"}),
                "MFO": frozenset(),
            },
            "P2": {
                "BPO": frozenset({"GO:0000001"}),
                "CCO": frozenset(),
                "MFO": frozenset(),
            },
        }
        direct1 = {
            "P1": {
                "BPO": frozenset({"GO:0000003"}),
                "CCO": frozenset({"GO:0000002"}),
                "MFO": frozenset(),
            },
            "P2": {
                "BPO": frozenset({"GO:0000004"}),
                "CCO": frozenset(),
                "MFO": frozenset(),
            },
        }
        closure0 = direct0
        closure1 = {
            **direct1,
            "P2": {
                "BPO": frozenset({"GO:0000001", "GO:0000004"}),
                "CCO": frozenset(),
                "MFO": frozenset(),
            },
        }
        rows = build_temporal_state_rows(
            ["P1", "P2", "P3"],
            direct0,
            direct1,
            closure0,
            closure1,
            {"P1", "P2"},
            {"P1", "P2", "P3"},
        )
        masks = descriptive_cohort_masks(rows, ["P1", "P2"], "BPO")
        self.assertEqual(masks["global_t0_empty"].tolist(), [False, False])
        self.assertEqual(masks["aspect_t0_empty"].tolist(), [True, False])
        self.assertEqual(masks["aspect_has_gain"].tolist(), [True, True])
        self.assertEqual(masks["cross_ontology_known"].tolist(), [True, False])
        self.assertEqual(masks["same_aspect_partial"].tolist(), [False, True])
        p2_bp = next(
            row for row in rows if row["protein_id"] == "P2" and row["aspect"] == "BPO"
        )
        self.assertEqual(p2_bp["retained_terms"], ("GO:0000001",))
        self.assertEqual(p2_bp["gained_terms"], ("GO:0000004",))
        p3_bp = next(
            row for row in rows if row["protein_id"] == "P3" and row["aspect"] == "BPO"
        )
        self.assertEqual(p3_bp["global_knowledge_state"], "unknown")

    def test_cli_publishes_hash_bound_policy_labelled_states(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            scope = root / "scope.tsv"
            t0_presence = root / "t0_presence.tsv"
            t1_presence = root / "t1_presence.tsv"
            t0_direct = root / "t0_direct.tsv"
            t1_direct = root / "t1_direct.tsv"
            t0_closure = root / "t0_closure.tsv"
            t1_closure = root / "t1_closure.tsv"
            exposure = root / "exposure.tsv"
            output = root / "ledger"
            write_tsv(scope, ["protein_id"], [["P1"], ["P2"], ["P3"]])
            write_tsv(t0_presence, ["protein_id"], [["P1"], ["P2"]])
            write_tsv(t1_presence, ["protein_id"], [["P1"], ["P2"], ["P3"]])
            write_tsv(
                t0_direct,
                ["protein_id", "aspect", "go_term"],
                [["P1", "CCO", "GO:0000002"], ["P2", "BPO", "GO:0000001"]],
            )
            write_tsv(
                t1_direct,
                ["protein_id", "aspect", "go_term"],
                [
                    ["P1", "CCO", "GO:0000002"],
                    ["P1", "BPO", "GO:0000003"],
                    ["P2", "BPO", "GO:0000004"],
                    ["P3", "BPO", "GO:0000003"],
                ],
            )
            write_tsv(
                t0_closure,
                ["protein_id", "aspect", "go_term"],
                [["P1", "CCO", "GO:0000002"], ["P2", "BPO", "GO:0000001"]],
            )
            write_tsv(
                t1_closure,
                ["protein_id", "aspect", "go_term"],
                [
                    ["P1", "CCO", "GO:0000002"],
                    ["P1", "BPO", "GO:0000003"],
                    ["P2", "BPO", "GO:0000001"],
                    ["P2", "BPO", "GO:0000004"],
                    ["P3", "BPO", "GO:0000003"],
                ],
            )
            write_tsv(
                exposure,
                [
                    "protein_id",
                    "train_id_member",
                    "valid_id_member",
                    "train_sequence_member",
                    "valid_sequence_member",
                    "train_homology_cluster_member",
                    "modality_availability",
                    "feature_temporal_policy",
                ],
                [
                    ["P1", "1", "0", "1", "0", "unknown", "sequence,text", "t0_frozen"],
                    ["P2", "0", "1", "0", "1", "unknown", "sequence,text", "t0_frozen"],
                ],
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(BUILDER),
                    "--t0-direct-annotations",
                    str(t0_direct),
                    "--t1-direct-annotations",
                    str(t1_direct),
                    "--t0-closure-annotations",
                    str(t0_closure),
                    "--t1-closure-annotations",
                    str(t1_closure),
                    "--t0-protein-presence",
                    str(t0_presence),
                    "--t1-protein-presence",
                    str(t1_presence),
                    "--exposure-table",
                    str(exposure),
                    "--protein-scope",
                    str(scope),
                    "--output-dir",
                    str(output),
                    "--t0-snapshot",
                    "2025-03-08",
                    "--t1-snapshot",
                    "2026-06-17",
                    "--evidence-policy-id",
                    "supervisor_snapshot_membership",
                    "--graph-policy-id",
                    "cafa_narrow_is_a_part_of",
                    "--relationship",
                    "is_a",
                    "--relationship",
                    "part_of",
                    "--benchmark-id",
                    "fixture",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(
                (output / "temporal_annotation_ledger.json").read_text()
            )
            self.assertEqual(report["status"], "complete")
            self.assertEqual(report["t0_snapshot"], "2025-03-08")
            self.assertEqual(report["t1_snapshot"], "2026-06-17")
            self.assertEqual(report["protein_count"], 3)
            self.assertEqual(report["cohort_counts"]["BPO"]["cross_ontology_known"], 1)
            self.assertEqual(
                report["gainer_cohort_counts"]["BPO"]["cross_ontology_known"],
                1,
            )
            self.assertEqual(
                report["gainer_cohort_counts"]["BPO"]["same_aspect_partial"],
                1,
            )
            self.assertEqual(report["global_knowledge_counts"]["unknown"], 1)
            self.assertEqual(report["transition_counts"]["retained_known"], 2)
            self.assertTrue((output / "output_manifest.json").is_file())
            self.assertTrue((output / "RUN_COMPLETE.json").is_file())
            with (output / "protein_cohorts.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            p1_bp = next(
                row
                for row in rows
                if row["protein_id"] == "P1" and row["aspect"] == "BPO"
            )
            self.assertEqual(p1_bp["global_t0_empty"], "0")
            self.assertEqual(p1_bp["aspect_t0_empty"], "1")
            self.assertEqual(p1_bp["aspect_knowledge_state"], "cross_ontology_known")
            self.assertEqual(p1_bp["gained_terms"], "GO:0000003")
            p2_bp = next(
                row
                for row in rows
                if row["protein_id"] == "P2" and row["aspect"] == "BPO"
            )
            self.assertEqual(p2_bp["retained_terms"], "GO:0000001")

    def test_duplicate_annotation_row_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            scope = root / "scope.tsv"
            presence = root / "presence.tsv"
            t0_direct = root / "t0_direct.tsv"
            t1_direct = root / "t1_direct.tsv"
            t0_closure = root / "t0_closure.tsv"
            t1_closure = root / "t1_closure.tsv"
            write_tsv(scope, ["protein_id"], [["P1"]])
            write_tsv(presence, ["protein_id"], [["P1"]])
            duplicate = [["P1", "BPO", "GO:0000001"]] * 2
            write_tsv(t0_direct, ["protein_id", "aspect", "go_term"], duplicate)
            for path in (t1_direct, t0_closure, t1_closure):
                write_tsv(path, ["protein_id", "aspect", "go_term"], [])
            result = subprocess.run(
                [
                    sys.executable,
                    str(BUILDER),
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
                    str(root / "output"),
                    "--t0-snapshot",
                    "2025-03-08",
                    "--t1-snapshot",
                    "2026-06-17",
                    "--evidence-policy-id",
                    "supervisor",
                    "--graph-policy-id",
                    "narrow",
                    "--relationship",
                    "is_a",
                    "--benchmark-id",
                    "fixture",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Duplicate annotation row", result.stderr)


if __name__ == "__main__":
    unittest.main()
