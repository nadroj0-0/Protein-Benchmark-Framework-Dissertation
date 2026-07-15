from __future__ import annotations

import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pfp_benchmark_reuse.benchmark import parse_benchmark  # noqa: E402
from pfp_benchmark_reuse.models import REGENERATE_MODALITIES  # noqa: E402
from pfp_benchmark_reuse.planner import (  # noqa: E402
    PlanningError,
    build_plan,
    validate_plan,
)

from helpers import rows_in, write_benchmark  # noqa: E402


class PlannerPolicyTests(unittest.TestCase):
    def test_self_comparison_is_all_reuse(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = write_benchmark(
                Path(tmp) / "benchmark",
                rows_in("bp-training.csv", ("P2", "BBBB"), ("P1", "AAAA")),
            )
            source = parse_benchmark("source", directory)
            target = parse_benchmark("target", directory)
            plan = build_plan((source,), target)
            self.assertEqual([record.protein_id for record in plan.reuse_records], ["P1", "P2"])
            self.assertEqual(plan.regenerate_records, ())

    def test_adding_one_new_protein_gives_one_regenerate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark("source", write_benchmark(root / "source"))
            target = parse_benchmark(
                "target",
                write_benchmark(
                    root / "target",
                    rows_in("bp-training.csv", ("P1", "AAAA"), ("P2", "BBBB")),
                ),
            )
            plan = build_plan((source,), target)
            self.assertEqual([record.protein_id for record in plan.regenerate_records], ["P2"])
            self.assertEqual(plan.regenerate_records[0].reason, "protein-id-absent")
            self.assertEqual(
                plan.regenerate_records[0].regenerate_modalities, REGENERATE_MODALITIES
            )

    def test_multiple_references_contribute_independently_and_aggregate_matches(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = parse_benchmark(
                "alpha",
                write_benchmark(
                    root / "first",
                    rows_in("bp-training.csv", ("P1", "AAAA"), ("P3", "CCCC")),
                ),
            )
            second = parse_benchmark(
                "zeta",
                write_benchmark(
                    root / "second",
                    rows_in("mf-test.csv", ("P2", "BBBB"), ("P3", "CCCC")),
                ),
            )
            target = parse_benchmark(
                "target",
                write_benchmark(
                    root / "target",
                    rows_in(
                        "cc-validation.csv",
                        ("P1", "AAAA"),
                        ("P2", "BBBB"),
                        ("P3", "CCCC"),
                    ),
                ),
            )
            plan = build_plan((second, first), target)
            self.assertEqual([record.protein_id for record in plan.reuse_records], ["P1", "P2", "P3"])
            p2 = next(record for record in plan.records if record.protein_id == "P2")
            p3 = next(record for record in plan.records if record.protein_id == "P3")
            self.assertEqual(p2.matching_embedded_benchmarks, ("zeta",))
            self.assertEqual(p3.matching_embedded_benchmarks, ("alpha", "zeta"))

    def test_changed_sequence_and_different_id_regenerate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark("source", write_benchmark(root / "source"))
            changed = parse_benchmark(
                "changed",
                write_benchmark(root / "changed", rows_in("bp-training.csv", ("P1", "BBBB"))),
            )
            other_id = parse_benchmark(
                "other",
                write_benchmark(root / "other", rows_in("bp-training.csv", ("P9", "AAAA"))),
            )
            changed_record = build_plan((source,), changed).regenerate_records[0]
            other_record = build_plan((source,), other_id).regenerate_records[0]
            self.assertEqual(changed_record.reason, "sequence-mismatch")
            self.assertEqual(other_record.reason, "protein-id-absent")

    def test_protein_ids_and_sequences_are_case_sensitive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark("source", write_benchmark(root / "source"))
            lowercase_id = parse_benchmark(
                "lowercase-id",
                write_benchmark(root / "lowercase-id", rows_in("bp-training.csv", ("p1", "AAAA"))),
            )
            lowercase_sequence = parse_benchmark(
                "lowercase-sequence",
                write_benchmark(
                    root / "lowercase-sequence",
                    rows_in("bp-training.csv", ("P1", "aaaa")),
                ),
            )
            self.assertEqual(
                build_plan((source,), lowercase_id).regenerate_records[0].reason,
                "protein-id-absent",
            )
            self.assertEqual(
                build_plan((source,), lowercase_sequence).regenerate_records[0].reason,
                "sequence-mismatch",
            )

    def test_split_ontology_and_label_changes_individually_still_reuse(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark(
                "source",
                write_benchmark(
                    root / "source", rows_in("bp-training.csv", ("P1", "AAAA", "1"))
                ),
            )
            variants = {
                "split-only": ("bp-test.csv", "1"),
                "ontology-only": ("mf-training.csv", "1"),
                "label-only": ("bp-training.csv", "0"),
            }
            for index, (case, (membership, label)) in enumerate(variants.items()):
                with self.subTest(case=case):
                    target = parse_benchmark(
                        "target",
                        write_benchmark(
                            root / ("target-%d" % index),
                            rows_in(membership, ("P1", "AAAA", label)),
                        ),
                    )
                    record = build_plan((source,), target).reuse_records[0]
                    self.assertEqual(record.target_memberships, (membership,))
                    self.assertEqual(record.matching_embedded_benchmarks, ("source",))

    def test_reference_conflict_fails_even_when_absent_from_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = parse_benchmark(
                "first",
                write_benchmark(root / "first", rows_in("bp-training.csv", ("CONFLICT", "AAAA"))),
            )
            second = parse_benchmark(
                "second",
                write_benchmark(root / "second", rows_in("bp-training.csv", ("CONFLICT", "BBBB"))),
            )
            target = parse_benchmark(
                "target",
                write_benchmark(root / "target", rows_in("bp-training.csv", ("P1", "CCCC"))),
            )
            with self.assertRaisesRegex(PlanningError, "conflicting sequences"):
                build_plan((first, second), target)

    def test_duplicate_names_and_target_name_collision_fail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = parse_benchmark("same", write_benchmark(root / "first"))
            second = parse_benchmark("same", write_benchmark(root / "second"))
            target = parse_benchmark("target", write_benchmark(root / "target"))
            with self.assertRaisesRegex(PlanningError, "must be unique"):
                build_plan((first, second), target)

            colliding_target = parse_benchmark("same", write_benchmark(root / "collision"))
            with self.assertRaisesRegex(PlanningError, "must be unique"):
                build_plan((first,), colliding_target)

    def test_known_union_includes_reference_only_proteins(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark(
                "source",
                write_benchmark(
                    root / "source",
                    rows_in("bp-training.csv", ("P1", "AAAA"), ("EXTRA", "EEEE")),
                ),
            )
            target = parse_benchmark("target", write_benchmark(root / "target"))
            plan = build_plan((source,), target)
            self.assertEqual(
                [record.protein_id for record in plan.known_embedded_proteins],
                ["EXTRA", "P1"],
            )
            self.assertEqual([record.protein_id for record in plan.records], ["P1"])

    def test_partition_validator_rejects_missing_and_third_action(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = parse_benchmark("source", write_benchmark(root / "source"))
            target = parse_benchmark("target", write_benchmark(root / "target"))
            plan = build_plan((source,), target)
            with self.assertRaisesRegex(PlanningError, "does not equal"):
                validate_plan(replace(plan, records=()))
            with self.assertRaisesRegex(PlanningError, "duplicate target protein IDs"):
                validate_plan(replace(plan, records=(plan.records[0], plan.records[0])))
            bad_record = replace(plan.records[0], action="third")
            with self.assertRaisesRegex(PlanningError, "reuse or regenerate"):
                validate_plan(replace(plan, records=(bad_record,)))
            forged = replace(
                plan.records[0],
                sequence="ZZZZ",
                sequence_sha256="0" * 64,
                target_memberships=("mf-test.csv",),
            )
            with self.assertRaisesRegex(PlanningError, "do not exactly match target"):
                validate_plan(replace(plan, records=(forged,)))


if __name__ == "__main__":
    unittest.main()
