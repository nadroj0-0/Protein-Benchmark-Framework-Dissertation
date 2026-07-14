from __future__ import annotations

from pathlib import Path
from dataclasses import replace
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.models import ClusterInfo
from homology_cluster_benchmark.splitting import assign_splits
from homology_cluster_benchmark.validation import (
    production_balance_within_tolerance,
    split_balance_metrics,
)
from tests.helpers import fixture_config


class SplittingTests(unittest.TestCase):
    def setUp(self):
        self.clusters = {
            f"C{index:02d}": ClusterInfo(f"C{index:02d}", weight, index % 4 + 1)
            for index, weight in enumerate([31, 17, 13, 11, 9, 8, 7, 6, 5, 4, 3, 2], start=1)
        }

    def test_both_policies_are_deterministic_and_cluster_disjoint(self):
        for policy in ("cluster-count-random", "sequence-balanced"):
            first = assign_splits(self.clusters, policy, seed=0)
            second = assign_splits(self.clusters, policy, seed=0)
            self.assertEqual(
                {key: value.split for key, value in first.items()},
                {key: value.split for key, value in second.items()},
            )
            self.assertEqual(set(first), set(self.clusters))
            self.assertEqual({item.split for item in first.values()}, {"training", "validation", "test"})

    def test_cluster_count_random_targets_cluster_counts(self):
        result = assign_splits(self.clusters, "cluster-count-random", seed=3)
        test_count = sum(item.split == "test" for item in result.values())
        self.assertEqual(test_count, 2)
        self.assertGreaterEqual(sum(item.split == "validation" for item in result.values()), 1)

    def test_sequence_balanced_uses_members_not_labelled_proteins(self):
        result = assign_splits(self.clusters, "sequence-balanced", seed=4)
        total = sum(item.member_count for item in result.values())
        test = sum(item.member_count for item in result.values() if item.split == "test")
        self.assertLess(abs(test / total - 0.20), 0.12)

    def test_giant_cluster_edge_case_remains_whole(self):
        clusters = {
            "GIANT": ClusterInfo("GIANT", 900, 1),
            "A": ClusterInfo("A", 40, 50),
            "B": ClusterInfo("B", 30, 50),
            "C": ClusterInfo("C", 20, 50),
            "D": ClusterInfo("D", 10, 50),
        }
        result = assign_splits(clusters, "sequence-balanced", seed=0)
        self.assertIn(result["GIANT"].split, {"training", "validation", "test"})
        self.assertEqual(len(result), len(clusters))

    def test_adversarial_vector_is_near_both_stage_targets_and_deterministic(self):
        weights = [1, 16, 226, 9, 299, 24, 76, 237, 276]
        clusters = {
            f"C{index}": ClusterInfo(f"C{index}", weight, 1)
            for index, weight in enumerate(weights)
        }
        first = assign_splits(clusters, "sequence-balanced", seed=0)
        second = assign_splits(clusters, "sequence-balanced", seed=0)
        self.assertEqual(first, second)
        total = sum(weights)
        development = sum(
            item.member_count for item in first.values() if item.split != "test"
        )
        training = sum(
            item.member_count for item in first.values() if item.split == "training"
        )
        self.assertLessEqual(abs(development / total - 0.80), 0.01)
        self.assertLessEqual(abs(training / development - 0.90), 0.01)

    def test_each_candidate_is_optimized_before_ranking(self):
        weights = [73, 494, 96, 307, 9, 430, 388, 105]
        clusters = {
            f"C{index}": ClusterInfo(f"C{index}", weight, 1)
            for index, weight in enumerate(weights)
        }
        assignments = assign_splits(clusters, "sequence-balanced", seed=0)
        total = sum(weights)
        development = sum(
            item.member_count
            for item in assignments.values()
            if item.split != "test"
        )
        self.assertLessEqual(abs(development / total - 0.80), 0.01)

    def test_impossible_giant_distribution_is_reported_for_reviewed_policy_gate(self):
        clusters = {
            "GIANT": ClusterInfo("GIANT", 900, 1),
            "A": ClusterInfo("A", 40, 1),
            "B": ClusterInfo("B", 30, 1),
            "C": ClusterInfo("C", 20, 1),
            "D": ClusterInfo("D", 10, 1),
        }
        assignments = assign_splits(clusters, "sequence-balanced", seed=0)
        with tempfile.TemporaryDirectory() as tmp:
            fixture = fixture_config(Path(tmp) / "out", Path(tmp) / "temp")
            production = replace(
                fixture, fixture_mode=False, cluster_assignments=None, min_count=50
            )
            metrics = split_balance_metrics(assignments, production)
            self.assertGreater(
                max(metrics["development_deviation"], metrics["training_deviation"]), 0.05
            )
            self.assertTrue(production_balance_within_tolerance(assignments, production))


if __name__ == "__main__":
    unittest.main()
