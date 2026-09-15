from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


FRAMEWORK = Path(__file__).parents[3]
DIAGNOSTICS = FRAMEWORK / "scripts" / "diagnostics"
sys.path.insert(0, str(DIAGNOSTICS))
ANALYSIS = DIAGNOSTICS / "analyze_pfp_specificity_support.py"


def load_analysis():
    specification = importlib.util.spec_from_file_location(
        "test_specificity_support_analysis", ANALYSIS
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot import {ANALYSIS}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def row(group: str, support: int, test: int = 5, f1: float = 0.0):
    return {
        "root": 0,
        "xu_bin": group,
        "training_positive_proteins": support,
        "test_positive_proteins": test,
        "term_f1": f1,
    }


class SpecificitySupportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module = load_analysis()

    def test_term_metrics_include_scores_equal_to_threshold(self) -> None:
        truth = np.asarray([[1, 0], [0, 1]], dtype=np.uint8)
        scores = np.asarray([[0.5, 0.5], [0.1, 0.4]], dtype=np.float64)
        result = self.module.term_metrics(truth, scores, 0.5)
        self.assertEqual(result["tp"].tolist(), [1, 0])
        self.assertEqual(result["fp"].tolist(), [0, 1])
        self.assertEqual(result["fn"].tolist(), [0, 1])
        self.assertAlmostEqual(result["f1"][0], 1.0)
        self.assertAlmostEqual(result["f1"][1], 0.0)

    def test_common_support_contract_is_tie_preserving_and_adequate(self) -> None:
        rows = []
        for support in (10, 20, 30, 40):
            rows.extend(row("specificity_q1", support) for _ in range(3))
            rows.extend(row("specificity_q4", support) for _ in range(3))
        contract = self.module.build_support_contract(rows, 4, 3)
        self.assertEqual(contract["status"], "complete")
        self.assertEqual(contract["common_support_min"], 10)
        self.assertEqual(contract["common_support_max"], 40)
        self.assertTrue(contract["bins"])
        for specification in contract["bins"]:
            self.assertGreaterEqual(specification["q1_term_count"], 3)
            self.assertGreaterEqual(specification["q4_term_count"], 3)

    def test_sparse_adjacent_bins_are_merged_deterministically(self) -> None:
        rows = []
        for support in range(1, 9):
            rows.append(row("specificity_q1", support))
            rows.append(row("specificity_q4", support))
        contract = self.module.build_support_contract(rows, 4, 4)
        self.assertEqual(contract["status"], "complete")
        self.assertLessEqual(len(contract["bins"]), 2)
        self.assertEqual(contract["bins"][0]["lower_inclusive"], 1)
        self.assertEqual(contract["bins"][-1]["upper_inclusive"], 8)

    def test_no_common_support_returns_explicit_failure_state(self) -> None:
        rows = [row("specificity_q1", value) for value in range(1, 6)]
        rows.extend(row("specificity_q4", value) for value in range(10, 15))
        contract = self.module.build_support_contract(rows, 4, 2)
        self.assertEqual(contract["status"], "insufficient_overlap")
        self.assertIn("do not overlap", contract["reason"])

    def test_sensitivity_filter_reuses_fixed_band_and_reports_sparse_cell(self) -> None:
        rows = [row("specificity_q1", 10, test=5, f1=0.8) for _ in range(3)]
        rows.extend(row("specificity_q4", 10, test=1, f1=0.3) for _ in range(3))
        summary = self.module.summarize_q1_q4(
            rows, "BPO", 5, "support_bin_1", 10, 10, 3
        )
        self.assertEqual(summary["q1_broad_term_count"], 3)
        self.assertEqual(summary["q4_specific_term_count"], 0)
        self.assertFalse(summary["meets_original_minimum_cell_rule"])

        report = {
            "aspects": {
                aspect: {"support_contract": {"status": "complete"}}
                for aspect in self.module.ASPECTS
            }
        }
        rendered = self.module.markdown_report(report, [], [summary])
        self.assertIn("| NA |", rendered)


if __name__ == "__main__":
    unittest.main()
