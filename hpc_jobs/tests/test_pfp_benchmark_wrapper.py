#!/usr/bin/env python3
"""Static regression tests for the generic PFP Grid Engine wrapper."""

from __future__ import annotations

import subprocess
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
WRAPPER = REPOSITORY_ROOT / "hpc_jobs" / "active" / "hpc_pfp_benchmark.sh"


class PfpBenchmarkWrapperTests(unittest.TestCase):
    def test_wrapper_is_valid_bash(self):
        completed = subprocess.run(
            ["bash", "-n", str(WRAPPER)],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_optional_evidence_arrays_are_guarded_before_size_estimation(self):
        source = WRAPPER.read_text(encoding="utf-8")
        self.assertIn(
            'if [[ "$BENCHMARK_EVIDENCE_COUNT" -gt 0 ]]; then\n'
            '  for evidence in "${BENCHMARK_EVIDENCE[@]}"; '
            'do add_input_kb "$evidence"; done\n'
            "fi",
            source,
        )
        self.assertIn(
            'if [[ "$EMBEDDING_EVIDENCE_COUNT" -gt 0 ]]; then\n'
            '  for evidence in "${EMBEDDING_EVIDENCE[@]}"; '
            'do add_input_kb "$evidence"; done\n'
            "fi",
            source,
        )


if __name__ == "__main__":
    unittest.main()
