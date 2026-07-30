#!/usr/bin/env python3
"""Static contract tests for the contemporary follow-up Grid Engine jobs."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CAPTURE = (
    REPOSITORY_ROOT
    / "hpc_jobs"
    / "active"
    / "hpc_contemporary_followup_prediction_capture.sh"
)
ANALYSIS = (
    REPOSITORY_ROOT
    / "hpc_jobs"
    / "active"
    / "hpc_contemporary_followup_analysis.sh"
)
CENSUS = (
    REPOSITORY_ROOT
    / "hpc_jobs"
    / "active"
    / "hpc_contemporary_knowledge_cohort_census.sh"
)


class ContemporaryFollowupWrapperTests(unittest.TestCase):
    def test_wrappers_are_executable_and_valid_bash(self):
        for path in (CAPTURE, ANALYSIS):
            self.assertTrue(os.access(path, os.X_OK), path)
            completed = subprocess.run(
                ["bash", "-n", str(path)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_capture_is_inference_only_and_pinned_to_canonical_variant(self):
        source = CAPTURE.read_text(encoding="utf-8")
        self.assertIn("text-cutoff-2025-03-08__ppi-paper-faithful", source)
        self.assertNotIn("widened", source.lower())
        self.assertNotIn("train.py", source)
        self.assertIn("for split in valid test", source)
        self.assertIn('--evaluation-split "$split"', source)
        self.assertIn("scores_content_sha256", source)
        self.assertIn("Fresh test capture differs from accepted arrays", source)

    def test_capture_uses_gpu_and_zeus_but_analysis_is_cpu_only(self):
        capture = CAPTURE.read_text(encoding="utf-8")
        analysis = ANALYSIS.read_text(encoding="utf-8")
        self.assertIn("#$ -q gpu.q@zeus1.local,gpu.q@zeus2.local", capture)
        self.assertIn("#$ -l gpu=true", capture)
        self.assertIn("#$ -pe gpu 1", capture)
        self.assertNotIn("#$ -l gpu=true", analysis)
        self.assertIn("#$ -pe smp 2", analysis)

    def test_analysis_separates_specificity_and_calibration_inputs(self):
        source = ANALYSIS.read_text(encoding="utf-8")
        self.assertIn("--specificity-measure all_separate", source)
        self.assertIn("--bootstrap-replicates 2000", source)
        self.assertIn("--validation-prediction-manifest", source)
        self.assertIn("--test-prediction-manifest", source)
        self.assertIn("Paired capture is incomplete", source)
        self.assertIn("text-cutoff-2025-03-08__ppi-paper-faithful", source)
        self.assertNotIn("widened", source.lower())

    def test_census_uses_only_the_outer_transport_manifest_verifier(self):
        source = CENSUS.read_text(encoding="utf-8")
        completed = subprocess.run(
            ["bash", "-n", str(CENSUS)],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(source.count("manage_output_manifest.py verify"), 1)
        self.assertIn('cp -a "$ANALYSIS_OUTPUT" "$PUBLISH_STAGE/analysis"', source)
        self.assertIn(
            '--root "$PUBLISH_STAGE" --include-nested-control-files', source
        )


if __name__ == "__main__":
    unittest.main()
