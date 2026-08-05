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
    REPOSITORY_ROOT / "hpc_jobs" / "active" / "hpc_contemporary_followup_analysis.sh"
)
CENSUS = (
    REPOSITORY_ROOT
    / "hpc_jobs"
    / "active"
    / "hpc_contemporary_knowledge_cohort_census.sh"
)
CAFA3_CENSUS = (
    REPOSITORY_ROOT / "hpc_jobs" / "active" / "hpc_cafa3_knowledge_state_census.sh"
)
FORENSICS = REPOSITORY_ROOT / "hpc_jobs" / "active" / "hpc_benchmark_forensics.sh"
SPECIFICITY_COMPARISON = (
    REPOSITORY_ROOT / "hpc_jobs" / "active" / "hpc_specificity_comparison.sh"
)


class ContemporaryFollowupWrapperTests(unittest.TestCase):
    def test_wrappers_are_executable_and_valid_bash(self):
        for path in (
            CAPTURE,
            ANALYSIS,
            CAFA3_CENSUS,
            FORENSICS,
            SPECIFICITY_COMPARISON,
        ):
            self.assertTrue(os.access(path, os.X_OK), path)
            completed = subprocess.run(
                ["bash", "-n", str(path)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_nk_lk_forensics_preserves_per_ontology_split_contract(self):
        source = FORENSICS.read_text(encoding="utf-8")
        self.assertIn("contemporary-nk-lk", source)
        self.assertIn('SPLIT_OVERLAP_POLICY="per-ontology-disjoint"', source)
        self.assertIn("2025_01_to_2026_02_supervisor_nk_lk", source)
        self.assertIn("/outputs", source)

    def test_specificity_comparison_is_cpu_only_and_manifest_bound(self):
        source = SPECIFICITY_COMPARISON.read_text(encoding="utf-8")
        self.assertIn("compare_pfp_specificity_runs.py", source)
        self.assertIn("#$ -pe smp 1", source)
        self.assertNotIn("gpu=true", source)
        self.assertIn("activate_or_create_mmfp_env", source)
        self.assertIn('PYTHON_BIN="$(command -v python)"', source)
        self.assertIn('"$PYTHON_BIN" scripts/diagnostics/compare_pfp_specificity_runs.py', source)
        self.assertIn("RUN_COMPLETE.json", source)
        self.assertIn("output_manifest.json", source)

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

    def test_analysis_can_run_combined_prediction_diagnostics(self):
        source = ANALYSIS.read_text(encoding="utf-8")
        self.assertIn('"$ANALYSIS" == "diagnostics"', source)
        self.assertIn("evaluate_pfp_label_sensitivity.py", source)
        self.assertIn("evaluate_pfp_information_content.py", source)
        self.assertIn('$ANALYSIS_OUTPUT/root_exclusion', source)
        self.assertIn('$ANALYSIS_OUTPUT/specificity', source)

    def test_analysis_allows_provenance_safe_specificity_source_overrides(self):
        source = ANALYSIS.read_text(encoding="utf-8")
        self.assertIn("--source-run", source)
        self.assertIn("--source-label", source)
        self.assertIn("--obo", source)
        self.assertIn('"prediction_manifest_sha256"', source)
        self.assertIn('"obo_sha256"', source)
        self.assertIn('echo "Source label     : $SOURCE_LABEL"', source)

    def test_analysis_preserves_failure_diagnostics_before_scratch_cleanup(self):
        source = ANALYSIS.read_text(encoding="utf-8")
        self.assertIn('FAILURE_OUTPUT="${OUTPUT_DIR}.failed-${JOB_TOKEN}"', source)
        self.assertIn("publish_failure()", source)
        self.assertIn('cp -p "$LOG_FILE" "$failure_stage/logs/analysis.log"', source)
        self.assertIn(
            'cp -a "$ANALYSIS_OUTPUT" "$failure_stage/partial_analysis"', source
        )
        self.assertIn('mv "$failure_stage" "$FAILURE_OUTPUT"', source)
        self.assertEqual(source.count('2>&1 | tee "$LOG_FILE"'), 3)
        self.assertLess(
            source.index('publish_failure "$status"'),
            source.index('rm -rf -- "$WORK"'),
        )

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
        self.assertIn('--root "$PUBLISH_STAGE" --include-nested-control-files', source)

    def test_cafa3_census_uses_official_lists_and_published_csvs(self):
        source = CAFA3_CENSUS.read_text(encoding="utf-8")
        self.assertIn("canonical_cafa3", source)
        self.assertIn("cafa3_official/benchmark20171115.tar", source)
        self.assertIn(
            "d41dd38436461f4aa8072fca0e3c7476f36475e68cf1dc2555e78ccbdb15d70c",
            source,
        )
        self.assertIn("Official CAFA3 archive SHA-256 mismatch", source)
        self.assertNotIn("deepgoplus/data-cafa.tar.gz", source)
        self.assertIn("build_cafa3_knowledge_state_census.py", source)
        self.assertNotIn("goa_uniprot", source)
        self.assertNotIn("train.py", source)
        self.assertIn("#$ -pe smp 1", source)
        self.assertIn("activate_or_create_mmfp_env", source)
        self.assertIn('PYTHON_BIN="$(command -v python)"', source)
        self.assertIn("sys.version_info >= (3, 9)", source)
        self.assertIn('2>&1 | tee "$LOG_FILE"', source)


if __name__ == "__main__":
    unittest.main()
