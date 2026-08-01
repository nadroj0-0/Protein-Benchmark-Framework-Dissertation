from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = ROOT.parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cafa_benchmark_builder.snapshot import _git_commit


WRAPPER = (
    WORKSPACE_ROOT / "hpc_jobs" / "active"
    / "hpc_contemporary_temporal_benchmark.sh"
)


class SnapshotProvenanceTests(unittest.TestCase):
    def test_explicit_verified_revision_does_not_require_git(self):
        revision = "a" * 40
        with mock.patch.dict(
            os.environ, {"CAFA_BUILDER_FRAMEWORK_REVISION": revision}, clear=False
        ), mock.patch("subprocess.run") as run:
            self.assertEqual(_git_commit(), revision)
        run.assert_not_called()

    def test_invalid_explicit_revision_fails_loudly(self):
        with mock.patch.dict(
            os.environ, {"CAFA_BUILDER_FRAMEWORK_REVISION": "main"}, clear=False
        ):
            with self.assertRaisesRegex(ValueError, "complete lowercase"):
                _git_commit()

    def test_missing_git_returns_unknown_instead_of_crashing(self):
        with mock.patch.dict(os.environ, {}, clear=False), mock.patch(
            "subprocess.run", side_effect=FileNotFoundError("git")
        ):
            os.environ.pop("CAFA_BUILDER_FRAMEWORK_REVISION", None)
            self.assertIsNone(_git_commit())

    def test_non_repository_git_result_returns_unknown(self):
        completed = subprocess.CompletedProcess(
            args=["git", "rev-parse", "HEAD"], returncode=128, stdout=""
        )
        with mock.patch.dict(os.environ, {}, clear=False), mock.patch(
            "subprocess.run", return_value=completed
        ):
            os.environ.pop("CAFA_BUILDER_FRAMEWORK_REVISION", None)
            self.assertIsNone(_git_commit())

    def test_hpc_wrapper_exports_checked_out_revision_to_builder(self):
        text = WRAPPER.read_text(encoding="utf-8")
        self.assertIn(
            'export CAFA_BUILDER_FRAMEWORK_REVISION="$ACTUAL_FRAMEWORK_REVISION"',
            text,
        )
        self.assertIn('ACTUAL_FRAMEWORK_REVISION="$(git rev-parse HEAD)"', text)


if __name__ == "__main__":
    unittest.main()
