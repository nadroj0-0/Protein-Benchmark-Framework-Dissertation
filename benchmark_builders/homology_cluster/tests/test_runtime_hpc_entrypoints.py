from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PACKAGE_ROOT.parents[1]
DRIVER = (
    WORKSPACE_ROOT / "scripts" / "benchmark_generation"
    / "run_homology_cluster_runtime_hpc.sh"
)
PILOT = WORKSPACE_ROOT / "hpc_jobs" / "active" / "hpc_homology_cluster_runtime_pilot.sh"
ARRAY = WORKSPACE_ROOT / "hpc_jobs" / "active" / "hpc_homology_cluster_runtime_array.sh"


class RuntimeHPCEntrypointTests(unittest.TestCase):
    def _environment(self, root: Path, kind: str, task: str) -> tuple[dict[str, str], Path]:
        scratch = root / "scratch"
        results = root / "results"
        scratch.mkdir()
        results.mkdir()
        build = root / "fake-build.sh"
        build.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            "mkdir -p \"$1/benchmark\"\n"
            "printf 'complete\\n' > \"$1/benchmark/result.txt\"\n"
        )
        build.chmod(0o755)
        env = os.environ.copy()
        env.update({
            "HOMOLOGY_RUNTIME_TEST_MODE": "1",
            "HOMOLOGY_RUNTIME_TEST_BUILD_COMMAND": str(build),
            "HOMOLOGY_RUNTIME_KIND": kind,
            "SGE_TASK_ID": task,
            "JOB_ID": "fixture-job",
            "RUN_ID": f"{kind}-fixture",
            "UNIPROT_SOURCE_SCOPE": "sprot-only",
            "WORK_BASE": str(scratch),
            "RESULTS_ROOT": str(results),
            "FRAMEWORK_SOURCE_ROOT": str(WORKSPACE_ROOT),
        })
        return env, scratch

    @staticmethod
    def _run(env: dict[str, str]) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", str(DRIVER)],
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

    def test_pilot_copies_results_and_deletes_scratch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, scratch = self._environment(root, "pilot", "1")
            completed = self._run(env)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            self.assertEqual(len(list((root / "results").rglob("benchmark/result.txt"))), 1)
            summary = next((root / "results").rglob("logs/disk_usage_summary.tsv"))
            self.assertIn("peak_work_bytes", summary.read_text())
            samples = next((root / "results").rglob("logs/disk_usage.tsv"))
            self.assertIn("scratch-created", samples.read_text())
            self.assertFalse(list(scratch.iterdir()))

    def test_array_task_runs_without_any_pilot_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, scratch = self._environment(root, "array", "6")
            for name in list(env):
                if name.startswith("PILOT_") or name.startswith("EXPECTED_PILOT_"):
                    env.pop(name)
            completed = self._run(env)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            summary = next((root / "results").rglob("FINAL_RESULT_PATH.txt"))
            self.assertIn("task_6_identity_5", summary.read_text())
            self.assertFalse(list(scratch.iterdir()))

    def test_copy_failure_is_nonzero_and_still_deletes_scratch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, scratch = self._environment(root, "array", "2")
            fail_copy = root / "fail-copy.sh"
            fail_copy.write_text("#!/usr/bin/env bash\nexit 23\n")
            fail_copy.chmod(0o755)
            env["HOMOLOGY_RUNTIME_TEST_COPY_COMMAND"] = str(fail_copy)
            completed = self._run(env)
            self.assertEqual(completed.returncode, 74, completed.stdout)
            self.assertIn("Copy-back failed", completed.stdout)
            self.assertFalse(list(scratch.iterdir()))
            self.assertFalse(list((root / "results").rglob("*.partial-*")))

    def test_wrappers_are_thin_and_lock_the_expected_task_ranges(self):
        pilot = PILOT.read_text()
        array = ARRAY.read_text()
        self.assertIn("#$ -t 1\n", pilot)
        self.assertIn("HOMOLOGY_RUNTIME_KIND=pilot", pilot)
        self.assertIn("#$ -l tmem=8G", pilot)
        self.assertIn("#$ -l tscratch=38G", pilot)
        self.assertIn("#$ -l scratch0free=300G", pilot)
        self.assertIn("#$ -t 1-6\n", array)
        self.assertIn("#$ -tc 2\n", array)
        self.assertIn("HOMOLOGY_RUNTIME_KIND=array", array)
        for text in (pilot, array):
            self.assertIn("#$ -pe smp 8", text)
            self.assertIn("#$ -l tmem=8G", text)
            self.assertIn("#$ -l tscratch=38G", text)
            self.assertIn("#$ -l scratch0free=300G", text)
            self.assertIn("run_homology_cluster_runtime_hpc.sh", text)
            self.assertNotIn("wget ", text)
            self.assertNotIn("rm -rf", text)

    def test_runtime_driver_supports_the_clusters_legacy_git(self):
        driver = DRIVER.read_text()
        self.assertIn("git_in_dir()", driver)
        self.assertNotIn("git -C", driver)


if __name__ == "__main__":
    unittest.main()
