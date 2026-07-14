from __future__ import annotations

import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PACKAGE_ROOT.parents[1]
HPC_WRAPPER = WORKSPACE_ROOT / "hpc_jobs" / "active" / "hpc_homology_cluster_benchmark.sh"
INNER_WRAPPER = (
    WORKSPACE_ROOT / "scripts" / "benchmark_generation"
    / "run_homology_cluster_benchmark.sh"
)


class HPCWrapperTests(unittest.TestCase):
    def _fixture(self, root: Path) -> tuple[dict[str, str], Path, Path]:
        framework = root / "framework"
        environment = root / "mmfp-env"
        home = root / "home"
        results = root / "results"
        work_base = root / "scratch"
        for path in (framework, environment, home, results, work_base):
            path.mkdir(parents=True)
        input_path = root / "u.fasta"
        input_path.write_bytes(b"AAAA")
        ready = root / "builder-ready"
        validation_log = root / "validations.log"
        signal_log = root / "signals.log"

        builder = root / "fake-builder.sh"
        builder.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            "trap 'printf TERM >> \"$SIGNAL_LOG\"; exit 143' TERM\n"
            "printf 'ready\\n' > \"$FAKE_READY\"\n"
            "if [[ \"${FAKE_MODE:-success}\" == fail ]]; then exit 7; fi\n"
            "if [[ \"${FAKE_MODE:-success}\" == sleep ]]; then while true; do sleep 1; done; fi\n"
            "run=\"$OUTPUT_ROOT/identity_30/$SPLIT_POLICY/$TRAINING_POPULATION/seed_$SEED/min_count_$MIN_COUNT\"\n"
            "mkdir -p \"$run\"\n"
            "printf '{\"complete\":true}\\n' > \"$run/RUN_COMPLETE.json\"\n"
            "printf 'payload\\n' > \"$run/payload.txt\"\n"
        )
        builder.chmod(0o755)
        validator = root / "fake-validator.sh"
        validator.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            "[[ -f \"$1/RUN_COMPLETE.json\" ]]\n"
            "printf '%s\\n' \"$1\" >> \"$VALIDATION_LOG\"\n"
        )
        validator.chmod(0o755)

        env = os.environ.copy()
        env.update({
            "HOME": str(home),
            "IDENTITY": "30",
            "JOB_ID": "fixture-job",
            "SPLIT_POLICY": "sequence-balanced",
            "TRAINING_POPULATION": "annotated-only",
            "SEED": "0",
            "MIN_COUNT": "1",
            "FIXTURE_MODE": "1",
            "WORK_BASE": str(work_base),
            "RESULTS_ROOT": str(results),
            "PERSISTENT_RESULTS_ROOT": str(results),
            "HOMOLOGY_WRAPPER_TEST_MODE": "1",
            "HOMOLOGY_FRAMEWORK_DIR": str(framework),
            "HOMOLOGY_SKIP_CONDA": "1",
            "HOMOLOGY_BUILDER_COMMAND": str(builder),
            "HOMOLOGY_VALIDATE_COMMAND": str(validator),
            "MMFP_ENV_DIR": str(environment),
            "CONDA_PREFIX": str(environment),
            "PYTHON_BIN": "python3",
            "UNIREF90_FASTA": str(input_path),
            "CLUSTER_ASSIGNMENTS": str(input_path),
            "FAKE_READY": str(ready),
            "VALIDATION_LOG": str(validation_log),
            "SIGNAL_LOG": str(signal_log),
            # Two staged four-byte inputs + 8x four MMseq estimate = exact 40-byte boundary.
            "HOMOLOGY_SCRATCH_FREE_BYTES_OVERRIDE": "40",
            # Two-times the eight staged input bytes = exact 16-byte boundary.
            "HOMOLOGY_PERSISTENT_FREE_BYTES_OVERRIDE": "16",
            "HOMOLOGY_PRECOPY_FREE_BYTES_OVERRIDE": "1000000000",
            "HOMOLOGY_SIGNAL_GRACE_SECONDS": "1",
        })
        return env, results, work_base / "homology_cluster_fixture-job_30"

    def test_success_validates_scratch_and_copy_then_atomically_publishes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, results, work = self._fixture(root)
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout)
            finals = [path for path in results.iterdir() if ".partial-" not in path.name]
            self.assertEqual(len(finals), 1)
            self.assertFalse(finals[0].name.endswith(".failed"))
            self.assertTrue(list(finals[0].rglob("RUN_COMPLETE.json")))
            self.assertEqual(len((root / "validations.log").read_text().splitlines()), 2)
            self.assertFalse(work.exists())
            self.assertFalse(list(results.glob("*.partial-*")))
            capacity = json.loads(
                next(finals[0].glob("logs/hpc_capacity_preflight.json")).read_text()
            )
            self.assertEqual(capacity["local_input_bytes"], 8)
            self.assertEqual(capacity["scratch_required_bytes"], 40)

    def test_builder_failure_publishes_failed_diagnostics_without_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, results, work = self._fixture(root)
            env["FAKE_MODE"] = "fail"
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertEqual(completed.returncode, 7, completed.stdout)
            failed = list(results.glob("*.failed"))
            self.assertEqual(len(failed), 1)
            self.assertTrue((failed[0] / "FAILURE.json").is_file())
            self.assertFalse(list(failed[0].rglob("RUN_COMPLETE.json")))
            failure = json.loads((failed[0] / "FAILURE.json").read_text())
            self.assertEqual(failure["exit_code"], 7)
            self.assertEqual(failure["failure_stage"], "benchmark-builder")
            self.assertFalse(work.exists())
            self.assertFalse(list(results.glob("*.partial-*")))

    def test_sigterm_and_sigint_terminate_child_and_publish_only_failed_diagnostics(self):
        for sent_signal, expected_status in ((signal.SIGTERM, 143), (signal.SIGINT, 130)):
            with self.subTest(signal=sent_signal), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                env, results, _ = self._fixture(root)
                env["FAKE_MODE"] = "sleep"
                process = subprocess.Popen(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                )
                deadline = time.time() + 10
                while not (root / "builder-ready").exists() and time.time() < deadline:
                    time.sleep(0.05)
                self.assertTrue((root / "builder-ready").exists())
                process.send_signal(sent_signal)
                output, _ = process.communicate(timeout=15)
                self.assertEqual(process.returncode, expected_status, output)
                self.assertIn("TERM", (root / "signals.log").read_text())
                failed = list(results.glob("*.failed"))
                self.assertEqual(len(failed), 1)
                failure = json.loads((failed[0] / "FAILURE.json").read_text())
                self.assertEqual(failure["signal"], signal.Signals(sent_signal).name[3:])
                self.assertFalse(list(failed[0].rglob("RUN_COMPLETE.json")))

    def test_copy_failure_preserves_exact_scratch_path_and_never_publishes_final(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, results, work = self._fixture(root)
            copier = root / "fail-copy.sh"
            copier.write_text(
                "#!/usr/bin/env bash\n"
                "set -e\n"
                "cp -a \"$1/.\" \"$2/\"\n"
                "exit 9\n"
            )
            copier.chmod(0o755)
            env["HOMOLOGY_COPY_COMMAND"] = str(copier)
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertTrue(work.exists())
            self.assertIn(str(work), completed.stdout)
            failed = list(results.glob("*.failed"))
            self.assertEqual(len(failed), 1)
            self.assertFalse(list(failed[0].rglob("RUN_COMPLETE.json")))
            failure = json.loads((failed[0] / "FAILURE.json").read_text())
            self.assertEqual(failure["exit_code"], 9)
            self.assertEqual(failure["failure_stage"], "copy-to-persistent-partial")
            self.assertFalse(list(results.glob("*.partial-*")))

    def test_capacity_boundaries_and_wrong_conda_prefix_fail_early(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for label, variable, value, expected in (
                ("scratch", "HOMOLOGY_SCRATCH_FREE_BYTES_OVERRIDE", "39", "Scratch"),
                ("persistent", "HOMOLOGY_PERSISTENT_FREE_BYTES_OVERRIDE", "15", "Persistent"),
            ):
                with self.subTest(label=label):
                    case = root / label
                    case.mkdir()
                    env, _, _ = self._fixture(case)
                    env[variable] = value
                    completed = subprocess.run(
                        ["bash", str(HPC_WRAPPER)], env=env, text=True,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                    )
                    self.assertNotEqual(completed.returncode, 0)
                    self.assertIn(expected, completed.stdout)
                    self.assertFalse((case / "builder-ready").exists())
                    failed = list((case / "results").glob("*.failed"))
                    self.assertEqual(len(failed), 1)
                    capacity = json.loads(
                        (failed[0] / "logs" / "hpc_capacity_preflight.json").read_text()
                    )
                    self.assertEqual(capacity["local_input_bytes"], 8)
                    self.assertEqual(capacity["scratch_required_bytes"], 40)

            case = root / "wrong-prefix"
            case.mkdir()
            env, _, _ = self._fixture(case)
            other = case / "other-env"
            other.mkdir()
            env["CONDA_PREFIX"] = str(other)
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("CONDA_PREFIX mismatch", completed.stdout)
            self.assertFalse((case / "builder-ready").exists())

    def test_exact_precopy_capacity_failure_records_observation_and_no_success_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, results, _ = self._fixture(root)
            env["HOMOLOGY_PRECOPY_FREE_BYTES_OVERRIDE"] = "0"
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertNotEqual(completed.returncode, 0, completed.stdout)
            self.assertTrue((root / "builder-ready").is_file())
            failed = list(results.glob("*.failed"))
            self.assertEqual(len(failed), 1)
            capacity_path = failed[0] / "logs" / "hpc_capacity_precopy.json"
            self.assertTrue(capacity_path.is_file())
            capacity = json.loads(capacity_path.read_text())
            self.assertEqual(capacity["persistent_free_bytes"], 0)
            self.assertGreater(capacity["persistent_required_bytes"], 0)
            self.assertFalse(list(failed[0].rglob("RUN_COMPLETE.json")))

    def test_destination_roots_must_resolve_to_the_same_filesystem_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, _, _ = self._fixture(root)
            different = root / "different-results"
            different.mkdir()
            env["PERSISTENT_RESULTS_ROOT"] = str(different)
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertEqual(completed.returncode, 2, completed.stdout)
            self.assertIn("must resolve to the same destination", completed.stdout)
            self.assertFalse((root / "builder-ready").exists())

    def test_real_conda_hook_activates_exact_mmfp_environment_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, _, _ = self._fixture(root)
            environment = Path(env["MMFP_ENV_DIR"])
            bin_dir = environment / "bin"
            bin_dir.mkdir()
            (bin_dir / "python").symlink_to(Path(sys.executable).resolve())
            activation_log = root / "conda-activation.log"
            fake_conda = root / "fake-conda"
            fake_conda.write_text(
                "#!/usr/bin/env bash\n"
                "if [[ \"$1\" == shell.bash && \"$2\" == hook ]]; then\n"
                "  printf '%s\\n' \\\n"
                "    'conda() {' \\\n"
                "    '  printf \"%s\\n\" \"$*\" >> \"$FAKE_CONDA_LOG\"' \\\n"
                "    '  if [[ \"$1\" == activate ]]; then' \\\n"
                "    '    export CONDA_PREFIX=\"$2\"' \\\n"
                "    '    export PATH=\"$2/bin:$PATH\"' \\\n"
                "    '  fi' \\\n"
                "    '}'\n"
                "else\n"
                "  exit 2\n"
                "fi\n"
            )
            fake_conda.chmod(0o755)
            env.pop("HOMOLOGY_SKIP_CONDA")
            env.pop("CONDA_PREFIX")
            env["CONDA_EXE"] = str(fake_conda)
            env["FAKE_CONDA_LOG"] = str(activation_log)
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout)
            self.assertEqual(
                activation_log.read_text().strip(), f"activate {environment}"
            )

    def test_production_rejects_test_mode_and_fixture_rejects_unenabled_overrides(self):
        production = os.environ.copy()
        production.update({
            "IDENTITY": "30", "FIXTURE_MODE": "0", "HOMOLOGY_WRAPPER_TEST_MODE": "1",
        })
        completed = subprocess.run(
            ["bash", str(HPC_WRAPPER)], env=production, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("permitted only with FIXTURE_MODE=1", completed.stdout)

        for override in (
            "HOMOLOGY_BUILDER_COMMAND", "HOMOLOGY_VALIDATE_COMMAND",
            "HOMOLOGY_COPY_COMMAND", "HOMOLOGY_FRAMEWORK_DIR", "HOMOLOGY_SKIP_CONDA",
            "HOMOLOGY_SCRATCH_FREE_BYTES_OVERRIDE",
            "HOMOLOGY_PERSISTENT_FREE_BYTES_OVERRIDE",
            "HOMOLOGY_PRECOPY_FREE_BYTES_OVERRIDE", "HOMOLOGY_SIGNAL_GRACE_SECONDS",
            "CAPACITY_PYTHON",
        ):
            with self.subTest(override=override):
                env = os.environ.copy()
                env.update({"IDENTITY": "30", "FIXTURE_MODE": "1", override: "1"})
                completed = subprocess.run(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertEqual(completed.returncode, 2, completed.stdout)
                self.assertIn("test-only override", completed.stdout)

    def test_inner_shell_sigint_terminates_its_background_process(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ready = root / "ready"
            signal_log = root / "signals"
            fake_python = root / "fake-python"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "trap 'printf TERM >> \"$SIGNAL_LOG\"; exit 143' TERM\n"
                "printf ready > \"$FAKE_READY\"\n"
                "while true; do sleep 1; done\n"
            )
            fake_python.chmod(0o755)
            input_path = root / "input"
            input_path.write_text("x")
            env = os.environ.copy()
            env.update({
                "PYTHON_BIN": str(fake_python), "IDENTITY": "30", "FIXTURE_MODE": "1",
                "UNIREF90_FASTA": str(input_path), "IDMAPPING": str(input_path),
                "UNIPROT_SEQUENCES": str(input_path), "GOA": str(input_path),
                "GO_OBO": str(input_path), "CLUSTER_ASSIGNMENTS": str(input_path),
                "OUTPUT_ROOT": str(root / "outputs"), "TEMP_DIR": str(root / "temp"),
                "FAKE_READY": str(ready), "SIGNAL_LOG": str(signal_log),
            })
            process = subprocess.Popen(
                ["bash", str(INNER_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            )
            deadline = time.time() + 10
            while not ready.exists() and time.time() < deadline:
                time.sleep(0.05)
            self.assertTrue(ready.exists())
            process.send_signal(signal.SIGINT)
            output, _ = process.communicate(timeout=15)
            self.assertEqual(process.returncode, 130, output)
            self.assertIn("TERM", signal_log.read_text())

    def test_all_members_fails_before_workspace_or_builder_work(self):
        env = os.environ.copy()
        env.update({"IDENTITY": "30", "TRAINING_POPULATION": "all-cluster-members"})
        completed = subprocess.run(
            ["bash", str(HPC_WRAPPER)], env=env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("unsupported", completed.stdout)


if __name__ == "__main__":
    unittest.main()
