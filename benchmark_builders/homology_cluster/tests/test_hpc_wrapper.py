from __future__ import annotations

import json
import hashlib
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
            "identity_dir=identity_$IDENTITY\n"
            "[[ \"$IDENTITY\" != 5 ]] || identity_dir=identity_05\n"
            "run=\"$OUTPUT_ROOT/source_sprot-only/framework_fixture/$identity_dir/"
            "$SPLIT_POLICY/$TRAINING_POPULATION/seed_$SEED/min_count_$MIN_COUNT\"\n"
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
            "count=$(awk 'END {print NR}' \"$VALIDATION_LOG\")\n"
            "if [[ -n \"${FAKE_VALIDATE_FAIL_AT:-}\" && \"$count\" == \"$FAKE_VALIDATE_FAIL_AT\" ]]; then exit 8; fi\n"
        )
        validator.chmod(0o755)

        env = os.environ.copy()
        env.update({
            "HOME": str(home),
            "IDENTITY": "30",
            "SGE_TASK_ID": "1",
            "JOB_ID": "fixture-job",
            "NSLOTS": "8",
            "SPLIT_POLICY": "sequence-balanced",
            "TRAINING_POPULATION": "annotated-only",
            "UNIPROT_SOURCE_SCOPE": "sprot-only",
            "RUN_ID": "fixture-run",
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
            "UNIPROT_SPROT_SEQUENCES": str(input_path),
            "CLUSTER_ASSIGNMENTS": str(input_path),
            "FAKE_READY": str(ready),
            "VALIDATION_LOG": str(validation_log),
            "SIGNAL_LOG": str(signal_log),
            # Three staged four-byte inputs + 8x four MMseq estimate = exact 44-byte boundary.
            "HOMOLOGY_SCRATCH_FREE_BYTES_OVERRIDE": "44",
            # Two-times the twelve staged input bytes = exact 24-byte boundary.
            "HOMOLOGY_PERSISTENT_FREE_BYTES_OVERRIDE": "24",
            "HOMOLOGY_PRECOPY_FREE_BYTES_OVERRIDE": "1000000000",
            "HOMOLOGY_SIGNAL_GRACE_SECONDS": "1",
        })
        return (
            env,
            results,
            work_base / "homology_cluster_fixture-job_1_30_sprot-only_fixture-run",
        )

    def _production_revision_fixture(self, root: Path) -> dict[str, str]:
        bin_dir = root / "bin"
        environment = root / "mmfp-env"
        environment_bin = environment / "bin"
        home = root / "home"
        results = root / "results"
        scratch = root / "scratch"
        for path in (bin_dir, environment_bin, home, results, scratch):
            path.mkdir(parents=True)
        (environment_bin / "python").symlink_to(Path(sys.executable).resolve())
        fake_conda = bin_dir / "conda"
        fake_conda.write_text(
            "#!/usr/bin/env bash\n"
            "if [[ \"$1\" == shell.bash && \"$2\" == hook ]]; then\n"
            "  printf '%s\\n' \\\n"
            "    'conda() {' \\\n"
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
        fake_git = bin_dir / "git"
        fake_git.write_text(
            "#!/usr/bin/env bash\n"
            "if [[ \"$1\" == clone ]]; then\n"
            "  mkdir -p \"$3/scripts\"\n"
            "  printf 'validate_mmfp_env() { :; }\\n' > \"$3/scripts/reproduction_common.sh\"\n"
            "  exit 0\n"
            "fi\n"
            "if [[ \"$1\" != -C ]]; then exit 2; fi\n"
            "case \"$3\" in\n"
            "  checkout) exit \"${FAKE_GIT_CHECKOUT_STATUS:-0}\" ;;\n"
            "  rev-parse) printf '%s\\n' \"$FAKE_GIT_HEAD\" ;;\n"
            "  symbolic-ref) exit 1 ;;\n"
            "  status) [[ -z \"${FAKE_GIT_DIRTY:-}\" ]] || printf ' M changed\\n' ;;\n"
            "  *) exit 2 ;;\n"
            "esac\n"
        )
        fake_git.chmod(0o755)
        source = root / "input"
        source.write_text("AAAA")
        manifest = root / "frozen.json"
        manifest.write_text('{"inputs":[]}\n')
        manifest_hash = hashlib.sha256(manifest.read_bytes()).hexdigest()
        env = os.environ.copy()
        env.update({
            "PATH": f"{bin_dir}{os.pathsep}{env.get('PATH', '')}",
            "HOME": str(home),
            "SGE_TASK_ID": "1",
            "JOB_ID": "12345",
            "NSLOTS": "8",
            "REQUESTED_SLOTS": "8",
            "UNIPROT_SOURCE_SCOPE": "sprot-only",
            "RUN_ID": "revision-fixture",
            "DIAGNOSTIC_PILOT": "1",
            "FIXTURE_MODE": "0",
            "NO_DOWNLOADS": "1",
            "FROZEN_INPUT_MANIFEST": str(manifest),
            "EXPECTED_FROZEN_INPUT_MANIFEST_SHA256": manifest_hash,
            "EXPECTED_MMSEQS_VERSION": "15-6f452",
            "UNIREF90_FASTA": str(source),
            "IDMAPPING": str(source),
            "UNIPROT_SPROT_SEQUENCES": str(source),
            "GOA": str(source),
            "GO_OBO": str(source),
            "RESULTS_ROOT": str(results),
            "PERSISTENT_RESULTS_ROOT": str(results),
            "WORK_BASE": str(scratch),
            "MMFP_ENV_DIR": str(environment),
            "CONDA_EXE": str(fake_conda),
            "FRAMEWORK_REPO_URL": "fixture://framework",
            "FAKE_GIT_HEAD": "a" * 40,
            "MMSEQS_BIN": str(root / "missing-mmseqs"),
        })
        for name in (
            "HOMOLOGY_WRAPPER_TEST_MODE", "HOMOLOGY_FRAMEWORK_DIR",
            "HOMOLOGY_SKIP_CONDA", "CONDA_PREFIX", "THREADS",
        ):
            env.pop(name, None)
        return env

    @staticmethod
    def _failed_roots(results: Path) -> list[Path]:
        return sorted(path for path in results.rglob("*.failed") if path.is_dir())

    @staticmethod
    def _partial_roots(results: Path) -> list[Path]:
        return sorted(path for path in results.rglob("*.partial-*") if path.is_dir())

    def test_success_validates_scratch_and_copy_then_atomically_publishes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, results, work = self._fixture(root)
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout)
            markers = list(results.rglob("RUN_COMPLETE.json"))
            self.assertEqual(len(markers), 1)
            final_root = next(results.rglob("logs/hpc_capacity_preflight.json")).parent.parent
            self.assertFalse(final_root.name.endswith(".failed"))
            self.assertTrue(markers[0].is_relative_to(final_root))
            self.assertEqual(len((root / "validations.log").read_text().splitlines()), 2)
            self.assertFalse(work.exists())
            self.assertFalse(self._partial_roots(results))
            capacity = json.loads(
                (final_root / "logs" / "hpc_capacity_preflight.json").read_text()
            )
            self.assertEqual(capacity["local_input_bytes"], 12)
            self.assertEqual(capacity["scratch_required_bytes"], 44)
            context = json.loads((final_root / "logs" / "hpc_task_context.json").read_text())
            self.assertEqual(context["sge_task_id"], 1)
            self.assertEqual(context["identity_percent"], 30)
            self.assertEqual(context["nslots"], 8)
            self.assertEqual(context["mmseqs_threads"], 8)

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
            failed = self._failed_roots(results)
            self.assertEqual(len(failed), 1)
            self.assertTrue((failed[0] / "FAILURE.json").is_file())
            self.assertFalse(list(failed[0].rglob("RUN_COMPLETE.json")))
            failure = json.loads((failed[0] / "FAILURE.json").read_text())
            self.assertEqual(failure["exit_code"], 7)
            self.assertEqual(failure["failure_stage"], "benchmark-builder")
            self.assertFalse(work.exists())
            self.assertFalse(self._partial_roots(results))

    def test_sigterm_and_sigint_terminate_child_and_publish_only_failed_diagnostics(self):
        for sent_signal, expected_status in ((signal.SIGTERM, 143), (signal.SIGINT, 130)):
            with self.subTest(signal=sent_signal), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                env, results, work = self._fixture(root)
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
                failed = self._failed_roots(results)
                self.assertEqual(len(failed), 1)
                failure = json.loads((failed[0] / "FAILURE.json").read_text())
                self.assertEqual(failure["signal"], signal.Signals(sent_signal).name[3:])
                self.assertFalse(list(failed[0].rglob("RUN_COMPLETE.json")))
                self.assertFalse(work.exists())

    def test_copy_failure_cleans_owned_scratch_and_never_publishes_final(self):
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
            self.assertFalse(work.exists())
            self.assertIn("mandatory scratch cleanup", completed.stdout)
            failed = self._failed_roots(results)
            self.assertEqual(len(failed), 1)
            self.assertFalse(list(failed[0].rglob("RUN_COMPLETE.json")))
            failure = json.loads((failed[0] / "FAILURE.json").read_text())
            self.assertEqual(failure["exit_code"], 9)
            self.assertEqual(failure["failure_stage"], "copy-to-persistent-partial")
            self.assertFalse(self._partial_roots(results))

    def test_capacity_boundaries_and_wrong_conda_prefix_fail_early(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for label, variable, value, expected in (
                ("scratch", "HOMOLOGY_SCRATCH_FREE_BYTES_OVERRIDE", "43", "Scratch"),
                ("persistent", "HOMOLOGY_PERSISTENT_FREE_BYTES_OVERRIDE", "23", "Persistent"),
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
                    failed = self._failed_roots(case / "results")
                    self.assertEqual(len(failed), 1)
                    capacity = json.loads(
                        (failed[0] / "logs" / "hpc_capacity_preflight.json").read_text()
                    )
                    self.assertEqual(capacity["local_input_bytes"], 12)
                    self.assertEqual(capacity["scratch_required_bytes"], 44)

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
            failed = self._failed_roots(results)
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
            "IDENTITY": "30", "SGE_TASK_ID": "1", "UNIPROT_SOURCE_SCOPE": "sprot-only",
            "FIXTURE_MODE": "0", "HOMOLOGY_WRAPPER_TEST_MODE": "1",
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
                env.update({
                    "IDENTITY": "30", "SGE_TASK_ID": "1", "FIXTURE_MODE": "1",
                    override: "1",
                })
                completed = subprocess.run(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertEqual(completed.returncode, 2, completed.stdout)
                self.assertIn("test-only override", completed.stdout)

    def test_production_requires_one_full_lowercase_framework_commit(self):
        invalid_revisions = (None, "main", "v1.0.0", "a" * 12, "A" * 40, "g" * 40)
        for revision in invalid_revisions:
            with self.subTest(revision=revision), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                env = self._production_revision_fixture(root)
                if revision is None:
                    env.pop("FRAMEWORK_REVISION", None)
                else:
                    env["FRAMEWORK_REVISION"] = revision
                completed = subprocess.run(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertEqual(completed.returncode, 2, completed.stdout)
                self.assertIn("exactly 40 lowercase hex", completed.stdout)
                self.assertFalse(list((root / "scratch").glob("homology_cluster_*")))

    def test_checkout_failure_head_mismatch_and_dirty_checkout_fail_before_mmseqs(self):
        cases = (
            ("checkout-failure", {"FAKE_GIT_CHECKOUT_STATUS": "5"}, None),
            ("head-mismatch", {"FAKE_GIT_HEAD": "b" * 40}, "differs"),
            ("dirty-checkout", {"FAKE_GIT_DIRTY": "1"}, "dirty"),
        )
        for label, changes, expected in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                env = self._production_revision_fixture(root)
                env["FRAMEWORK_REVISION"] = "a" * 40
                env.update(changes)
                completed = subprocess.run(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertNotEqual(completed.returncode, 0, completed.stdout)
                if expected is not None:
                    self.assertIn(expected, completed.stdout)
                self.assertNotIn("MMseqs2 is not executable", completed.stdout)
                self.assertFalse(list((root / "scratch").glob("homology_cluster_*")))

    def test_valid_detached_full_commit_passes_revision_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env = self._production_revision_fixture(root)
            env["FRAMEWORK_REVISION"] = "a" * 40
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertNotEqual(completed.returncode, 0, completed.stdout)
            self.assertIn("MMseqs2 is not executable", completed.stdout)
            self.assertNotIn("differs from FRAMEWORK_REVISION", completed.stdout)
            self.assertNotIn("checkout is dirty", completed.stdout)
            self.assertFalse(list((root / "scratch").glob("homology_cluster_*")))

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
                "UNIPROT_SOURCE_SCOPE": "sprot-only",
                "UNIREF90_FASTA": str(input_path), "IDMAPPING": str(input_path),
                "UNIPROT_SPROT_SEQUENCES": str(input_path), "GOA": str(input_path),
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
        env.update({
            "IDENTITY": "30", "SGE_TASK_ID": "1",
            "TRAINING_POPULATION": "all-cluster-members",
        })
        completed = subprocess.run(
            ["bash", str(HPC_WRAPPER)], env=env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("unsupported", completed.stdout)

    def test_array_task_mapping_and_slot_mismatch_fail_before_builder(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, _, _ = self._fixture(root)
            env["SGE_TASK_ID"] = "2"
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertEqual(completed.returncode, 2, completed.stdout)
            self.assertIn("conflicts with locked", completed.stdout)
            self.assertFalse((root / "builder-ready").exists())

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, _, _ = self._fixture(root)
            env["THREADS"] = "4"
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertNotEqual(completed.returncode, 0, completed.stdout)
            self.assertIn("must equal scheduler-provided NSLOTS=8", completed.stdout)
            self.assertFalse((root / "builder-ready").exists())

    def test_all_six_locked_array_mappings_and_invalid_task_ids(self):
        for task_id, identity in enumerate((30, 25, 20, 15, 10, 5), start=1):
            with self.subTest(task_id=task_id), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                env, results, _ = self._fixture(root)
                env["SGE_TASK_ID"] = str(task_id)
                env["IDENTITY"] = str(identity)
                completed = subprocess.run(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertEqual(completed.returncode, 0, completed.stdout)
                context_path = next(results.rglob("hpc_task_context.json"))
                context = json.loads(context_path.read_text())
                self.assertEqual(context["sge_task_id"], task_id)
                self.assertEqual(context["identity_percent"], identity)

        for task_id in (None, "0", "-1", "7", "not-an-integer"):
            with self.subTest(invalid=task_id), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                env, _, _ = self._fixture(root)
                if task_id is None:
                    env.pop("SGE_TASK_ID")
                else:
                    env["SGE_TASK_ID"] = task_id
                completed = subprocess.run(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertEqual(completed.returncode, 2, completed.stdout)
                self.assertIn("SGE_TASK_ID must be one integer from 1 through 6", completed.stdout)

    def test_validation_failures_clean_scratch_and_publish_only_failure_diagnostics(self):
        for fail_at, expected_stage in (("1", "scratch-validation"), ("2", "copied-validation")):
            with self.subTest(fail_at=fail_at), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                env, results, work = self._fixture(root)
                env["FAKE_VALIDATE_FAIL_AT"] = fail_at
                completed = subprocess.run(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertEqual(completed.returncode, 8, completed.stdout)
                failed = self._failed_roots(results)
                self.assertEqual(len(failed), 1)
                failure = json.loads((failed[0] / "FAILURE.json").read_text())
                self.assertEqual(failure["failure_stage"], expected_stage)
                self.assertFalse(list(failed[0].rglob("RUN_COMPLETE.json")))
                self.assertFalse(work.exists())

    def test_preexisting_scratch_and_symlink_are_never_seized_or_removed(self):
        for use_symlink in (False, True):
            with self.subTest(use_symlink=use_symlink), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                env, _, work = self._fixture(root)
                sentinel_root = root / "other-task"
                if use_symlink:
                    sentinel_root.mkdir()
                    (sentinel_root / "sentinel").write_text("keep")
                    work.symlink_to(sentinel_root, target_is_directory=True)
                else:
                    work.mkdir()
                    (work / "sentinel").write_text("keep")
                completed = subprocess.run(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertEqual(completed.returncode, 2, completed.stdout)
                self.assertIn("Refusing to reuse pre-existing task scratch path", completed.stdout)
                kept = sentinel_root / "sentinel" if use_symlink else work / "sentinel"
                self.assertEqual(kept.read_text(), "keep")

    def test_existing_persistent_publication_is_not_mutated_on_collision(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, results, work = self._fixture(root)
            final = (
                results / "source_sprot-only" / "framework_fixture" / "fixture-run"
                / "job_fixture-job" / "task_1"
                / "identity_30_sequence-balanced_annotated-only_seed_0_min_count_1"
            )
            final.mkdir(parents=True)
            marker = final / "RUN_COMPLETE.json"
            marker.write_text('{"complete":true,"sentinel":"prior-run"}\n')
            payload = final / "payload.bin"
            payload.write_bytes(b"prior-publication")
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertNotEqual(completed.returncode, 0, completed.stdout)
            self.assertIn("Refusing to overwrite pre-existing persistent task path", completed.stdout)
            self.assertTrue(marker.is_file())
            self.assertEqual(payload.read_bytes(), b"prior-publication")
            self.assertFalse(self._failed_roots(results))
            self.assertFalse(work.exists())

    def test_relative_and_root_work_bases_fail_before_creating_scratch(self):
        for work_base in ("relative-scratch", "/"):
            with self.subTest(work_base=work_base), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                env, _, _ = self._fixture(root)
                env["WORK_BASE"] = work_base
                completed = subprocess.run(
                    ["bash", str(HPC_WRAPPER)], env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertEqual(completed.returncode, 2, completed.stdout)
                self.assertIn("existing absolute non-root directory", completed.stdout)

    def test_cleanup_rejects_a_mismatched_ownership_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, results, work = self._fixture(root)
            corruptor = root / "corrupt-owner.sh"
            corruptor.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "work=$(dirname \"$(dirname \"$OUTPUT_ROOT\")\")\n"
                "printf 'not-the-task-owner\\n' > \"$work/.homology-task-owner\"\n"
                "exit 7\n"
            )
            corruptor.chmod(0o755)
            env["HOMOLOGY_BUILDER_COMMAND"] = str(corruptor)
            completed = subprocess.run(
                ["bash", str(HPC_WRAPPER)], env=env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            )
            self.assertEqual(completed.returncode, 7, completed.stdout)
            self.assertIn("Refusing scratch cleanup with mismatched ownership marker", completed.stdout)
            self.assertTrue(work.exists())
            self.assertEqual(len(self._failed_roots(results)), 1)


if __name__ == "__main__":
    unittest.main()
