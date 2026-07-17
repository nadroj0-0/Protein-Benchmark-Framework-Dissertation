from __future__ import annotations

import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PACKAGE_ROOT.parents[1]
PILOT_LAUNCHER = (
    WORKSPACE_ROOT / "hpc_jobs" / "launchers" / "submit_homology_cluster_pilot.sh"
)
ARRAY_LAUNCHER = (
    WORKSPACE_ROOT / "hpc_jobs" / "launchers" / "submit_homology_cluster_array.sh"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ArrayLauncherTests(unittest.TestCase):
    def _environment(self, root: Path) -> tuple[dict[str, str], Path, Path]:
        root.mkdir(parents=True, exist_ok=True)
        bin_dir = root / "bin"
        bin_dir.mkdir()
        qsub_log = root / "qsub.log"
        auth_log = root / "authorize.log"
        fake_qsub = bin_dir / "qsub"
        fake_qsub.write_text(
            "#!/usr/bin/env bash\n"
            "printf '%s\\n' \"$@\" > \"$FAKE_QSUB_LOG\"\n"
            "printf 'DIAGNOSTIC_PILOT=%s\\n' \"${DIAGNOSTIC_PILOT:-unset}\" >> \"$FAKE_QSUB_LOG\"\n"
            "printf 'FIXTURE_MODE=%s\\n' \"${FIXTURE_MODE:-unset}\" >> \"$FAKE_QSUB_LOG\"\n"
            "exit \"${FAKE_QSUB_STATUS:-0}\"\n"
        )
        fake_qsub.chmod(0o755)
        fake_python = bin_dir / "python-wrapper"
        fake_python.write_text(
            "#!/usr/bin/env bash\n"
            "if [[ \"${1:-}\" == -c ]]; then exec \"$REAL_PYTHON\" \"$@\"; fi\n"
            "printf '%s\\n' \"$@\" > \"$FAKE_AUTHORIZE_LOG\"\n"
            "exit \"${FAKE_AUTHORIZE_STATUS:-0}\"\n"
        )
        fake_python.chmod(0o755)

        files = {}
        for name in ("manifest", "uniref", "idmapping", "sprot", "goa", "obo"):
            path = root / name
            path.write_text(f"fixture-{name}\n")
            files[name] = path
        results = root / "results"
        results.mkdir()
        env = os.environ.copy()
        env.update({
            "PATH": f"{bin_dir}{os.pathsep}{env.get('PATH', '')}",
            "PYTHON_BIN": str(fake_python),
            "REAL_PYTHON": sys.executable,
            "FAKE_QSUB_LOG": str(qsub_log),
            "FAKE_AUTHORIZE_LOG": str(auth_log),
            "UNIPROT_SOURCE_SCOPE": "sprot-only",
            "FRAMEWORK_REVISION": "a" * 40,
            "FROZEN_INPUT_MANIFEST": str(files["manifest"]),
            "UNIREF90_FASTA": str(files["uniref"]),
            "UNIREF90_FASTA_SHA256": _sha256(files["uniref"]),
            "IDMAPPING": str(files["idmapping"]),
            "IDMAPPING_SHA256": _sha256(files["idmapping"]),
            "UNIPROT_SPROT_SEQUENCES": str(files["sprot"]),
            "UNIPROT_SPROT_SEQUENCES_SHA256": _sha256(files["sprot"]),
            "GOA": str(files["goa"]),
            "GOA_SHA256": _sha256(files["goa"]),
            "GO_OBO": str(files["obo"]),
            "GO_OBO_SHA256": _sha256(files["obo"]),
            "RESULTS_ROOT": str(results),
            "EXPECTED_MMSEQS_VERSION": "15-6f452",
            "MMSEQS_BIN": "/shared/mmseqs",
            "SPLIT_POLICY": "sequence-balanced",
            "TRAINING_POPULATION": "annotated-only",
            "SEED": "0",
            "MIN_COUNT": "50",
            "RUN_ID": "launcher-fixture",
        })
        for name in (
            "UNIPROT_TREMBL_SEQUENCES",
            "UNIPROT_TREMBL_SEQUENCES_SHA256",
            "ATTRITION_OVERRIDE",
        ):
            env.pop(name, None)
        return env, qsub_log, auth_log

    @staticmethod
    def _run(path: Path, env: dict[str, str]) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", str(path)],
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

    def test_pilot_dry_run_prints_exact_contract_without_calling_qsub(self):
        with tempfile.TemporaryDirectory() as tmp:
            env, qsub_log, _ = self._environment(Path(tmp))
            env["DRY_RUN"] = "1"
            env["THREADS"] = "99"
            env["HOMOLOGY_BUILDER_COMMAND"] = "/poisoned/inherited/value"
            completed = self._run(PILOT_LAUNCHER, env)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            self.assertFalse(qsub_log.exists())
            self.assertIn("Launcher mode       : diagnostic-pilot", completed.stdout)
            self.assertIn("Task range          : 1", completed.stdout)
            self.assertIn("qsub -t 1 -pe smp 2", completed.stdout)
            self.assertIn("DIAGNOSTIC_PILOT=1", completed.stdout)
            self.assertNotIn("HOMOLOGY_BUILDER_COMMAND", completed.stdout)
            exported_line = next(
                line for line in completed.stdout.splitlines()
                if line.startswith("Exported variables")
            )
            self.assertNotIn("THREADS", exported_line)

    def test_pilot_fake_qsub_argv_and_failure_status_are_propagated(self):
        with tempfile.TemporaryDirectory() as tmp:
            env, qsub_log, _ = self._environment(Path(tmp))
            env["DRY_RUN"] = "0"
            env["FAKE_QSUB_STATUS"] = "17"
            completed = self._run(PILOT_LAUNCHER, env)
            self.assertEqual(completed.returncode, 17, completed.stdout)
            lines = qsub_log.read_text().splitlines()
            self.assertEqual(lines[:5], ["-t", "1", "-pe", "smp", "2"])
            self.assertIn("DIAGNOSTIC_PILOT=1", lines)
            self.assertIn("FIXTURE_MODE=0", lines)

    def _add_full_array_evidence(self, root: Path, env: dict[str, str]) -> None:
        pilot_run = root / "pilot-run"
        pilot_run.mkdir()
        evidence = {
            "ATTRITION_POLICY": root / "policy.json",
            "PILOT_APPROVAL": root / "approval.json",
            "PILOT_COMPLETION_MARKER": pilot_run / "RUN_COMPLETE.json",
            "PILOT_ATTRITION_REPORT": pilot_run / "attrition_report.json",
            "PILOT_TASK_CONTEXT": root / "hpc_task_context.json",
            "PILOT_MEASUREMENT_EVIDENCE": root / "pilot_measurements.json",
        }
        for name, path in evidence.items():
            path.write_text(f'{{"fixture":"{name}"}}\n')
            env[name] = str(path)
        env["PILOT_RUN_DIR"] = str(pilot_run)

    def test_full_array_authorization_then_fake_qsub_uses_locked_range_and_pe(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, qsub_log, auth_log = self._environment(root)
            self._add_full_array_evidence(root, env)
            env["DRY_RUN"] = "0"
            env["DIAGNOSTIC_PILOT"] = "1"
            env["FIXTURE_MODE"] = "1"
            completed = self._run(ARRAY_LAUNCHER, env)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            self.assertTrue(auth_log.is_file())
            authorization_args = auth_log.read_text()
            self.assertIn("authorize-array", authorization_args)
            self.assertIn("--pilot-measurement-evidence", authorization_args)
            lines = qsub_log.read_text().splitlines()
            self.assertEqual(lines[:5], ["-t", "1-6", "-pe", "smp", "2"])
            self.assertIn("DIAGNOSTIC_PILOT=0", lines)
            self.assertIn("FIXTURE_MODE=0", lines)
            export_argument = lines[6]
            self.assertNotIn("THREADS", export_argument)
            self.assertNotIn("HOMOLOGY_", export_argument)

    def test_full_array_authorization_failure_prevents_qsub(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, qsub_log, auth_log = self._environment(root)
            self._add_full_array_evidence(root, env)
            env["DRY_RUN"] = "0"
            env["FAKE_AUTHORIZE_STATUS"] = "9"
            completed = self._run(ARRAY_LAUNCHER, env)
            self.assertEqual(completed.returncode, 9, completed.stdout)
            self.assertTrue(auth_log.is_file())
            self.assertFalse(qsub_log.exists())

    def test_hash_mismatch_and_array_wide_override_fail_before_qsub(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, qsub_log, _ = self._environment(root)
            env["UNIREF90_FASTA_SHA256"] = "0" * 64
            completed = self._run(PILOT_LAUNCHER, env)
            self.assertEqual(completed.returncode, 2, completed.stdout)
            self.assertIn("does not match reviewed", completed.stdout)
            self.assertFalse(qsub_log.exists())

            env, qsub_log, _ = self._environment(root / "override")
            self._add_full_array_evidence(root / "override", env)
            env["ATTRITION_OVERRIDE"] = str(root / "override.json")
            completed = self._run(ARRAY_LAUNCHER, env)
            self.assertEqual(completed.returncode, 2, completed.stdout)
            self.assertIn("cannot be applied array-wide", completed.stdout)
            self.assertFalse(qsub_log.exists())

    def test_path_traversal_run_id_fails_before_qsub(self):
        with tempfile.TemporaryDirectory() as tmp:
            env, qsub_log, _ = self._environment(Path(tmp))
            env["RUN_ID"] = ".."
            completed = self._run(PILOT_LAUNCHER, env)
            self.assertEqual(completed.returncode, 2, completed.stdout)
            self.assertIn("contain an alphanumeric", completed.stdout)
            self.assertFalse(qsub_log.exists())


if __name__ == "__main__":
    unittest.main()
