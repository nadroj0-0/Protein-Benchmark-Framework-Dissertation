from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
HELPER = REPOSITORY_ROOT / "hpc_jobs" / "qsub_with_notifications.sh"


class QsubWithNotificationsTests(unittest.TestCase):
    def _run(
        self,
        root: Path,
        email_file: Path,
        *,
        qsub_exit: int = 0,
    ) -> tuple[subprocess.CompletedProcess[str], list[str]]:
        capture = root / "qsub-arguments.txt"
        fake_qsub = root / "qsub"
        fake_qsub.write_text(
            "#!/usr/bin/env bash\n"
            "printf '%s\\n' \"$@\" > \"$QSUB_CAPTURE\"\n"
            "exit \"$QSUB_EXIT\"\n"
        )
        fake_qsub.chmod(0o755)

        env = os.environ.copy()
        env.update(
            {
                "GRID_ENGINE_NOTIFY_EMAIL_FILE": str(email_file),
                "GRID_ENGINE_QSUB_BIN": str(fake_qsub),
                "QSUB_CAPTURE": str(capture),
                "QSUB_EXIT": str(qsub_exit),
            }
        )
        completed = subprocess.run(
            [str(HELPER), "-N", "fixture", "worker.sh", "--flag"],
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        arguments = capture.read_text().splitlines()
        return completed, arguments

    def test_valid_email_enables_every_supported_mail_event(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            email_file = root / "email"
            email_file.write_text("person@example.org\n")

            completed, arguments = self._run(root, email_file)

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(
                arguments,
                [
                    "-m",
                    "beas",
                    "-M",
                    "person@example.org",
                    "-N",
                    "fixture",
                    "worker.sh",
                    "--flag",
                ],
            )
            self.assertIn("notifications enabled", completed.stderr)
            self.assertNotIn("person@example.org", completed.stderr)

    def test_missing_file_submits_normally(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)

            completed, arguments = self._run(root, root / "missing")

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(arguments, ["-N", "fixture", "worker.sh", "--flag"])
            self.assertIn("notifications disabled", completed.stderr)

    def test_non_file_path_submits_normally(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)

            completed, arguments = self._run(root, root)

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(arguments, ["-N", "fixture", "worker.sh", "--flag"])

    def test_unreadable_file_submits_normally(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            email_file = root / "email"
            email_file.write_text("person@example.org\n")
            email_file.chmod(0)
            if os.access(email_file, os.R_OK):
                email_file.chmod(0o600)
                self.skipTest("test user can read mode-000 files")
            try:
                completed, arguments = self._run(root, email_file)
            finally:
                email_file.chmod(0o600)

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(arguments, ["-N", "fixture", "worker.sh", "--flag"])

    def test_malformed_or_multiline_file_submits_normally(self):
        for content in ("not-an-address\n", "first@example.org\nsecond@example.org\n"):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                email_file = root / "email"
                email_file.write_text(content)

                completed, arguments = self._run(root, email_file)

                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertEqual(arguments, ["-N", "fixture", "worker.sh", "--flag"])
                self.assertIn("notifications disabled", completed.stderr)

    def test_qsub_exit_status_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            email_file = root / "email"
            email_file.write_text("person@example.org\n")

            completed, arguments = self._run(root, email_file, qsub_exit=23)

            self.assertEqual(completed.returncode, 23)
            self.assertEqual(arguments[:4], ["-m", "beas", "-M", "person@example.org"])


if __name__ == "__main__":
    unittest.main()
