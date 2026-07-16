from __future__ import annotations

import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from homology_cluster_benchmark.provenance import git_state


COMMIT = "a" * 40


class GitStateTests(unittest.TestCase):
    def _verified_environment(self, repository: Path) -> dict[str, str]:
        return {
            "HOMOLOGY_HOST_GIT_VERIFIED_COMMIT": COMMIT,
            "HOMOLOGY_HOST_GIT_VERIFIED_CLEAN": "1",
            "HOMOLOGY_HOST_GIT_VERIFIED_REPOSITORY": str(repository),
        }

    def test_host_verified_state_does_not_execute_git(self):
        with tempfile.TemporaryDirectory() as tmp:
            repository = Path(tmp).resolve()
            with mock.patch.dict(os.environ, self._verified_environment(repository), clear=False):
                with mock.patch("homology_cluster_benchmark.provenance.subprocess.run") as run:
                    observed = git_state(repository)

            run.assert_not_called()
            self.assertEqual(observed["commit"], COMMIT)
            self.assertIs(observed["dirty"], False)
            self.assertEqual(observed["status_porcelain"], [])
            self.assertEqual(observed["verification"], "hpc-wrapper-host-git")

    def test_partial_or_unclean_host_verification_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            repository = Path(tmp).resolve()
            cases = (
                {"HOMOLOGY_HOST_GIT_VERIFIED_COMMIT": COMMIT},
                {
                    **self._verified_environment(repository),
                    "HOMOLOGY_HOST_GIT_VERIFIED_CLEAN": "0",
                },
                {
                    **self._verified_environment(repository),
                    "HOMOLOGY_HOST_GIT_VERIFIED_COMMIT": "not-a-commit",
                },
                {
                    **self._verified_environment(repository),
                    "HOMOLOGY_HOST_GIT_VERIFIED_REPOSITORY": str(repository / "other"),
                },
            )
            for environment in cases:
                with self.subTest(environment=environment):
                    with mock.patch.dict(os.environ, environment, clear=True):
                        with self.assertRaises(ValueError):
                            git_state(repository)

    def test_direct_runs_still_execute_git(self):
        with tempfile.TemporaryDirectory() as tmp:
            repository = Path(tmp).resolve()
            commit = mock.Mock(returncode=0, stdout=f"{COMMIT}\n")
            status = mock.Mock(returncode=0, stdout="")
            with mock.patch.dict(os.environ, {}, clear=True):
                with mock.patch(
                    "homology_cluster_benchmark.provenance.subprocess.run",
                    side_effect=(commit, status),
                ) as run:
                    observed = git_state(repository)

            self.assertEqual(run.call_count, 2)
            self.assertEqual(observed["commit"], COMMIT)
            self.assertIs(observed["dirty"], False)
            self.assertEqual(observed["verification"], "git-executable")


if __name__ == "__main__":
    unittest.main()
