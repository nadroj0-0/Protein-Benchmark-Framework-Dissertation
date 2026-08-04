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
DANIEL_ARRAY = (
    WORKSPACE_ROOT / "hpc_jobs" / "active"
    / "hpc_homology_cluster_runtime_array_24core_daniel_aligned.sh"
)
UNIREF50_ARRAY = (
    WORKSPACE_ROOT / "hpc_jobs" / "active"
    / "hpc_homology_cluster_runtime_array_24core_uniref50.sh"
)
UNIREF50_ARRAY_12CORE = (
    WORKSPACE_ROOT / "hpc_jobs" / "active"
    / "hpc_homology_cluster_runtime_array_12core_uniref50.sh"
)
CACHED_RANDOM_RESPLIT = (
    WORKSPACE_ROOT / "hpc_jobs" / "active"
    / "hpc_homology_cluster_cached_resplit_30_25_20.sh"
)
UNIREF50_COMMON_CACHE = (
    WORKSPACE_ROOT / "hpc_jobs" / "active"
    / "hpc_homology_uniref50_common_cache.sh"
)
SUPERVISOR_DANIEL_30 = (
    WORKSPACE_ROOT / "hpc_jobs" / "active"
    / "hpc_homology_supervisor_daniel_30.sh"
)
GUARDED_WORKER = (
    WORKSPACE_ROOT / "hpc_jobs" / "active" / "hpc_homology_cluster_benchmark.sh"
)
SNAPSHOT = WORKSPACE_ROOT / "hpc_jobs" / "active" / "hpc_homology_progress_snapshot.sh"


class RuntimeHPCEntrypointTests(unittest.TestCase):
    def test_cached_random_resplit_is_small_and_cannot_recluster(self):
        worker = CACHED_RANDOM_RESPLIT.read_text()
        self.assertIn("#$ -pe smp 4", worker)
        self.assertIn("#$ -t 1-3", worker)
        self.assertIn("#$ -tc 3", worker)
        self.assertIn("export SPLIT_POLICY=cluster-count-random", worker)
        self.assertIn("export REQUIRE_HOMOLOGY_COMMON_CACHE=1", worker)
        self.assertIn("export REQUIRE_HOMOLOGY_CLUSTER_CACHE=1", worker)
        self.assertIn("uniref50_sensitivity_4", worker)
        self.assertIn("cached_random_resplit", worker)
        self.assertNotIn("qsub", worker)

    def test_uniref50_jobs_are_namespaced_and_block_on_the_new_shared_cache(self):
        array = UNIREF50_ARRAY.read_text()
        self.assertIn("#$ -pe smp 24", array)
        self.assertIn("export UNIREF_LEVEL=50", array)
        self.assertIn('MMSEQS_SENSITIVITY="${MMSEQS_SENSITIVITY:-4}"', array)
        self.assertIn("uniref50/2026_02/uniref50.fasta.gz", array)
        self.assertIn("uniref50/common_preprocessing", array)
        self.assertIn("export REQUIRE_HOMOLOGY_COMMON_CACHE=1", array)
        self.assertNotIn("qsub", array)

        array_12core = UNIREF50_ARRAY_12CORE.read_text()
        self.assertIn("#$ -pe smp 12", array_12core)
        self.assertIn("#$ -t 1-6", array_12core)
        self.assertIn("#$ -tc 6", array_12core)
        self.assertIn("export UNIREF_LEVEL=50", array_12core)
        self.assertIn('MMSEQS_SENSITIVITY="${MMSEQS_SENSITIVITY:-4}"', array_12core)
        self.assertIn("export MMSEQS_PROFILE=daniel-aligned-defaults", array_12core)
        self.assertIn("uniref50_sensitivity_4_daniel_aligned_12core", array_12core)
        self.assertIn('export MINIMUM_SCRATCH_GB="${MINIMUM_SCRATCH_GB:-300}"', array_12core)
        self.assertNotIn("qsub", array_12core)

        cache = UNIREF50_COMMON_CACHE.read_text()
        self.assertIn("--uniref-level 50", cache)
        self.assertIn("--uniref50-fasta", cache)
        self.assertIn("--full-hashes", cache)
        self.assertNotIn("qsub", cache)

    def test_daniel_aligned_array_requests_24_cores_and_locked_profile(self):
        worker = DANIEL_ARRAY.read_text()
        self.assertIn("#$ -pe smp 24", worker)
        self.assertIn("#$ -t 1-6", worker)
        self.assertIn("#$ -tc 6", worker)
        self.assertIn("export MMSEQS_PROFILE=daniel-aligned-defaults", worker)
        self.assertIn("export REQUIRE_HOMOLOGY_COMMON_CACHE=1", worker)
        self.assertIn("MINIMUM_CLUSTER_CACHE_FREE_GB", worker)
        self.assertIn("/SAN/bioinf/bmpfp/benchmarks/homology/2026_02", worker)
        self.assertIn("/SAN/bioinf/bmpfp/derived_inputs/homology/2026_02", worker)
        self.assertIn('export MINIMUM_SCRATCH_GB="${MINIMUM_SCRATCH_GB:-400}"', worker)
        self.assertNotIn("qsub", worker)

    def test_supervisor_external_job_is_provenance_isolated(self):
        worker = SUPERVISOR_DANIEL_30.read_text()
        self.assertIn("#$ -t 1\n", worker)
        self.assertIn("export HOMOLOGY_RUNTIME_KIND=array", worker)
        self.assertNotIn("export HOMOLOGY_RUNTIME_KIND=pilot", worker)
        self.assertIn("export UNIREF_LEVEL=50", worker)
        self.assertIn("export MMSEQS_SENSITIVITY=4", worker)
        self.assertIn("export SPLIT_POLICY=cluster-count-random", worker)
        self.assertIn("supervisor_daniel_buchan", worker)
        self.assertIn("EXTERNAL_CLUSTER_ASSIGNMENTS", worker)
        self.assertIn("EXTERNAL_CLUSTER_PROVENANCE", worker)
        self.assertIn("unset HOMOLOGY_CLUSTER_CACHE_ROOT", worker)
        self.assertNotIn("qsub", worker)

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

    def test_pilot_accepts_symlinked_scratch_base_and_only_deletes_owned_work(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, original_scratch = self._environment(root, "pilot", "1")
            resolved_scratch = root / "resolved-scratch"
            resolved_scratch.mkdir()
            scratch_link = root / "scratch-link"
            scratch_link.symlink_to(resolved_scratch, target_is_directory=True)
            original_scratch.rmdir()
            env["WORK_BASE"] = str(scratch_link)

            completed = self._run(env)

            self.assertEqual(completed.returncode, 0, completed.stdout)
            self.assertTrue(scratch_link.is_symlink())
            self.assertFalse(list(resolved_scratch.iterdir()))
            self.assertEqual(len(list((root / "results").rglob("benchmark/result.txt"))), 1)

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

    def test_runtime_driver_namespaces_uniref50_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env, scratch = self._environment(root, "array", "1")
            env["UNIREF_LEVEL"] = "50"
            env["MMSEQS_SENSITIVITY"] = "4"
            completed = self._run(env)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            summary = next((root / "results").rglob("FINAL_RESULT_PATH.txt"))
            self.assertIn("uniref50_sensitivity_4", summary.read_text())
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
        self.assertIn("#$ -pe smp 4", pilot)
        self.assertIn("#$ -l tmem=16G", pilot)
        self.assertIn("#$ -l tscratch=75G", pilot)
        self.assertIn("#$ -l scratch0free=300G", pilot)
        self.assertIn("#$ -l h_rt=96:0:0", pilot)
        self.assertIn("#$ -t 1-6\n", array)
        self.assertIn("#$ -tc 6\n", array)
        self.assertIn("HOMOLOGY_RUNTIME_KIND=array", array)
        self.assertIn("#$ -pe smp 8", array)
        self.assertIn("#$ -l tmem=15G", array)
        self.assertIn("#$ -l tscratch=38G", array)
        self.assertIn("#$ -l h_rt=168:0:0", array)
        for text in (pilot, array):
            self.assertIn("#$ -l scratch0free=300G", text)
            self.assertIn("run_homology_cluster_runtime_hpc.sh", text)
            self.assertNotIn("wget ", text)
            self.assertNotIn("rm -rf", text)

        guarded_worker = GUARDED_WORKER.read_text()
        self.assertIn("#$ -pe smp 2", guarded_worker)
        self.assertIn("#$ -l tmem=84G", guarded_worker)
        self.assertIn("#$ -l tscratch=150G", guarded_worker)
        self.assertIn("#$ -l scratch0free=300G", guarded_worker)
        self.assertIn("#$ -l h_rt=168:0:0", guarded_worker)

    def test_runtime_driver_supports_the_clusters_legacy_git(self):
        driver = DRIVER.read_text()
        self.assertIn("git_in_dir()", driver)
        self.assertNotIn("git -C", driver)

    def test_progress_snapshot_is_read_only_and_host_guarded(self):
        snapshot = SNAPSHOT.read_text()
        self.assertIn('--expected-host', snapshot)
        self.assertIn('--main-job-id', snapshot)
        self.assertIn('copy_log', snapshot)
        self.assertIn('stable_during_copy', snapshot)
        self.assertIn('homology_runtime_${MAIN_JOB_ID}_${task_id}_${identity}', snapshot)
        self.assertNotIn('rm -rf', snapshot)
        self.assertNotIn('qdel', snapshot)

    def test_progress_snapshot_copies_logs_without_modifying_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            scratch = root / "scratch"
            work = scratch / "homology_runtime_123_1_30_runtime-123"
            artifact_logs = work / "artifacts" / "logs"
            mmseqs_logs = work / "tmp" / "homology-fixture" / "logs" / "mmseqs"
            artifact_logs.mkdir(parents=True)
            mmseqs_logs.mkdir(parents=True)
            runtime_log = artifact_logs / "runtime.log"
            cluster_log = mmseqs_logs / "mmseqs_cluster.log"
            runtime_log.write_text("runtime-stage\n")
            cluster_log.write_text("cluster-stage\n")
            destination = root / "snapshots"
            environment = os.environ.copy()
            environment["HOMOLOGY_SNAPSHOT_SCRATCH_BASE"] = str(scratch)
            host = subprocess.check_output(["hostname", "-s"], text=True).strip()

            completed = subprocess.run(
                [
                    "bash",
                    str(SNAPSHOT),
                    "--main-job-id",
                    "123",
                    "--expected-host",
                    host,
                    "--task",
                    "1:30",
                    "--destination-root",
                    str(destination),
                ],
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stdout)
            latest = destination / "job_123" / f"host_{host}" / "latest"
            self.assertTrue(latest.is_symlink())
            self.assertEqual(
                (latest / "task_1_identity_30" / "logs" / "runtime.log").read_text(),
                "runtime-stage\n",
            )
            self.assertEqual(
                (
                    latest
                    / "task_1_identity_30"
                    / "logs"
                    / "mmseqs"
                    / "mmseqs_cluster.log"
                ).read_text(),
                "cluster-stage\n",
            )
            self.assertEqual(runtime_log.read_text(), "runtime-stage\n")
            self.assertEqual(cluster_log.read_text(), "cluster-stage\n")

    def test_runtime_driver_exports_host_verified_git_state_for_singularity(self):
        driver = DRIVER.read_text()
        self.assertIn(
            'export HOMOLOGY_HOST_GIT_VERIFIED_COMMIT="$FRAMEWORK_REVISION"', driver
        )
        self.assertIn('export HOMOLOGY_HOST_GIT_VERIFIED_CLEAN=1', driver)
        self.assertIn(
            'export HOMOLOGY_HOST_GIT_VERIFIED_REPOSITORY="$FRAMEWORK_DIR"', driver
        )
        self.assertIn(
            "SINGULARITYENV_HOMOLOGY_HOST_GIT_VERIFIED_COMMIT", driver
        )

    def test_runtime_driver_separates_release_tag_from_binary_identity(self):
        driver = DRIVER.read_text()
        self.assertIn(
            'MMSEQS_RELEASE_TAG="${MMSEQS_RELEASE_TAG:-18-8cc5c}"', driver
        )
        self.assertIn(
            "EXPECTED_MMSEQS_BINARY_VERSION=\"${EXPECTED_MMSEQS_BINARY_VERSION:-"
            "8cc5ce367b5638c4306c2d7cfc652dd099a4643f}\"",
            driver,
        )
        self.assertIn('echo "release_tag=$MMSEQS_RELEASE_TAG"', driver)
        self.assertIn(
            'echo "expected_binary_version=$EXPECTED_MMSEQS_BINARY_VERSION"',
            driver,
        )
        self.assertNotIn(
            '[[ "$observed_mmseqs_version" == "$EXPECTED_MMSEQS_VERSION" ]]',
            driver,
        )

    def test_runtime_driver_freezes_goa_234_to_immutable_archive(self):
        driver = DRIVER.read_text()
        self.assertIn(
            'GOA_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/'
            'goa_uniprot_all.gaf.234.gz"',
            driver,
        )
        self.assertIn(
            'PINNED_GOA_SHA256="f315375b07946a0649142b2f4de2e15e282316989677a04e7a561203186dd2ff"',
            driver,
        )
        self.assertNotIn("GOA_RELEASES_URL", driver)
        self.assertNotIn("GOA_MD5_URL", driver)
        self.assertNotIn("goa_current_release_numbers", driver)
        self.assertIn(
            'artifact_catalog_configure "$FRAMEWORK_DIR" "${ARTIFACT_CATALOG:-}"',
            driver,
        )
        self.assertIn("use_catalog_input_if_available GOA GOA_SHA256", driver)
        self.assertIn('unset "$path_variable" "$hash_variable"', driver)
        self.assertIn("    goa_t1 \\", driver)
        self.assertNotIn("SAN_INPUT_ROOT=", driver)

    def test_runtime_driver_prefers_portable_common_cache_and_keeps_raw_fallback(self):
        driver = DRIVER.read_text()
        self.assertIn("homology_common_preprocessing_2026_02", driver)
        self.assertIn("HOMOLOGY_COMMON_PREPROCESSING_CACHE", driver)
        self.assertIn("REQUIRE_HOMOLOGY_COMMON_CACHE", driver)
        self.assertIn("Staging common preprocessing cache into job-owned scratch", driver)
        self.assertIn(
            'if [[ -z "$HOMOLOGY_COMMON_PREPROCESSING_CACHE" ]]; then', driver
        )
        self.assertIn("stage_or_download idmapping", driver)
        self.assertIn("homology_cluster_benchmark.runtime_contract policy", driver)
        self.assertIn(
            'HOMOLOGY_COMMON_PREPROCESSING_CACHE="$STAGED_COMMON_CACHE"', driver
        )
        self.assertNotIn("/SAN/bioinf/bmpfp", driver)

    def test_runtime_driver_uses_persistent_cluster_cache_without_scratch_copy(self):
        driver = DRIVER.read_text()
        self.assertIn("homology_mmseqs_cluster_cache_root_2026_02", driver)
        self.assertIn("HOMOLOGY_CLUSTER_CACHE_ROOT", driver)
        self.assertIn("REQUIRE_HOMOLOGY_CLUSTER_CACHE", driver)
        self.assertIn(
            'HOMOLOGY_CLUSTER_CACHE_ROOT="$HOMOLOGY_CLUSTER_CACHE_ROOT"', driver
        )
        self.assertNotIn("STAGED_CLUSTER_CACHE", driver)
        self.assertNotIn("cp -a \"$HOMOLOGY_CLUSTER_CACHE_ROOT\"", driver)

    def test_runtime_driver_allows_validated_external_assignments_without_cluster_cache(self):
        driver = DRIVER.read_text()
        self.assertIn("EXTERNAL_CLUSTER_MODE=0", driver)
        self.assertIn("EXTERNAL_CLUSTER_MODE=1", driver)
        self.assertIn(
            "Using provenance-paired external cluster assignments; "
            "framework cluster cache is intentionally disabled",
            driver,
        )
        self.assertIn('&& "$EXTERNAL_CLUSTER_MODE" != "1"', driver)

    def test_runtime_driver_binds_and_probes_persistent_cache_before_initializing_it(self):
        driver = DRIVER.read_text()
        bind_call = 'add_mmfp_singularity_bind "$bind_root"'
        probe_message = "Cluster-cache container read/write preflight passed"
        initialize_call = "homology_cluster_benchmark.cluster_cache init-root"

        self.assertIn("configure_cluster_cache_container_access()", driver)
        self.assertIn(bind_call, driver)
        self.assertIn(probe_message, driver)
        self.assertIn("HOMOLOGY_CLUSTER_CACHE_PREFLIGHT_ONLY", driver)
        self.assertIn("CLUSTER_CACHE_PREFLIGHT_COMPLETE.txt", driver)
        self.assertLess(driver.index(bind_call), driver.index(initialize_call))
        self.assertLess(driver.index(probe_message), driver.index(initialize_call))


if __name__ == "__main__":
    unittest.main()
