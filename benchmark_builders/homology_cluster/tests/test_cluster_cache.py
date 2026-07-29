from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.cluster_cache import (
    ALIGNMENT_STATISTICS_FILE,
    ASSIGNMENTS_FILE,
    CHECKPOINT_ASSIGNMENTS_FILE,
    CHECKPOINT_MARKER,
    CACHE_ROOT_MARKER,
    CLUSTER_FASTA_FILE,
    cluster_cache_contract,
    cluster_cache_directory,
    initialize_cluster_cache_root,
    import_publication_cluster_cache,
    inspect_cluster_cache,
    inspect_cluster_cache_root,
    inspect_cluster_checkpoint,
    load_cluster_checkpoint,
    load_cluster_cache,
    publish_cluster_checkpoint,
    publish_cluster_cache,
)
from homology_cluster_benchmark.common_cache import build_common_preprocessing_cache
from homology_cluster_benchmark.frozen_inputs import write_synthetic_fixture_manifest
from homology_cluster_benchmark.inputs import resolve_input, sha256_file
from homology_cluster_benchmark.mmseqs import (
    ClusterIndex,
    MMseqsRuntime,
    write_command_manifest,
    build_mmseqs_commands,
)
from homology_cluster_benchmark.pipeline import _input_specs, build_benchmark
from homology_cluster_benchmark.uniref import UniRefIndex

from tests.helpers import FIXTURES, fixture_config, uniref50_fixture_config


VERSION = "18-8cc5c"


def _fake_mmseqs(path: Path) -> None:
    path.write_text(f"#!/bin/sh\nprintf '%s\\n' '{VERSION}'\n", encoding="utf-8")
    path.chmod(0o755)


class ClusterCacheTests(unittest.TestCase):
    def _runtime(self, executable: Path) -> MMseqsRuntime:
        return MMseqsRuntime(
            requested_executable=str(executable),
            resolved_executable=str(executable),
            observed_version=VERSION,
            version_token=VERSION,
            version_exit_code=0,
            executable_sha256=sha256_file(executable),
        )

    def test_uniref50_contract_has_a_distinct_cache_namespace(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "mmseqs"
            _fake_mmseqs(executable)
            runtime = self._runtime(executable)
            legacy = fixture_config(
                root / "legacy-output", root / "legacy-temp",
                mmseqs_bin=str(executable), expected_mmseqs_version=VERSION,
            )
            uniref50 = uniref50_fixture_config(
                root / "u50-output", root / "u50-temp",
                mmseqs_bin=str(executable), expected_mmseqs_version=VERSION,
            )
            legacy_contract = cluster_cache_contract(
                legacy, runtime, sha256_file(FIXTURES / "uniref90.fasta")
            )
            u50_contract = cluster_cache_contract(
                uniref50, runtime, sha256_file(FIXTURES / "uniref50.fasta")
            )
            self.assertIn("uniref90_2026_02", str(cluster_cache_directory(root, legacy_contract)))
            self.assertIn("uniref50_2026_02", str(cluster_cache_directory(root, u50_contract)))
            self.assertNotEqual(legacy_contract, u50_contract)

    def test_publish_load_and_tamper_detection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "mmseqs"
            _fake_mmseqs(executable)
            config = fixture_config(
                root / "output",
                root / "temp",
                mmseqs_bin=str(executable),
                expected_mmseqs_version=VERSION,
            )
            uniref = UniRefIndex.build(
                FIXTURES / "uniref90.fasta", root / "uniref.sqlite"
            )
            clusters = ClusterIndex.build(
                FIXTURES / "clusters.tsv", uniref, root / "clusters.sqlite"
            )
            commands = root / "mmseqs_commands.tsv"
            write_command_manifest(
                commands,
                build_mmseqs_commands(
                    config, FIXTURES / "uniref90.fasta", root / "mmseqs-work"
                ),
            )
            contract = cluster_cache_contract(
                config,
                self._runtime(executable),
                sha256_file(FIXTURES / "uniref90.fasta"),
            )
            cache = publish_cluster_cache(
                root / "cache",
                contract,
                clusters,
                commands,
                producer={"run_id": "fixture"},
            )
            self.assertEqual(cache.payload["counts"], {"members": 7, "clusters": 6})
            self.assertTrue((root / "cache" / CACHE_ROOT_MARKER).is_file())
            self.assertEqual(load_cluster_cache(root / "cache", contract).root, cache.root)
            assignments = cache.root / ASSIGNMENTS_FILE
            assignments.write_bytes(assignments.read_bytes() + b"tamper")
            with self.assertRaisesRegex(ValueError, "file-size mismatch"):
                inspect_cluster_cache(cache.root, verify_file_hashes=True)

    def test_post_createtsv_checkpoint_is_atomic_resumable_and_tamper_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "mmseqs"
            _fake_mmseqs(executable)
            config = fixture_config(
                root / "output",
                root / "temp",
                mmseqs_bin=str(executable),
                expected_mmseqs_version=VERSION,
            )
            commands = root / "mmseqs_commands.tsv"
            write_command_manifest(
                commands,
                build_mmseqs_commands(
                    config, FIXTURES / "uniref90.fasta", root / "mmseqs-work"
                ),
            )
            contract = cluster_cache_contract(
                config,
                self._runtime(executable),
                sha256_file(FIXTURES / "uniref90.fasta"),
            )
            checkpoint = publish_cluster_checkpoint(
                root / "cache",
                contract,
                FIXTURES / "clusters.tsv",
                commands,
                producer={"run_id": "fixture"},
            )
            self.assertTrue((checkpoint.root / CHECKPOINT_MARKER).is_file())
            self.assertEqual(
                load_cluster_checkpoint(root / "cache", contract).root,
                checkpoint.root,
            )
            inspect_cluster_checkpoint(checkpoint.root, verify_file_hashes=True)
            assignments = checkpoint.root / CHECKPOINT_ASSIGNMENTS_FILE
            assignments.write_bytes(assignments.read_bytes() + b"tamper")
            with self.assertRaisesRegex(ValueError, "file-size mismatch"):
                inspect_cluster_checkpoint(checkpoint.root, verify_file_hashes=True)

    def test_profile_locked_exports_are_published_and_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "mmseqs"
            _fake_mmseqs(executable)
            config = fixture_config(
                root / "output",
                root / "temp",
                mmseqs_bin=str(executable),
                expected_mmseqs_version=VERSION,
                mmseqs_profile="daniel-aligned-defaults",
                createdb_shuffle=None,
                cluster_reassign=0,
                evalue=None,
                export_alignment_statistics=True,
                export_cluster_fasta=True,
                cluster_assignments=None,
                cluster_cache_root=root / "cache",
            )
            uniref = UniRefIndex.build(
                FIXTURES / "uniref90.fasta", root / "uniref.sqlite"
            )
            clusters = ClusterIndex.build(
                FIXTURES / "clusters.tsv", uniref, root / "clusters.sqlite"
            )
            commands = root / "mmseqs_commands.tsv"
            write_command_manifest(
                commands,
                build_mmseqs_commands(
                    config, FIXTURES / "uniref90.fasta", root / "mmseqs-work"
                ),
            )
            statistics = root / "statistics.tsv"
            statistics.write_text("query\ttarget\tevalue\tbits\traw\tpident\nU1\tU1\t0\t1\t1\t100\n")
            fasta = root / "clusters.faa"
            fasta.write_text(">U1\nAAAA\n")
            cache = publish_cluster_cache(
                root / "cache",
                cluster_cache_contract(
                    config, self._runtime(executable),
                    sha256_file(FIXTURES / "uniref90.fasta"),
                ),
                clusters,
                commands,
                producer={"run_id": "fixture"},
                derived_artifacts={
                    "alignment_statistics": statistics,
                    "cluster_fasta": fasta,
                },
            )
            self.assertTrue((cache.root / ALIGNMENT_STATISTICS_FILE).is_file())
            self.assertTrue((cache.root / CLUSTER_FASTA_FILE).is_file())
            self.assertEqual(
                set(cache.payload["exports"]),
                {"alignment_statistics", "cluster_fasta"},
            )
            (cache.root / CLUSTER_FASTA_FILE).write_text(">U1\nTAMPER\n")
            with self.assertRaisesRegex(ValueError, "file-size mismatch|file hash mismatch"):
                inspect_cluster_cache(cache.root, verify_file_hashes=True)

    def test_contract_excludes_downstream_and_operational_choices(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "mmseqs"
            _fake_mmseqs(executable)
            base = fixture_config(
                root / "output",
                root / "temp",
                mmseqs_bin=str(executable),
                expected_mmseqs_version=VERSION,
            )
            changed = replace(
                base,
                split_policy="cluster-count-random",
                training_population="all-cluster-members",
                seed=987,
                threads=6,
                requested_slots=6,
                allocated_slots=6,
                framework_revision="a" * 40,
            )
            runtime = self._runtime(executable)
            digest = sha256_file(FIXTURES / "uniref90.fasta")
            self.assertEqual(
                cluster_cache_contract(base, runtime, digest),
                cluster_cache_contract(changed, runtime, digest),
            )
            different_identity = replace(base, identity=0.20)
            self.assertNotEqual(
                cluster_cache_directory(
                    root / "cache", cluster_cache_contract(base, runtime, digest)
                ),
                cluster_cache_directory(
                    root / "cache",
                    cluster_cache_contract(different_identity, runtime, digest),
                ),
            )

    def test_root_marker_must_match_exactly(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "cache"
            marker = initialize_cluster_cache_root(root)
            self.assertEqual(marker, root.resolve() / CACHE_ROOT_MARKER)
            inspect_cluster_cache_root(root)
            payload = json.loads(marker.read_text(encoding="utf-8"))
            payload["role"] = "wrong"
            marker.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "incompatible"):
                inspect_cluster_cache_root(root)

    def test_second_pipeline_run_reuses_cache_without_mmseqs_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "mmseqs"
            _fake_mmseqs(executable)
            cache_root = root / "cluster-cache"
            base = fixture_config(
                root / "first-output",
                root / "first-temp",
                cluster_assignments=None,
                cluster_cache_root=cache_root,
                mmseqs_bin=str(executable),
                expected_mmseqs_version=VERSION,
            )

            def execute_fixture(commands, log_dir):
                log_dir.mkdir(parents=True, exist_ok=True)
                (log_dir / "mmseqs_fixture.log").write_text("fixture\n", encoding="utf-8")
                target = next(
                    Path(command.argv[-1])
                    for command in commands if command.stage == "createtsv"
                )
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(FIXTURES / "clusters.tsv", target)

            with mock.patch(
                "homology_cluster_benchmark.pipeline.execute_commands",
                side_effect=execute_fixture,
            ):
                first = build_benchmark(base)

            second_config = replace(
                base,
                output_dir=root / "second-output",
                temp_dir=root / "second-temp",
                require_cluster_cache=True,
                threads=6,
                requested_slots=6,
                allocated_slots=6,
            )
            with mock.patch(
                "homology_cluster_benchmark.pipeline.execute_commands",
                side_effect=AssertionError("MMseqs must not run on a cache hit"),
            ):
                second = build_benchmark(second_config)

            scientific_names = {
                *(f"{aspect}-{split}.csv" for aspect in ("bp", "cc", "mf")
                  for split in ("training", "validation", "test")),
                "cluster_split_assignments.tsv",
                "protein_cluster_assignments.tsv",
                "retained_clusters.tsv",
            }
            self.assertEqual(
                {name: sha256_file(first.output_dir / name) for name in scientific_names},
                {name: sha256_file(second.output_dir / name) for name in scientific_names},
            )
            manifest = json.loads(
                (second.output_dir / "input_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["cluster_cache"]["action"], "reused")
            self.assertTrue((second.output_dir / "cluster_cache_manifest.json").is_file())
            self.assertEqual(
                manifest["cluster_cache"]["assignment_sha256"],
                json.loads(
                    (second.output_dir / "cluster_cache_manifest.json").read_text(
                        encoding="utf-8"
                    )
                )["assignment_sha256"],
            )

    def test_pipeline_resumes_checkpoint_after_post_createtsv_validation_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "mmseqs"
            _fake_mmseqs(executable)
            cache_root = root / "cluster-cache"
            first = fixture_config(
                root / "failed-output",
                root / "failed-temp",
                cluster_assignments=None,
                cluster_cache_root=cache_root,
                mmseqs_bin=str(executable),
                expected_mmseqs_version=VERSION,
            )

            def execute_fixture(commands, log_dir):
                log_dir.mkdir(parents=True, exist_ok=True)
                (log_dir / "mmseqs_fixture.log").write_text(
                    "fixture\n", encoding="utf-8"
                )
                target = next(
                    Path(command.argv[-1])
                    for command in commands if command.stage == "createtsv"
                )
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(FIXTURES / "clusters.tsv", target)

            with mock.patch(
                "homology_cluster_benchmark.pipeline.execute_commands",
                side_effect=execute_fixture,
            ), mock.patch(
                "homology_cluster_benchmark.pipeline.ClusterIndex.build",
                side_effect=RuntimeError("synthetic validation interruption"),
            ), self.assertRaisesRegex(RuntimeError, "synthetic validation interruption"):
                build_benchmark(first)

            checkpoints = list(cache_root.rglob(CHECKPOINT_MARKER))
            self.assertEqual(len(checkpoints), 1)
            inspect_cluster_checkpoint(
                checkpoints[0].parent, verify_file_hashes=True
            )

            resumed = replace(
                first,
                output_dir=root / "resumed-output",
                temp_dir=root / "resumed-temp",
            )
            with mock.patch(
                "homology_cluster_benchmark.pipeline.execute_commands",
                side_effect=AssertionError("MMseqs clustering must not repeat"),
            ):
                result = build_benchmark(resumed)

            self.assertTrue((result.output_dir / "RUN_COMPLETE.json").is_file())
            self.assertFalse(list(cache_root.rglob(CHECKPOINT_MARKER)))
            self.assertEqual(
                json.loads(
                    (result.output_dir / "input_manifest.json").read_text(
                        encoding="utf-8"
                    )
                )["cluster_cache"]["action"],
                "built",
            )

    def test_completed_publication_can_be_imported_without_mmseqs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = fixture_config(root / "publication", root / "publication-work")
            specs = _input_specs(config)
            resolved = {
                name: resolve_input(spec, root / "downloads", allow_downloads=False)
                for name, spec in specs.items()
            }
            manifest = write_synthetic_fixture_manifest(
                root / "frozen-inputs.json",
                specs,
                resolved,
                config.uniprot_source_scope,
            )
            common = build_common_preprocessing_cache(
                root / "common-cache",
                root / "common-work",
                manifest.path,
                {name: item.resolved_path for name, item in resolved.items()},
                source_scope=config.uniprot_source_scope,
                fixture_mode=True,
            )
            publication = build_benchmark(
                replace(config, frozen_input_manifest=manifest.path)
            ).output_dir
            metadata_path = publication / "publication_metadata.json"
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata.update({
                "fixture_mode": False,
                "expected_mmseqs_version": VERSION,
                "observed_mmseqs_version": VERSION,
                "mmseqs_executable_sha256": "a" * 64,
            })
            metadata_path.write_text(
                json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            with mock.patch(
                "homology_cluster_benchmark.pipeline.validate_publication"
            ):
                imported = import_publication_cluster_cache(
                    publication,
                    common,
                    root / "cluster-cache",
                    root / "import-work",
                )
            self.assertEqual(imported.payload["counts"], {"members": 7, "clusters": 6})
            self.assertEqual(
                imported.payload["producer"]["imported_from_publication"]["run_id"],
                metadata["run_id"],
            )
            ClusterIndex.build(
                imported.assignments,
                UniRefIndex(common / "uniref90.sqlite"),
                root / "imported-clusters.sqlite",
                has_header=True,
            )


if __name__ == "__main__":
    unittest.main()
