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
    ASSIGNMENTS_FILE,
    CACHE_ROOT_MARKER,
    cluster_cache_contract,
    cluster_cache_directory,
    initialize_cluster_cache_root,
    import_publication_cluster_cache,
    inspect_cluster_cache,
    inspect_cluster_cache_root,
    load_cluster_cache,
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

from tests.helpers import FIXTURES, fixture_config


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
