from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.common_cache import (
    CACHE_MARKER,
    CACHE_SCHEMA_VERSION,
    CommonPreprocessingState,
    SCHEMA_V2_PREPROCESSING_SOURCE_SHA256,
    _load_common_preprocessing_state,
    build_common_preprocessing_cache,
    inspect_common_preprocessing_cache,
)
from homology_cluster_benchmark.frozen_inputs import write_synthetic_fixture_manifest
from homology_cluster_benchmark.inputs import resolve_input, sha256_file
from homology_cluster_benchmark.pipeline import _input_specs, build_benchmark

from tests.helpers import fixture_config, uniref50_fixture_config


class CommonPreprocessingCacheTests(unittest.TestCase):
    def test_uniref50_common_cache_is_distinct_and_reusable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = uniref50_fixture_config(root / "unused-output", root / "unused-temp")
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
                uniref_level=50,
            )
            cache = build_common_preprocessing_cache(
                root / "common-cache",
                root / "cache-work",
                manifest.path,
                {name: spec.path for name, spec in specs.items() if spec.path is not None},
                source_scope=config.uniprot_source_scope,
                fixture_mode=True,
                uniref_level=50,
            )
            payload = inspect_common_preprocessing_cache(
                cache, expected_uniref_level=50, verify_file_hashes=True
            )
            self.assertEqual(payload["policy"]["uniref_level"], 50)
            self.assertEqual(payload["counts"]["uniref50_entries"], 6)
            result = build_benchmark(replace(
                config,
                output_dir=root / "cached-output",
                temp_dir=root / "cached-temp",
                frozen_input_manifest=manifest.path,
                common_preprocessing_cache=cache,
            ))
            self.assertTrue((result.output_dir / "RUN_COMPLETE.json").is_file())

    def _fixture_manifest(self, root: Path):
        config = fixture_config(root / "unused-output", root / "unused-temp")
        specs = _input_specs(config)
        resolved = {
            name: resolve_input(spec, root / "downloads", allow_downloads=False)
            for name, spec in specs.items()
        }
        return config, write_synthetic_fixture_manifest(
            root / "frozen-inputs.json", specs, resolved, config.uniprot_source_scope
        )

    def test_cached_and_raw_preprocessing_produce_identical_scientific_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config, manifest = self._fixture_manifest(root)
            specs = _input_specs(config)
            cache = build_common_preprocessing_cache(
                root / "common-cache",
                root / "cache-work",
                manifest.path,
                {name: spec.path for name, spec in specs.items() if spec.path is not None},
                source_scope=config.uniprot_source_scope,
                fixture_mode=True,
            )

            raw = build_benchmark(
                replace(
                    config,
                    output_dir=root / "raw-output",
                    temp_dir=root / "raw-temp",
                    frozen_input_manifest=manifest.path,
                )
            )
            cached = build_benchmark(
                replace(
                    config,
                    output_dir=root / "cached-output",
                    temp_dir=root / "cached-temp",
                    frozen_input_manifest=manifest.path,
                    common_preprocessing_cache=cache,
                )
            )

            scientific_names = {
                *(f"{aspect}-{split}.csv" for aspect in ("bp", "cc", "mf")
                  for split in ("training", "validation", "test")),
                "train_data.pkl",
                "train_data_train.pkl",
                "train_data_valid.pkl",
                "test_data.pkl",
                "terms.pkl",
                "uniprot_to_uniref90.tsv",
                "protein_cluster_assignments.tsv",
                "cluster_split_assignments.tsv",
            }
            self.assertEqual(
                {name: sha256_file(raw.output_dir / name) for name in scientific_names},
                {name: sha256_file(cached.output_dir / name) for name in scientific_names},
            )
            input_manifest = json.loads(
                (cached.output_dir / "input_manifest.json").read_text(encoding="utf-8")
            )
            self.assertTrue(input_manifest["common_preprocessing_cache"]["used"])

    def test_tampered_cache_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config, manifest = self._fixture_manifest(root)
            specs = _input_specs(config)
            cache = build_common_preprocessing_cache(
                root / "common-cache",
                root / "cache-work",
                manifest.path,
                {name: spec.path for name, spec in specs.items() if spec.path is not None},
                source_scope=config.uniprot_source_scope,
                fixture_mode=True,
            )
            state = cache / "preprocessing_state.pkl.gz"
            state.write_bytes(state.read_bytes() + b"tamper")
            with self.assertRaisesRegex(ValueError, "file-size mismatch"):
                inspect_common_preprocessing_cache(cache, verify_file_hashes=True)

    def test_cache_marker_path_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config, manifest = self._fixture_manifest(root)
            specs = _input_specs(config)
            cache = build_common_preprocessing_cache(
                root / "common-cache",
                root / "cache-work",
                manifest.path,
                {name: spec.path for name, spec in specs.items() if spec.path is not None},
                source_scope=config.uniprot_source_scope,
                fixture_mode=True,
            )
            payload = inspect_common_preprocessing_cache(cache / CACHE_MARKER)
            self.assertEqual(payload["uniprot_source_scope"], "sprot-only")
            self.assertEqual(payload["schema_version"], CACHE_SCHEMA_VERSION)

    def test_schema_v3_ignores_legacy_nonproducer_source_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config, manifest = self._fixture_manifest(root)
            specs = _input_specs(config)
            cache = build_common_preprocessing_cache(
                root / "common-cache",
                root / "cache-work",
                manifest.path,
                {name: spec.path for name, spec in specs.items() if spec.path is not None},
                source_scope=config.uniprot_source_scope,
                fixture_mode=True,
            )
            marker = cache / CACHE_MARKER
            payload = json.loads(marker.read_text(encoding="utf-8"))
            payload["preprocessing_source_sha256"]["config.py"] = "a" * 64
            marker.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

            inspected = inspect_common_preprocessing_cache(cache)
            self.assertEqual(inspected["schema_version"], CACHE_SCHEMA_VERSION)

    def test_schema_v3_still_rejects_changed_preprocessing_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config, manifest = self._fixture_manifest(root)
            specs = _input_specs(config)
            cache = build_common_preprocessing_cache(
                root / "common-cache",
                root / "cache-work",
                manifest.path,
                {name: spec.path for name, spec in specs.items() if spec.path is not None},
                source_scope=config.uniprot_source_scope,
                fixture_mode=True,
            )
            marker = cache / CACHE_MARKER
            payload = json.loads(marker.read_text(encoding="utf-8"))
            payload["preprocessing_source_sha256"]["goa.py"] = "0" * 64
            marker.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "unsupported preprocessing code"):
                inspect_common_preprocessing_cache(cache)

    def test_schema_v2_main_module_state_is_loaded_compatibly(self):
        with tempfile.TemporaryDirectory() as tmp:
            state_path = Path(tmp) / "legacy.pkl.gz"
            script = (
                "from dataclasses import dataclass\n"
                "import gzip,pickle,sys\n"
                "@dataclass(frozen=True)\n"
                "class CommonPreprocessingState:\n"
                "  goa: object\n"
                "  catalog: object\n"
                "  decisions: list\n"
                "  requested_raw: set\n"
                "with gzip.open(sys.argv[1], 'wb') as handle:\n"
                "  pickle.dump(CommonPreprocessingState('goa','catalog',['decision'],{'P1'}), handle)\n"
            )
            subprocess.run(
                [sys.executable, "-c", script, str(state_path)], check=True
            )
            state = _load_common_preprocessing_state(state_path, 2)
            self.assertIsInstance(state, CommonPreprocessingState)
            self.assertEqual(state.decisions, ["decision"])
            self.assertEqual(state.requested_raw, {"P1"})

    def test_schema_v2_producer_fingerprint_is_exactly_pinned(self):
        self.assertEqual(
            SCHEMA_V2_PREPROCESSING_SOURCE_SHA256["common_cache.py"],
            "07eb91fe7cfa8fd3bb8c23f62d633c56ebf5e4ce1905c755dd2e6006cf146994",
        )


if __name__ == "__main__":
    unittest.main()
