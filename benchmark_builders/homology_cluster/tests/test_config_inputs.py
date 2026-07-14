from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from homology_cluster_benchmark.inputs import resolve_input
from homology_cluster_benchmark.frozen_inputs import (
    bind_frozen_inputs,
    load_frozen_input_manifest,
    write_synthetic_fixture_manifest,
)
from homology_cluster_benchmark.models import InputSpec

from tests.helpers import FIXTURES, fixture_config


class ConfigAndInputTests(unittest.TestCase):
    def test_binding_80_20_fraction_and_frozen_releases_are_locked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaisesRegex(ValueError, "0.80/0.20"):
                fixture_config(
                    root / "out", root / "temp", development_fraction=0.5
                ).validate()
            with self.assertRaisesRegex(ValueError, "frozen"):
                fixture_config(
                    root / "out", root / "temp", release_goa="999"
                ).validate()
            for field, value, message in (
                ("training_fraction_within_development", 0.8, "0.90/0.10"),
                ("sensitivity", 6.0, "sensitivity"),
                ("evalue", 1e-3, "E-value"),
            ):
                with self.subTest(field=field):
                    with self.assertRaisesRegex(ValueError, message):
                        replace(
                            fixture_config(root / field, root / "temp"),
                            **{field: value},
                        ).validate()

    def test_precomputed_assignments_and_low_min_count_require_fixture_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = fixture_config(root / "out", root / "temp")
            with self.assertRaisesRegex(ValueError, "fixture-only"):
                replace(config, fixture_mode=False, min_count=50).validate()
            with self.assertRaisesRegex(ValueError, "at least 50"):
                replace(
                    config, fixture_mode=False, cluster_assignments=None, min_count=1
                ).validate(require_pinned_inputs=False)

            # A future representative-subset smoke may execute MMseqs2 directly while
            # remaining explicitly non-production.
            replace(config, cluster_assignments=None).validate()

    def test_production_requires_manifest_exact_version_and_expected_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = replace(
                fixture_config(root / "out", root / "temp"),
                fixture_mode=False,
                cluster_assignments=None,
                min_count=50,
            )
            with self.assertRaisesRegex(ValueError, "frozen-input-manifest"):
                config.validate()
            with self.assertRaisesRegex(ValueError, "expected SHA-256"):
                replace(
                    config,
                    frozen_input_manifest=root / "declared.json",
                    expected_mmseqs_version="15-6f452",
                ).validate()

    def test_frozen_manifest_schema_binding_and_fixture_boundary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = fixture_config(root / "out", root / "temp")
            specs = {
                name: getattr(config, name)
                for name in (
                    "uniref90_fasta", "idmapping", "uniprot_sequences", "goa", "go_obo"
                )
            }
            resolved = {
                name: resolve_input(spec, root / "downloads", allow_downloads=False)
                for name, spec in specs.items()
            }
            manifest_path = root / "fixture-manifest.json"
            manifest = write_synthetic_fixture_manifest(manifest_path, specs, resolved)
            self.assertFalse(manifest.authoritative_origin_recorded)
            self.assertTrue(all(bind_frozen_inputs(manifest, specs, resolved).values()) is False)

            payload = json.loads(manifest_path.read_text())
            payload["inputs"][1]["name"] = payload["inputs"][0]["name"]
            duplicate = root / "duplicate.json"
            duplicate.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                load_frozen_input_manifest(duplicate, fixture_mode=True)

            payload = json.loads(manifest_path.read_text())
            payload["inputs"][0]["size_bytes"] += 1
            mismatch = root / "mismatch.json"
            mismatch.write_text(json.dumps(payload))
            mismatched = load_frozen_input_manifest(mismatch, fixture_mode=True)
            with self.assertRaisesRegex(ValueError, "byte-size mismatch"):
                bind_frozen_inputs(mismatched, specs, resolved)

            payload = json.loads(manifest_path.read_text())
            payload["review"] = {
                "status": "self-hashed",
                "authoritative_origin": False,
                "reviewed_by": "local",
                "evidence": "computed from the same bytes",
            }
            self_hashed = root / "self-hashed.json"
            self_hashed.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "self-hashes alone"):
                load_frozen_input_manifest(self_hashed, fixture_mode=False)

    def test_manifest_template_placeholders_cannot_validate(self):
        template = ROOT / "frozen_input_manifest.template.json"
        with self.assertRaisesRegex(ValueError, "placeholder|positive"):
            load_frozen_input_manifest(template, fixture_mode=False)

    def test_local_input_precedes_url_and_hash_is_enforced_without_network(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "input.txt"
            path.write_text("frozen\n")
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            resolved = resolve_input(
                InputSpec(
                    "example", path=path, url="https://example.invalid/current",
                    expected_sha256=digest, release="frozen",
                ),
                root / "downloads",
                allow_downloads=False,
            )
            self.assertEqual(resolved.acquisition, "local")
            self.assertEqual(resolved.source_url, "https://example.invalid/current")
            self.assertEqual(resolved.sha256, digest)
            self.assertFalse((root / "downloads").exists())
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                resolve_input(
                    InputSpec(
                        "example", path=path, expected_sha256="0" * 64, release="frozen"
                    ),
                    root / "downloads",
                    allow_downloads=False,
                )

    def test_fixture_paths_exist(self):
        self.assertTrue((FIXTURES / "uniref90.fasta").is_file())


if __name__ == "__main__":
    unittest.main()
