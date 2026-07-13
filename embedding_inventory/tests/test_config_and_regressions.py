import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SRC = PACKAGE_ROOT / "src"
REPO_ROOT = PACKAGE_ROOT.parent
sys.path.insert(0, str(SRC))

from pfp_embedding_inventory.config import ConfigError, load_config  # noqa: E402
from helpers import unique_rows_by_file, write_nine_csvs  # noqa: E402


class ConfigTests(unittest.TestCase):
    def _base_config(self):
        modalities = {}
        for name, directory, dimension in (
            ("prott5", "prott5", 1024),
            ("text", "exp_text_embeddings_temporal", 768),
            ("structure", "IF1", 512),
            ("ppi", "ppi", 512),
        ):
            modalities[name] = {
                "directory": directory,
                "expected_dim": dimension,
                "sequence_dependent": name in {"prott5", "structure"},
                "allow_sequence_hash_reuse": name == "prott5",
                "missing_action": "generate",
                "invalid_action": "generate",
                "provenance": {
                    "compatibility": "compatible",
                    "label": name,
                    "source_identity": name + "-v1",
                    "target_identity": name + "-v1",
                    "evidence": "test evidence",
                    "requires_mapping_evidence": name != "prott5",
                },
            }
        return {
            "schema_version": 2,
            "name": "test",
            "target_benchmark_contract": {
                "id_overlap": "allow", "sequence_overlap": "allow",
                "protein_id_pattern": "^[^\\s/\\\\]+$",
                "sequence_pattern": "^[A-Za-z*.-]+$",
            },
            "source_benchmark_contract": {
                "id_overlap": "allow", "sequence_overlap": "allow",
                "protein_id_pattern": "^[^\\s/\\\\]+$",
                "sequence_pattern": "^[A-Za-z*.-]+$",
            },
            "artifact_scope": {"mode": "none"},
            "modalities": modalities,
        }

    def test_schema_one_is_rejected_instead_of_silently_reinterpreted(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            payload = self._base_config()
            payload["schema_version"] = 1
            payload["benchmark_contract"] = payload.pop("target_benchmark_contract")
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ConfigError, "schema_version 1 is not accepted"):
                load_config(path)

    def test_non_object_config_root_is_rejected_cleanly(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            path.write_text("[]")
            with self.assertRaisesRegex(ConfigError, "root must be a JSON object"):
                load_config(path)

    def test_source_and_target_contracts_are_independent(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            payload = self._base_config()
            payload["target_benchmark_contract"]["id_overlap"] = "global-evaluation-disjoint"
            payload["source_benchmark_contract"]["id_overlap"] = "allow"
            path.write_text(json.dumps(payload))
            config = load_config(path)
            self.assertEqual(config.target_benchmark_contract.id_overlap, "global-evaluation-disjoint")
            self.assertEqual(config.source_benchmark_contract.id_overlap, "allow")

    def test_contract_typo_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            payload = self._base_config()
            payload["target_benchmark_contract"]["id_overalp"] = payload[
                "target_benchmark_contract"
            ].pop("id_overlap")
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ConfigError, "missing required keys: id_overlap"):
                load_config(path)

    def test_text_source_must_be_singular(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            payload = self._base_config()
            payload["modalities"]["text"]["directory"] = [
                "exp_text_embeddings",
                "exp_text_embeddings_temporal",
            ]
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ConfigError, "singular"):
                load_config(path)

    def test_claimed_compatible_identity_must_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            payload = self._base_config()
            payload["modalities"]["ppi"]["provenance"]["target_identity"] = "STRING-v13"
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ConfigError, "source and target identities differ"):
                load_config(path)

    def test_modality_directory_cannot_escape_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            payload = self._base_config()
            payload["modalities"]["structure"]["directory"] = "../outside"
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ConfigError, "cache-relative"):
                load_config(path)

    def test_required_modality_invariants_cannot_be_disabled(self):
        mutations = [
            ("prott5", "expected_dim", 512, "PFP dimension 1024"),
            ("structure", "sequence_dependent", False, "sequence_dependent must be true"),
            ("text", "allow_sequence_hash_reuse", True, "allow_sequence_hash_reuse must be false"),
            ("prott5", "sequence_dependent", "yes", "must be a JSON boolean"),
        ]
        for modality, key, value, message in mutations:
            with self.subTest(modality=modality, key=key), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "config.json"
                payload = self._base_config()
                payload["modalities"][modality][key] = value
                path.write_text(json.dumps(payload))
                with self.assertRaisesRegex(ConfigError, message):
                    load_config(path)


class LegacyRegressionTests(unittest.TestCase):
    def test_inventory_cli_writes_compact_provenance_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target"
            source = root / "source"
            write_nine_csvs(target, rows_by_file=unique_rows_by_file())
            write_nine_csvs(source, rows_by_file=unique_rows_by_file())
            cache = root / "cache"
            cache.mkdir()
            artifact_root = root / "artifact"
            artifact_root.mkdir()
            output = root / "output"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "verification" / "inventory_embeddings.py"),
                    "--benchmark-dir", str(target),
                    "--source-benchmark-dir", str(source),
                    "--embedding-cache", str(cache),
                    "--artifact-root", str(artifact_root),
                    "--config", str(REPO_ROOT / "configs" / "embedding_inventory.contemporary.json"),
                    "--policy", "maximize-coverage",
                    "--report-level", "compact",
                    "--output-dir", str(output),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
            provenance = json.loads((output / "run_provenance.json").read_text())
            self.assertEqual(provenance["run"]["report_level"], "compact")
            self.assertEqual(len(provenance["inputs"]["target_csvs"]), 9)
            self.assertEqual(len(provenance["inputs"]["source_csvs"]), 9)
            self.assertIn("--benchmark-dir", provenance["command"])
            self.assertIs(provenance["software"]["dirty_worktree"], True)
            self.assertEqual(provenance["software"]["git_status_error"], "")
            self.assertTrue((output / "embedding_inventory.tsv.gz").exists())
            self.assertFalse((output / "benchmark_proteins_full.tsv.gz").exists())
            self.assertTrue(json.loads((output / "RUN_COMPLETE.json").read_text())["complete"])

    def test_generate_embeddings_fasta_still_works(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp)
            (data / "BPO_train_sequences.json").write_text(
                json.dumps({"P1": "ACDE", "P2": "FGHI"})
            )
            output = data / "proteins.fasta"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "embeddings" / "generate_embeddings_fasta.py"),
                    "--data-dir",
                    str(data),
                    "--out",
                    str(output),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
            self.assertEqual(output.read_text(), ">P1\nACDE\n>P2\nFGHI\n")

    def test_verify_embeddings_still_works(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp)
            np.save(data / "BPO_train_names.npy", np.array(["P1"], dtype=object))
            emb_dir = data / "embedding_cache" / "prott5"
            emb_dir.mkdir(parents=True)
            np.save(emb_dir / "P1.npy", np.zeros(4, dtype=np.float32))
            config = data / "legacy_config.json"
            config.write_text(
                json.dumps(
                    {
                        "aspects": ["BPO"],
                        "splits": ["train"],
                        "cache_dir": "embedding_cache",
                        "sample_size": 10,
                        "catastrophic_factor": 0.5,
                        "modalities": {
                            "prott5": {
                                "dirs": ["prott5"],
                                "dim": 4,
                                "min_coverage": 1.0,
                            }
                        },
                    }
                )
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "verification" / "verify_embeddings.py"),
                    "--data-dir",
                    str(data),
                    "--config",
                    str(config),
                    "--strict",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
            self.assertIn("RESULT: PASS", completed.stdout)


if __name__ == "__main__":
    unittest.main()
