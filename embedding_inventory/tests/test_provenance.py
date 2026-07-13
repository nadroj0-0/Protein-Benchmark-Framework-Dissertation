import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from dataclasses import replace
from pathlib import Path


SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from pfp_embedding_inventory.benchmark import parse_benchmark  # noqa: E402
from pfp_embedding_inventory.inventory import build_inventory  # noqa: E402
from pfp_embedding_inventory.models import (  # noqa: E402
    ArchiveSpec,
    AliasEntry,
    ArtifactScopeSpec,
    ReferenceFileSpec,
)
from pfp_embedding_inventory.provenance import (  # noqa: E402
    HashCache,
    assert_cache_catalog_unchanged,
    build_run_provenance,
    compute_cache_catalog,
    verify_artifact_scope,
)
from pfp_embedding_inventory.reports import write_reports  # noqa: E402

from helpers import DIRS, make_cache, make_config, write_nine_csvs  # noqa: E402


class ProvenanceAndArtifactTests(unittest.TestCase):
    def _verified_fixture(self, root: Path):
        benchmark_dir = root / "benchmark"
        write_nine_csvs(benchmark_dir)
        config = make_config({"text": "artifact-scoped", "structure": "artifact-scoped", "ppi": "artifact-scoped"})
        benchmark = parse_benchmark(benchmark_dir, config.target_benchmark_contract)
        artifact_root = root / "artifact"
        cache = artifact_root / "data" / "embedding_cache"
        make_cache(cache, ["P1"])
        archive = artifact_root / "published.tar.gz"
        archive.parent.mkdir(parents=True, exist_ok=True)
        archive.write_bytes(b"published archive bytes")
        reference = artifact_root / "reference.txt"
        reference.write_text("pinned source\n")
        subprocess.run(("git", "init", "-q", str(artifact_root)), check=True)
        subprocess.run(("git", "-C", str(artifact_root), "config", "user.email", "test@example.invalid"), check=True)
        subprocess.run(("git", "-C", str(artifact_root), "config", "user.name", "Inventory Test"), check=True)
        subprocess.run(("git", "-C", str(artifact_root), "add", "reference.txt"), check=True)
        subprocess.run(("git", "-C", str(artifact_root), "commit", "-q", "-m", "pin"), check=True)
        commit = subprocess.check_output(("git", "-C", str(artifact_root), "rev-parse", "HEAD"), text=True).strip()
        hashes = HashCache()
        catalog = compute_cache_catalog(cache, config, hashes)
        scope = ArtifactScopeSpec(
            mode="verified-published-cache",
            artifact_id="same-name-is-not-proof",
            metadata_url="https://example.invalid/artifact",
            expected_benchmark_fingerprint=benchmark.fingerprint,
            expected_cache_catalog_fingerprint=catalog.fingerprint,
            expected_modality_counts=catalog.modality_counts,
            expected_total_files=catalog.total_files,
            expected_total_bytes=catalog.total_bytes,
            archives=(ArchiveSpec("published.tar.gz", _sha(archive)),),
            expected_reference_commit=commit,
            reference_files=(ReferenceFileSpec("reference.txt", _sha(reference)),),
        )
        return benchmark, replace(config, artifact_scope=scope), artifact_root, cache, catalog

    def test_exact_benchmark_and_cache_pass_artifact_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            benchmark, config, artifact_root, _, catalog = self._verified_fixture(Path(tmp))
            proof = verify_artifact_scope(
                config, benchmark, benchmark, catalog, artifact_root, "paper-faithful"
            )
            self.assertTrue(proof.verified, proof.reasons)
            self.assertTrue(all(proof.checks.values()))
            aliases = {
                ("P1", "text"): [
                    AliasEntry(
                        "P1", "P1", "text", "unneeded-alias", "test-text-v1",
                        "description-sha256:" + ("a" * 64) + ";temporal-context:test-text-v1",
                    )
                ]
            }
            result = build_inventory(
                benchmark, benchmark, self._cache_from_root(artifact_root), config,
                "paper-faithful", aliases=aliases, artifact_verification=proof,
            )
            text = next(r for r in result.records if r.modality == "text")
            self.assertEqual(text.requested_action, "reuse")
            self.assertEqual(text.match_route, "exact-id")

            maximum_proof = verify_artifact_scope(
                config, benchmark, benchmark, catalog, artifact_root, "maximize-coverage"
            )
            self.assertFalse(maximum_proof.verified)
            maximum = build_inventory(
                benchmark, benchmark, self._cache_from_root(artifact_root), config,
                "maximize-coverage", artifact_verification=maximum_proof,
            )
            self.assertEqual(
                next(r for r in maximum.records if r.modality == "text").requested_action,
                "manual-review",
            )

    def test_changed_benchmark_fails_and_alias_cannot_establish_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, config, artifact_root, cache, catalog = self._verified_fixture(root)
            changed_dir = root / "changed"
            rows = {(o, s): [("P1", "ACDE", "1")] for o in ("bp", "cc", "mf") for s in ("training", "validation", "test")}
            rows[("bp", "training")] = [("P1", "ACDE", "1"), ("P2", "AAAA", "0")]
            write_nine_csvs(changed_dir, rows_by_file=rows)
            target = parse_benchmark(changed_dir, config.target_benchmark_contract)
            proof = verify_artifact_scope(
                config, target, source, catalog, artifact_root, "paper-faithful"
            )
            self.assertFalse(proof.verified)
            self.assertFalse(proof.checks["target_benchmark_fingerprint"])
            aliases = {
                ("P1", "text"): [
                    AliasEntry(
                        "P1", "P1", "text", "attempted-scope-alias", "test-text-v1",
                        "description-sha256:" + ("b" * 64) + ";temporal-context:test-text-v1",
                    )
                ]
            }
            result = build_inventory(
                target, source, cache, config, "paper-faithful",
                aliases=aliases,
                artifact_verification=proof,
            )
            text = next(r for r in result.records if r.protein_id == "P1" and r.modality == "text")
            self.assertEqual(text.requested_action, "manual-review")
            self.assertIn("config labels and aliases are not proof", text.reason)
            self.assertTrue(text.match_route.startswith("explicit-alias:"))

            canonical_proof = verify_artifact_scope(
                config, source, source, catalog, artifact_root, "paper-faithful"
            )
            with self.assertRaisesRegex(ValueError, "not bound to this target"):
                build_inventory(
                    target, source, cache, config, "paper-faithful",
                    artifact_verification=canonical_proof,
                )

    def test_archive_and_reference_corruption_invalidate_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            benchmark, config, artifact_root, _, catalog = self._verified_fixture(root)
            (artifact_root / "published.tar.gz").write_bytes(b"corrupt archive")
            archive_proof = verify_artifact_scope(
                config, benchmark, benchmark, catalog, artifact_root, "paper-faithful"
            )
            self.assertFalse(archive_proof.checks["archive:published.tar.gz"])
            (artifact_root / "reference.txt").write_text("changed source\n")
            reference_proof = verify_artifact_scope(
                config, benchmark, benchmark, catalog, artifact_root, "paper-faithful"
            )
            self.assertFalse(reference_proof.checks["reference_file:reference.txt"])

    def test_cache_catalog_changes_when_bytes_or_catalog_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, config, _, cache, catalog = self._verified_fixture(root)
            path = cache / DIRS["ppi"] / "P1.npy"
            data = bytearray(path.read_bytes())
            data[-1] ^= 1
            path.write_bytes(bytes(data))
            changed_bytes = compute_cache_catalog(cache, config)
            self.assertNotEqual(catalog.fingerprint, changed_bytes.fingerprint)
            extra = cache / DIRS["ppi"] / "EXTRA.npy"
            extra.write_bytes(path.read_bytes())
            changed_catalog = compute_cache_catalog(cache, config)
            self.assertNotEqual(changed_bytes.fingerprint, changed_catalog.fingerprint)
            self.assertEqual(changed_catalog.total_files, changed_bytes.total_files + 1)

    def test_mutation_after_artifact_proof_aborts_manifest_snapshot(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            benchmark, config, artifact_root, cache, catalog = self._verified_fixture(root)
            proof = verify_artifact_scope(
                config, benchmark, benchmark, catalog, artifact_root, "paper-faithful"
            )
            self.assertTrue(proof.verified)
            build_inventory(
                benchmark, benchmark, cache, config, "paper-faithful",
                artifact_verification=proof,
            )
            import numpy as np
            np.save(cache / DIRS["ppi"] / "P1.npy", np.ones(512, dtype=np.float32))
            with self.assertRaisesRegex(ValueError, "cache changed"):
                assert_cache_catalog_unchanged(
                    catalog, compute_cache_catalog(cache, config)
                )

    def test_run_provenance_hashes_shared_target_source_csvs_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            benchmark, config, artifact_root, cache, catalog = self._verified_fixture(root)
            proof = verify_artifact_scope(config, benchmark, benchmark, catalog, artifact_root, "paper-faithful")
            config_path = root / "config.json"
            config_path.write_text(json.dumps({"version": 1}))
            hashes = HashCache()
            first = build_run_provenance(
                command=("inventory", "--test"), repository=root,
                config_path=config_path, alias_path=None,
                target=benchmark, source=benchmark, embedding_cache=cache,
                artifact_root=artifact_root, catalog=catalog, verification=proof,
                policy="paper-faithful", report_level="compact",
                runtime_options={"test": True}, hash_cache=hashes,
            )
            for name in ("bp-training.csv", "cc-validation.csv", "mf-test.csv"):
                self.assertEqual(hashes.compute_count[(benchmark.directory / name).resolve()], 1)
            self.assertEqual(
                first["inputs"]["available_archives"]["published.tar.gz"]["sha256"],
                _sha(artifact_root / "published.tar.gz"),
            )
            old_hash = first["inputs"]["config"]["sha256"]
            config_path.write_text(json.dumps({"version": 2}))
            second = build_run_provenance(
                command=("inventory", "--test"), repository=root,
                config_path=config_path, alias_path=None,
                target=benchmark, source=benchmark, embedding_cache=cache,
                artifact_root=artifact_root, catalog=catalog, verification=proof,
                policy="paper-faithful", report_level="compact",
                runtime_options={"test": True}, hash_cache=HashCache(),
            )
            self.assertNotEqual(old_hash, second["inputs"]["config"]["sha256"])

            unverified_config = make_config()
            unverified_catalog = compute_cache_catalog(cache, unverified_config)
            unverified_proof = verify_artifact_scope(
                unverified_config, benchmark, benchmark, unverified_catalog,
                artifact_root, "paper-faithful",
            )
            unverified = build_run_provenance(
                command=("inventory", "--test"), repository=root,
                config_path=config_path, alias_path=None,
                target=benchmark, source=benchmark, embedding_cache=cache,
                artifact_root=artifact_root, catalog=unverified_catalog,
                verification=unverified_proof, policy="paper-faithful",
                report_level="compact", runtime_options={}, hash_cache=HashCache(),
            )
            self.assertEqual(unverified["inputs"]["available_archives"], {})

    def test_compact_has_full_inventory_fields_and_full_adds_sequences(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            benchmark, config, artifact_root, cache, catalog = self._verified_fixture(root)
            proof = verify_artifact_scope(config, benchmark, benchmark, catalog, artifact_root, "paper-faithful")
            result = build_inventory(benchmark, benchmark, cache, config, "paper-faithful", artifact_verification=proof)
            write_reports(result, root / "compact", cache, report_level="compact")
            self.assertTrue((root / "compact" / "embedding_inventory.tsv.gz").exists())
            self.assertFalse((root / "compact" / "benchmark_proteins_full.tsv.gz").exists())
            self.assertTrue((root / "compact" / "RUN_COMPLETE.json").exists())
            self.assertTrue((root / "compact" / "output_manifest.json").exists())
            write_reports(result, root / "full", cache, report_level="full")
            detailed = root / "full" / "embedding_inventory.tsv.gz"
            self.assertTrue(detailed.is_file())
            import gzip
            with gzip.open(detailed, "rt") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
                self.assertEqual(len(rows), 4)
                for field in (
                    "observed_shape", "expected_shape", "dtype", "finite",
                    "scientifically_eligible",
                ):
                    self.assertIn(field, reader.fieldnames)
            self.assertTrue((root / "full" / "benchmark_proteins_full.tsv.gz").exists())

    def test_output_cannot_mutate_artifact_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            benchmark, config, artifact_root, cache, catalog = self._verified_fixture(root)
            proof = verify_artifact_scope(
                config, benchmark, benchmark, catalog, artifact_root, "paper-faithful"
            )
            result = build_inventory(
                benchmark, benchmark, cache, config, "paper-faithful",
                artifact_verification=proof,
            )
            with self.assertRaisesRegex(ValueError, "cannot be inside"):
                write_reports(
                    result, artifact_root / "reports", cache,
                    protected_roots=(artifact_root,),
                )

    def test_failed_report_write_is_atomic_and_leaves_no_destination(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            benchmark, config, artifact_root, cache, catalog = self._verified_fixture(root)
            proof = verify_artifact_scope(
                config, benchmark, benchmark, catalog, artifact_root, "paper-faithful"
            )
            result = build_inventory(
                benchmark, benchmark, cache, config, "paper-faithful",
                artifact_verification=proof,
            )
            output = root / "failed"
            with patch(
                "pfp_embedding_inventory.reports._write_modality_lists",
                side_effect=KeyboardInterrupt("injected report interruption"),
            ):
                with self.assertRaisesRegex(KeyboardInterrupt, "injected report interruption"):
                    write_reports(result, output, cache)
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".failed.staging-*")), [])

    @staticmethod
    def _cache_from_root(artifact_root: Path) -> Path:
        return artifact_root / "data" / "embedding_cache"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__ == "__main__":
    unittest.main()
