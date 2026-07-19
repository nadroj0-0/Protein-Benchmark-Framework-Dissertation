from __future__ import annotations

import hashlib
import io
import csv
import gzip
import json
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np


FRAMEWORK_ROOT = Path(__file__).resolve().parents[3]
EMBEDDING_SCRIPTS = FRAMEWORK_ROOT / "scripts" / "embeddings"
sys.path.insert(0, str(EMBEDDING_SCRIPTS))

from finalize_embedding_state import (  # noqa: E402
    build_final_evidence,
    canonical_sha256,
    retire_source_embeddings,
    sha256_file,
    source_snapshot,
)
from manage_embedding_archive import create_archive, extract_archive  # noqa: E402
from build_embedding_baseline_archive import build_baseline  # noqa: E402


CONFIG = FRAMEWORK_ROOT / "configs" / "pfp_benchmark_run.temporal.json"
HPC_RUNNER = FRAMEWORK_ROOT / "hpc_jobs" / "active" / "hpc_pfp_benchmark.sh"
STATE_MANAGER = EMBEDDING_SCRIPTS / "manage_resumable_embedding_state.py"
FINALIZER = EMBEDDING_SCRIPTS / "finalize_embedding_state.py"

FAKE_PREPARER = r'''#!/usr/bin/env python3
import argparse, csv, json
from pathlib import Path
import numpy as np
import scipy.sparse as ssp

p = argparse.ArgumentParser()
p.add_argument("--cafa3-dir", required=True)
p.add_argument("--output-dir", required=True)
a = p.parse_args()
source, output = Path(a.cafa3_dir), Path(a.output_dir)
output.mkdir(parents=True, exist_ok=True)
for aspect, prefix in (("BPO", "bp"), ("CCO", "cc"), ("MFO", "mf")):
    terms = None
    for csv_split, pfp_split in (("training", "train"), ("validation", "valid"), ("test", "test")):
        with (source / f"{prefix}-{csv_split}.csv").open(newline="") as handle:
            rows = list(csv.reader(handle))
        header, values = rows[0], rows[1:]
        if terms is None:
            terms = header[2:]
        lookup = [header.index(term) for term in terms]
        names = np.array([row[0] for row in values], dtype=object)
        labels = np.array([[float(row[index]) for index in lookup] for row in values], dtype=np.float32)
        np.save(output / f"{aspect}_{pfp_split}_names.npy", names)
        ssp.save_npz(output / f"{aspect}_{pfp_split}_labels.npz", ssp.csr_matrix(labels))
        (output / f"{aspect}_{pfp_split}_sequences.json").write_text(
            json.dumps({row[0]: row[1] for row in values})
        )
    (output / f"{aspect}_go_terms.json").write_text(json.dumps(terms))
    (output / f"{aspect}_info.json").write_text(json.dumps({"aspect": aspect}))
'''


class EmbeddingFinalizationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def make_cache(self) -> Path:
        cache = self.root / "cache"
        for index, directory in enumerate(
            ("prott5", "exp_text_embeddings_temporal", "IF1", "ppi"), start=1
        ):
            path = cache / directory
            path.mkdir(parents=True)
            np.save(path / "P1.npy", np.asarray([index, index + 1], dtype=np.float32))
        return cache

    def test_archive_roundtrip_preserves_every_array_byte(self) -> None:
        cache = self.make_cache()
        archive = self.root / "cache.tar.gz"
        created = create_archive(cache, archive, CONFIG)
        extracted = self.root / "extracted"
        observed = extract_archive(archive, extracted, CONFIG)
        self.assertEqual(created["archive_sha256"], observed["archive_sha256"])
        self.assertEqual(created["member_count"], 4)
        self.assertEqual(
            created["member_content_sha256"], observed["member_content_sha256"]
        )
        for source in cache.glob("*/*.npy"):
            destination = extracted / source.parent.name / source.name
            self.assertEqual(source.read_bytes(), destination.read_bytes())

    def test_generated_cache_becomes_one_archive_backed_baseline(self) -> None:
        data = self.root / "data"
        data.mkdir()
        (data / "BPO_train_sequences.json").write_text(
            json.dumps({"P1": "ACDE", "P2": "FGHI"}), encoding="utf-8"
        )
        policy = self.root / "policy.json"
        policy.write_text(json.dumps({
            "schema_version": 1,
            "modalities": {
                "sequence": {"cache_directory": "prott5", "dimension": 2,
                             "min_accepted_count": 2},
                "text": {"cache_directory": "exp_text_embeddings_temporal", "dimension": 2,
                         "min_accepted_count": 1},
                "structure": {"cache_directory": "IF1", "dimension": 2,
                              "min_accepted_count": 1},
                "ppi": {"cache_directory": "ppi", "dimension": 2,
                        "min_accepted_count": 1},
            },
        }))
        cache = self.root / "generated"
        for directory in ("prott5", "exp_text_embeddings_temporal", "IF1", "ppi"):
            (cache / directory).mkdir(parents=True)
            np.save(cache / directory / "P1.npy", np.ones(2, dtype=np.float32))
        np.save(cache / "prott5/P2.npy", np.full(2, 2, dtype=np.float32))
        # Real text generation leaves auxiliary directories; they must not be
        # persisted into the scientific embedding baseline.
        (cache / "uniprot_text").mkdir()
        (cache / "uniprot_text/descriptions.tsv").write_text("raw\n")

        archive = self.root / "baseline/cache.tar.gz"
        assembly = self.root / "baseline/embedding_assembly.tsv.gz"
        summary = build_baseline(cache, data, policy, archive, assembly)

        self.assertEqual(summary["archive_member_count"], 5)
        self.assertEqual(summary["available_by_modality"]["sequence"], 2)
        self.assertEqual(summary["missing_pairs"], 3)
        with gzip.open(assembly, "rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 8)
        self.assertEqual(
            sum(row["status"] == "available" for row in rows), 5
        )
        with tarfile.open(archive, "r:gz") as handle:
            names = {member.name for member in handle if member.isfile()}
        self.assertNotIn(
            "data/embedding_cache/uniprot_text/descriptions.tsv", names
        )

    def test_baseline_rejects_unreduced_text_hidden_state(self) -> None:
        data = self.root / "data"
        data.mkdir()
        (data / "BPO_train_sequences.json").write_text(
            json.dumps({"P1": "ACDE"}), encoding="utf-8"
        )
        policy = self.root / "policy.json"
        policy.write_text(json.dumps({
            "schema_version": 1,
            "modalities": {
                "sequence": {"cache_directory": "prott5", "dimension": 1024,
                             "min_accepted_count": 1},
                "text": {"cache_directory": "exp_text_embeddings_temporal",
                         "dimension": 768, "min_accepted_count": 1},
                "structure": {"cache_directory": "IF1", "dimension": 512,
                              "min_accepted_count": 1},
                "ppi": {"cache_directory": "ppi", "dimension": 512,
                        "min_accepted_count": 1},
            },
        }))
        cache = self.root / "generated"
        specifications = {
            "prott5": np.ones(1024, dtype=np.float32),
            "exp_text_embeddings_temporal": np.ones(
                (1, 512, 768), dtype=np.float32
            ),
            "IF1": np.ones(512, dtype=np.float32),
            "ppi": np.ones(512, dtype=np.float32),
        }
        for directory, array in specifications.items():
            (cache / directory).mkdir(parents=True)
            np.save(cache / directory / "P1.npy", array)

        archive = self.root / "baseline/cache.tar.gz"
        assembly = self.root / "baseline/embedding_assembly.tsv.gz"
        with self.assertRaisesRegex(ValueError, r"expected=\(768,\)"):
            build_baseline(cache, data, policy, archive, assembly)
        self.assertFalse(archive.exists())
        self.assertFalse(assembly.exists())

    def test_extraction_rejects_path_traversal_and_removes_partial_output(self) -> None:
        archive = self.root / "unsafe.tar.gz"
        payload = b"unsafe"
        with tarfile.open(archive, "w:gz") as handle:
            member = tarfile.TarInfo("data/embedding_cache/prott5/../../outside.npy")
            member.size = len(payload)
            handle.addfile(member, io.BytesIO(payload))
        output = self.root / "unsafe-output"
        with self.assertRaisesRegex(ValueError, "outside data/embedding_cache"):
            extract_archive(archive, output, CONFIG)
        self.assertFalse(output.exists())

    def make_source_evidence(self) -> Path:
        state = self.root / "state"
        state.mkdir()
        baseline = self.root / "baseline.tar.gz"
        baseline.write_bytes(b"baseline")
        contract = {
            "schema_version": 1,
            "benchmark_id": "fixture",
            "benchmark_csvs": [],
            "targets": {"count": 1, "manifest_sha256": "a" * 64},
            "pfp_commit": "1e04fd6d6d3c40458fd41ec1a881ed6e24de768e",
            "framework_commit": "b" * 40,
            "policy": {"modalities": {}},
            "policy_sha256": "c" * 64,
            "environment": None,
            "source_files": [],
            "runtime": {},
            "baseline": {
                "archive": {"path": str(baseline), "sha256": sha256_file(baseline)},
                "assembly_report": {"path": str(self.root / "assembly.tsv")},
            },
        }
        contract["contract_sha256"] = canonical_sha256(contract)
        (state / "contract.json").write_text(json.dumps(contract) + "\n")
        coverage = {
            "contract_sha256": contract["contract_sha256"],
            "target_count": 1,
            "embedding_gate_passed": False,
            "coverage": {"sequence": {"accepted": 1}},
        }
        (state / "coverage.json").write_text(json.dumps(coverage) + "\n")
        (state / "targets.tsv").write_text("protein_id\tsequence_sha256\nP1\t" + "d" * 64 + "\n")
        (state / "pair_status.tsv").write_text(
            "protein_id\tmodality\tstate\tsequence_sha256\tembedding_sha256\t"
            "attempts\tlatest_reason\tlatest_detail\n"
            "P1\tsequence\taccepted\t"
            + "d" * 64
            + "\t"
            + "e" * 64
            + "\t0\t\t\n"
        )
        (state / "baseline_accepted.tsv").write_text(
            "protein_id\tmodality\tarchive_member\tembedding_sha256\n"
        )
        (state / "EVIDENCE_HASHES_COMPLETE.json").write_text("{}\n")
        delta = state / "cache" / "prott5"
        delta.mkdir(parents=True)
        np.save(delta / "P1.npy", np.ones(2, dtype=np.float32))
        return state

    def test_final_evidence_is_self_contained_and_bound_to_archive(self) -> None:
        state = self.make_source_evidence()
        archive_report = {
            "archive_sha256": "f" * 64,
            "archive_size_bytes": 123,
            "member_count": 1,
            "member_content_sha256": "0" * 64,
        }
        final_root = self.root / "final"
        evidence = self.root / "evidence"
        summary = build_final_evidence(
            state, evidence, final_root, "cache.tar.gz", archive_report
        )
        contract = json.loads((evidence / "contract.json").read_text())
        recorded = contract.pop("contract_sha256")
        self.assertEqual(recorded, canonical_sha256(contract))
        self.assertNotIn("baseline", contract)
        self.assertEqual(
            contract["finalized_embedding_archive"]["path"],
            str((final_root / "cache.tar.gz").resolve()),
        )
        coverage = json.loads((evidence / "coverage.json").read_text())
        self.assertEqual(coverage["contract_sha256"], recorded)
        self.assertEqual(summary["accepted_counts"], {"sequence": 1})

    def test_source_retirement_requires_published_validated_archive(self) -> None:
        state = self.make_source_evidence()
        snapshot = source_snapshot(state)
        baseline = Path(json.loads((state / "contract.json").read_text())["baseline"]["archive"]["path"])
        final = self.root / "final"
        final.mkdir()
        archive = final / "cache.tar.gz"
        archive.write_bytes(b"final")
        with self.assertRaisesRegex(ValueError, "validation marker"):
            retire_source_embeddings(state, snapshot, final)
        self.assertTrue(baseline.is_file())
        self.assertTrue((state / "cache").is_dir())

        (final / "CACHE_ARCHIVE_VALIDATED.json").write_text(
            json.dumps(
                {
                    "archive_name": archive.name,
                    "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
                }
            )
            + "\n"
        )
        result = retire_source_embeddings(state, snapshot, final)
        self.assertFalse(baseline.exists())
        self.assertFalse((state / "cache").exists())
        self.assertTrue(archive.is_file())
        self.assertTrue(result["source_evidence_retained"])

    def test_generic_hpc_runner_accepts_archive_without_removing_directory_mode(self) -> None:
        content = HPC_RUNNER.read_text(encoding="utf-8")
        self.assertIn("--embedding-cache-archive", content)
        self.assertIn("manage_embedding_archive.py", content)
        self.assertIn("--embedding-cache-root", content)

    def test_generic_hpc_runner_activates_mmfp_before_archive_extraction(self) -> None:
        content = HPC_RUNNER.read_text(encoding="utf-8")
        activation = content.index("activate_or_create_mmfp_env")
        extraction = content.index(
            '"$PYTHON_BIN" "$FRAMEWORK_DIR/scripts/embeddings/manage_embedding_archive.py"'
        )
        self.assertLess(activation, extraction)

    def test_generic_hpc_runner_binds_results_and_exports_verified_git_state(self) -> None:
        content = HPC_RUNNER.read_text(encoding="utf-8")
        activation = content.index("activate_or_create_mmfp_env")
        self.assertLess(
            content.index('add_mmfp_singularity_bind "$RESULTS_ROOT"'), activation
        )
        self.assertIn(
            'export PFP_HOST_GIT_VERIFIED_COMMIT="$PFP_COMMIT"', content
        )
        self.assertIn(
            'export FRAMEWORK_HOST_GIT_VERIFIED_COMMIT="$FRAMEWORK_COMMIT"',
            content,
        )

    def test_incompatible_ontology_fails_before_embedding_state_is_touched(self) -> None:
        benchmark = self.root / "benchmark"
        benchmark.mkdir()
        proteins = {
            "training": ("TRAIN", "AAAA"),
            "validation": ("VALID", "CCCC"),
            "test": ("TEST", "GGGG"),
        }
        terms = {
            "bp": ("GO:0000001", "biological_process"),
            "cc": ("GO:0000002", "cellular_component"),
            "mf": ("GO:0000003", "molecular_function"),
        }
        for prefix, (term, _) in terms.items():
            for split, (protein_id, sequence) in proteins.items():
                with (benchmark / f"{prefix}-{split}.csv").open(
                    "w", encoding="utf-8", newline=""
                ) as handle:
                    writer = csv.writer(handle)
                    writer.writerow(("proteins", "sequences", term))
                    writer.writerow((protein_id, sequence, "1"))

        incompatible_obo = self.root / "incompatible.obo"
        incompatible_obo.write_text(
            "format-version: 1.2\n\n"
            "[Term]\nid: GO:0000001\nname: obsolete fixture process\n"
            "namespace: biological_process\nis_obsolete: true\n\n"
            "[Term]\nid: GO:0000002\nnamespace: cellular_component\n\n"
            "[Term]\nid: GO:0000003\nnamespace: molecular_function\n",
            encoding="utf-8",
        )
        pfp = self.root / "PFP"
        (pfp / "scripts").mkdir(parents=True)
        (pfp / "scripts/prepare_cafa3_data.py").write_text(
            FAKE_PREPARER, encoding="utf-8"
        )
        state = self.root / "retry_state"
        state.mkdir()
        sentinel = state / "untouched.txt"
        sentinel.write_text("original\n", encoding="utf-8")

        work = self.root / "work"
        report = self.root / "report"
        result = subprocess.run(
            [
                sys.executable,
                str(FINALIZER),
                "--state-root",
                str(state),
                "--benchmark-dir",
                str(benchmark),
                "--obo-file",
                str(incompatible_obo),
                "--pfp-root",
                str(pfp),
                "--config",
                str(CONFIG),
                "--work-dir",
                str(work),
                "--final-root",
                str(self.root / "final"),
                "--report-dir",
                str(report),
                "--confirm-retries-finished",
                "--retire-source-embeddings",
            ],
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("BPO GO/OBO contract failed", result.stdout + result.stderr)
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "original\n")
        self.assertFalse((state / ".state.lock").exists())
        self.assertFalse((work / "hydrated_cache").exists())
        self.assertFalse((report / "reports/evidence_hash_upgrade.json").exists())
        self.assertTrue((report / "FINALIZATION_FAILED.json").is_file())

    def test_tiny_end_to_end_finalization_publishes_before_retirement(self) -> None:
        benchmark = self.root / "benchmark"
        benchmark.mkdir()
        proteins = {
            "training": ("TRAIN", "AAAA"),
            "validation": ("VALID", "CCCC"),
            "test": ("TEST", "GGGG"),
        }
        terms = {
            "bp": ("GO:0000001", "biological_process"),
            "cc": ("GO:0000002", "cellular_component"),
            "mf": ("GO:0000003", "molecular_function"),
        }
        for prefix, (term, _) in terms.items():
            for split, (protein_id, sequence) in proteins.items():
                with (benchmark / f"{prefix}-{split}.csv").open(
                    "w", encoding="utf-8", newline=""
                ) as handle:
                    writer = csv.writer(handle)
                    writer.writerow(("proteins", "sequences", term))
                    writer.writerow((protein_id, sequence, "1"))

        obo = self.root / "go.obo"
        obo.write_text(
            "format-version: 1.2\n\n"
            + "\n\n".join(
                f"[Term]\nid: {term}\nnamespace: {namespace}"
                for term, namespace in terms.values()
            )
            + "\n",
            encoding="utf-8",
        )
        pfp = self.root / "PFP"
        (pfp / "scripts").mkdir(parents=True)
        (pfp / "scripts/prepare_cafa3_data.py").write_text(
            FAKE_PREPARER, encoding="utf-8"
        )

        targets = self.root / "targets.tsv"
        with targets.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("protein_id", "sequence", "sequence_sha256"))
            for protein_id, sequence in (value for value in proteins.values()):
                writer.writerow(
                    (protein_id, sequence, hashlib.sha256(sequence.encode()).hexdigest())
                )
        policy = self.root / "policy.json"
        policy.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "modalities": {
                        "sequence": {
                            "cache_directory": "prott5",
                            "dimension": 1024,
                            "min_accepted_count": 3,
                        },
                        "text": {
                            "cache_directory": "exp_text_embeddings_temporal",
                            "dimension": 768,
                            "min_accepted_count": 3,
                        },
                        "structure": {
                            "cache_directory": "IF1",
                            "dimension": 512,
                            "min_accepted_count": 3,
                        },
                        "ppi": {
                            "cache_directory": "ppi",
                            "dimension": 512,
                            "min_accepted_count": 3,
                        },
                    },
                }
            ),
            encoding="utf-8",
        )
        modality_specs = {
            "prott5": ("prott5", 1024),
            "text": ("exp_text_embeddings_temporal", 768),
            "structure": ("IF1", 512),
            "ppi": ("ppi", 512),
        }
        package = self.root / "package/data/embedding_cache"
        for protein_id, _ in proteins.values():
            for report_modality, (directory, dimension) in modality_specs.items():
                if protein_id == "TEST" and report_modality == "ppi":
                    continue
                destination = package / directory
                destination.mkdir(parents=True, exist_ok=True)
                np.save(destination / f"{protein_id}.npy", np.ones(dimension, dtype=np.float32))
        assembly = self.root / "assembly.tsv.gz"
        with gzip.open(assembly, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("protein_id", "modality", "status", "dimension"),
                delimiter="\t",
            )
            writer.writeheader()
            for protein_id, _ in proteins.values():
                for report_modality, (_, dimension) in modality_specs.items():
                    writer.writerow(
                        {
                            "protein_id": protein_id,
                            "modality": report_modality,
                            "status": (
                                "missing"
                                if protein_id == "TEST" and report_modality == "ppi"
                                else "available"
                            ),
                            "dimension": dimension,
                        }
                    )
        baseline = self.root / "baseline.tar.gz"
        with tarfile.open(baseline, "w:gz") as handle:
            handle.add(self.root / "package/data", arcname="data")

        state = self.root / "retry_state"
        subprocess.run(
            [
                sys.executable,
                str(STATE_MANAGER),
                "initialize",
                "--state-root",
                str(state),
                "--benchmark-id",
                "fixture",
                "--benchmark-dir",
                str(benchmark),
                "--target-table",
                str(targets),
                "--policy",
                str(policy),
                "--pfp-commit",
                "1e04fd6d6d3c40458fd41ec1a881ed6e24de768e",
                "--framework-commit",
                "a" * 40,
                "--baseline-archive",
                str(baseline),
                "--baseline-assembly-report",
                str(assembly),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        generated = self.root / "generated/ppi"
        generated.mkdir(parents=True)
        np.save(generated / "TEST.npy", np.full(512, 2, dtype=np.float32))
        requested = self.root / "requested.tsv"
        requested.write_text("protein_id\tmodality\nTEST\tppi\n", encoding="utf-8")
        subprocess.run(
            [
                sys.executable,
                str(STATE_MANAGER),
                "merge",
                "--state-root",
                str(state),
                "--generated-cache-root",
                str(self.root / "generated"),
                "--attempt-id",
                "fixture-retry",
                "--requested-pairs",
                str(requested),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        final = self.root / "final"
        result = subprocess.run(
            [
                sys.executable,
                str(FINALIZER),
                "--state-root",
                str(state),
                "--benchmark-dir",
                str(benchmark),
                "--obo-file",
                str(obo),
                "--pfp-root",
                str(pfp),
                "--config",
                str(CONFIG),
                "--work-dir",
                str(self.root / "work"),
                "--final-root",
                str(final),
                "--report-dir",
                str(self.root / "report"),
                "--confirm-retries-finished",
                "--retire-source-embeddings",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((final / "FINAL_CACHE_COMPLETE.json").is_file())
        self.assertTrue((final / "contemporary_embedding_cache.tar.gz").is_file())
        self.assertFalse(baseline.exists())
        self.assertFalse((state / "cache").exists())
        roundtrip = self.root / "final_extract"
        extract_archive(
            final / "contemporary_embedding_cache.tar.gz", roundtrip, CONFIG
        )
        np.testing.assert_array_equal(
            np.load(roundtrip / "ppi/TEST.npy"), np.full(512, 2, dtype=np.float32)
        )


if __name__ == "__main__":
    unittest.main()
