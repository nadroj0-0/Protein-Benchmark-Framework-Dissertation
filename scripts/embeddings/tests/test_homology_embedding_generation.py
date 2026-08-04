from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import io
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = ROOT / "scripts" / "embeddings"
sys.path.insert(0, str(SCRIPT_DIR))


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    assert specification.loader is not None
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


PREPARE = load_module(
    "prepare_homology_embedding_workspace_test",
    SCRIPT_DIR / "prepare_homology_embedding_workspace.py",
)
ASSEMBLE = load_module(
    "assemble_pair_resolved_embedding_cache_test",
    SCRIPT_DIR / "assemble_pair_resolved_embedding_cache.py",
)
RESOLVER = sys.modules["resolve_embedding_reuse_sources"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def npy_bytes(dimension: int, value: float) -> bytes:
    handle = io.BytesIO()
    np.save(handle, np.full(dimension, value, dtype=np.float32))
    return handle.getvalue()


def write_archive(path: Path, entries: dict[str, bytes]) -> None:
    with tarfile.open(path, "w:gz") as archive:
        for name, data in entries.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))


class HomologyEmbeddingGenerationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.benchmark = self.root / "benchmark"
        self.ledger = self.root / "ledger"
        self.benchmark.mkdir()
        self.ledger.mkdir()
        self.sequences = {"TRAIN": "AAAA", "VALID": "CCCC", "TEST": "GGGG"}
        for aspect in ("bp", "cc", "mf"):
            for split, protein_id in (
                ("training", "TRAIN"),
                ("validation", "VALID"),
                ("test", "TEST"),
            ):
                (self.benchmark / f"{aspect}-{split}.csv").write_text(
                    "proteins,sequences,GO:0000001\n"
                    f"{protein_id},{self.sequences[protein_id]},1\n",
                    encoding="utf-8",
                )
        self.memberships = {
            "TRAIN": [f"{aspect}-training.csv" for aspect in ("bp", "cc", "mf")],
            "VALID": [f"{aspect}-validation.csv" for aspect in ("bp", "cc", "mf")],
            "TEST": [f"{aspect}-test.csv" for aspect in ("bp", "cc", "mf")],
        }
        for names in self.memberships.values():
            names.sort()

        self.source = self.root / "source.tar.gz"
        source_entries = {}
        for protein_id, modalities in {
            "TRAIN": ASSEMBLE.MODALITIES,
            "TEST": ("sequence", "ppi"),
        }.items():
            for modality in modalities:
                directory = ASSEMBLE.MODALITY_DIRECTORIES[modality]
                source_entries[
                    f"data/embedding_cache/{directory}/{protein_id}.npy"
                ] = npy_bytes(ASSEMBLE.EXPECTED_DIMENSIONS[modality], 1.0)
        write_archive(self.source, source_entries)
        self._write_ledger(source_entries)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _action_row(self, protein_id: str, regenerate: list[str]) -> dict[str, str]:
        sequence = self.sequences[protein_id]
        return {
            "protein_id": protein_id,
            "sequence": sequence,
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "action": "regenerate" if regenerate else "reuse",
            "reason": "one-or-more-modalities-require-regeneration" if regenerate else "all-modalities-reusable",
            "reuse_modalities": json.dumps(
                sorted(set(ASSEMBLE.MODALITIES) - set(regenerate)), separators=(",", ":")
            ),
            "regenerate_modalities": json.dumps(sorted(regenerate), separators=(",", ":")),
            "coarse_action": "reuse",
            "target_memberships": json.dumps(self.memberships[protein_id], separators=(",", ":")),
        }

    def _write_ledger(self, source_entries: dict[str, bytes]) -> None:
        protein_columns = (
            "protein_id",
            "sequence",
            "sequence_sha256",
            "action",
            "reason",
            "reuse_modalities",
            "regenerate_modalities",
            "coarse_action",
            "target_memberships",
        )
        actions = {
            "TRAIN": [],
            "VALID": list(ASSEMBLE.MODALITIES),
            "TEST": ["structure", "text"],
        }
        for filename, action in (
            ("reuse_proteins.tsv", "reuse"),
            ("regenerate_proteins.tsv", "regenerate"),
        ):
            with (self.ledger / filename).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=protein_columns, delimiter="\t")
                writer.writeheader()
                for protein_id in sorted(actions):
                    row = self._action_row(protein_id, actions[protein_id])
                    if row["action"] == action:
                        writer.writerow(row)

        pair_path = self.ledger / "resolved_embedding_pairs.tsv.gz"
        with gzip.open(pair_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=RESOLVER.PAIR_COLUMNS, delimiter="\t")
            writer.writeheader()
            for protein_id in sorted(actions):
                for modality in ASSEMBLE.MODALITIES:
                    regenerate = modality in actions[protein_id]
                    member = (
                        f"data/embedding_cache/{ASSEMBLE.MODALITY_DIRECTORIES[modality]}/"
                        f"{protein_id}.npy"
                    )
                    data = source_entries.get(member)
                    row = {field: "" for field in RESOLVER.PAIR_COLUMNS}
                    row.update(
                        {
                            "protein_id": protein_id,
                            "modality": modality,
                            "action": "regenerate" if regenerate else "reuse",
                            "reason": "conflicting-source-arrays" if regenerate else "single-valid-source",
                            "sequence_sha256": hashlib.sha256(
                                self.sequences[protein_id].encode("ascii")
                            ).hexdigest(),
                            "coarse_action": "reuse",
                        }
                    )
                    if not regenerate:
                        array = np.load(io.BytesIO(data), allow_pickle=False)
                        row.update(
                            {
                                "selected_archive": str(self.source.resolve()),
                                "selected_member": member,
                                "file_sha256": hashlib.sha256(data).hexdigest(),
                                "array_sha256": RESOLVER.canonical_array_sha256(array),
                            }
                        )
                    writer.writerow(row)

        summary = {
            "comparison_policy": {"granularity": "protein-modality"},
            "sources": [
                {
                    "archive": str(self.source.resolve()),
                    "archive_sha256": sha256(self.source),
                }
            ],
        }
        (self.ledger / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
        (self.ledger / "run_manifest.json").write_text("{}\n", encoding="utf-8")
        files = []
        for path in sorted(self.ledger.iterdir()):
            files.append({"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
        manifest = {"schema_version": 1, "payload_file_count": len(files), "files": files}
        (self.ledger / "output_manifest.json").write_text(
            json.dumps(manifest, sort_keys=True), encoding="utf-8"
        )
        (self.ledger / "RUN_COMPLETE.json").write_text(
            json.dumps(
                {
                    "complete": True,
                    "output_manifest_sha256": sha256(self.ledger / "output_manifest.json"),
                }
            ),
            encoding="utf-8",
        )

    def _generated_archives(self) -> dict[str, Path]:
        actions = {
            "sequence": ["VALID"],
            "text": ["TEST", "VALID"],
            "structure": ["TEST", "VALID"],
            "ppi": ["VALID"],
        }
        result = {}
        for index, modality in enumerate(ASSEMBLE.MODALITIES, start=2):
            path = self.root / f"generated_{modality}.tar.gz"
            entries = {
                f"data/embedding_cache/{ASSEMBLE.MODALITY_DIRECTORIES[modality]}/{protein_id}.npy":
                npy_bytes(ASSEMBLE.EXPECTED_DIMENSIONS[modality], float(index))
                for protein_id in actions[modality]
            }
            write_archive(path, entries)
            result[modality] = path
        return result

    def _policy(self) -> Path:
        path = self.root / "policy.json"
        path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "modalities": {
                        modality: {
                            "cache_directory": ASSEMBLE.MODALITY_DIRECTORIES[modality],
                            "dimension": ASSEMBLE.EXPECTED_DIMENSIONS[modality],
                            "min_accepted_count": 1,
                        }
                        for modality in ASSEMBLE.MODALITIES
                    },
                }
            ),
            encoding="utf-8",
        )
        return path

    def test_workspace_selects_only_requested_modality_pairs(self) -> None:
        data = self.root / "data"
        report = self.root / "workspace.json"
        ledger_manifest = PREPARE.verify_ledger(self.ledger)
        sequences, memberships, _ = PREPARE.load_benchmark(self.benchmark)
        rows = PREPARE.load_plan_rows(self.ledger)
        PREPARE.validate_against_benchmark(rows, sequences, memberships)
        selected = PREPARE.select_rows(rows, "text", None)
        result = PREPARE.write_workspace(selected, data)
        self.assertEqual(ledger_manifest["payload_file_count"], 5)
        self.assertEqual(result["protein_count"], 2)
        fasta = (data / "proteins.fasta").read_text(encoding="ascii")
        self.assertNotIn("TRAIN", fasta)
        self.assertIn("VALID", fasta)
        self.assertIn("TEST", fasta)
        self.assertFalse(report.exists())

    def test_assembly_streams_reuse_and_generated_pairs_into_one_cache(self) -> None:
        output = self.root / "cache.tar.gz"
        reports = self.root / "reports"
        result = ASSEMBLE.publish_cache(
            self.ledger,
            self._generated_archives(),
            self._policy(),
            output,
            reports,
        )
        self.assertEqual(result["target_proteins"], 3)
        self.assertEqual(result["available_pairs"], 12)
        self.assertEqual(result["missing_pairs"], 0)
        with tarfile.open(output, "r:gz") as archive:
            files = [member for member in archive if member.isfile()]
        self.assertEqual(len(files), 12)
        self.assertTrue((reports / "RUN_COMPLETE.json").is_file())

    def test_assembly_accepts_hash_verified_local_source_override(self) -> None:
        original = self.source
        staged = self.root / "staged-source.tar.gz"
        staged.write_bytes(original.read_bytes())
        original.unlink()
        output = self.root / "overridden-cache.tar.gz"
        result = ASSEMBLE.publish_cache(
            self.ledger,
            self._generated_archives(),
            self._policy(),
            output,
            self.root / "overridden-reports",
            {str(original.resolve()): staged},
        )
        self.assertEqual(result["available_pairs"], 12)
        reuse_inputs = [item for item in result["inputs"] if item["role"] == "reuse"]
        self.assertEqual(reuse_inputs[0]["path"], str(original.resolve()))
        self.assertEqual(reuse_inputs[0]["read_path"], str(staged.resolve()))

    def test_assembly_rejects_tampered_local_source_override(self) -> None:
        staged = self.root / "tampered-source.tar.gz"
        staged.write_bytes(self.source.read_bytes() + b"tampered")
        with self.assertRaisesRegex(ASSEMBLE.AssemblyError, "Source archive hash mismatch"):
            ASSEMBLE.publish_cache(
                self.ledger,
                self._generated_archives(),
                self._policy(),
                self.root / "tampered-cache.tar.gz",
                self.root / "tampered-reports",
                {str(self.source.resolve()): staged},
            )

    def test_assembly_rejects_generated_pair_not_selected_by_ledger(self) -> None:
        generated = self._generated_archives()
        write_archive(
            generated["ppi"],
            {
                "data/embedding_cache/ppi/VALID.npy": npy_bytes(512, 4.0),
                "data/embedding_cache/ppi/TRAIN.npy": npy_bytes(512, 4.0),
            },
        )
        with self.assertRaisesRegex(ASSEMBLE.AssemblyError, "unrequested pair"):
            ASSEMBLE.publish_cache(
                self.ledger,
                generated,
                self._policy(),
                self.root / "bad.tar.gz",
                self.root / "bad-reports",
            )

    def test_hpc_wrappers_keep_gpu_and_ppi_resources_separate(self) -> None:
        gpu = (ROOT / "hpc_jobs/active/hpc_homology_embedding_gpu_array.sh").read_text()
        ppi = (ROOT / "hpc_jobs/active/hpc_homology_embedding_ppi.sh").read_text()
        finalizer = (
            ROOT / "hpc_jobs/active/hpc_homology_embedding_finalize.sh"
        ).read_text()
        self.assertIn("#$ -t 1-3", gpu)
        self.assertIn("#$ -pe gpu 1", gpu)
        self.assertIn("FRAMEWORK_JOB_ROOT", gpu)
        self.assertNotIn('$(dirname "$0")', gpu)
        self.assertIn("--modality ppi", ppi)
        self.assertNotIn("-l gpu=true", ppi)
        self.assertIn("FRAMEWORK_JOB_ROOT", ppi)
        common = (
            ROOT / "hpc_jobs/lib/run_homology_embedding_modality_job.sh"
        ).read_text()
        self.assertIn('cp -a "$LEDGER_DIR/." "$LEDGER_STAGE/"', common)
        self.assertIn('--ledger-dir "$LEDGER_STAGE"', common)
        self.assertIn('--benchmark-dir "$BENCHMARK_STAGE"', common)
        self.assertIn('--text-cutoff-date) TEXT_CUTOFF_DATE="$2"', common)
        self.assertIn('[[ "$MODALITY" == "text"', common)
        self.assertLess(
            finalizer.index('load_framework_paths "$FRAMEWORK_DIR"'),
            finalizer.index("activate_or_create_mmfp_env"),
        )

    def test_targeted_text_runner_uses_temporal_recipe_only_with_a_cutoff(self) -> None:
        runner = (SCRIPT_DIR / "run_homology_embedding_modality.sh").read_text()
        self.assertIn('--text-cutoff-date) TEXT_CUTOFF_DATE="$2"', runner)
        self.assertIn('bash "$HERE/generate_embeddings_text_temporal_cls.sh"', runner)
        self.assertIn('"text_cutoff_date": sys.argv[5] or None', runner)
        self.assertIn('scripts/extract_uniprot_text.py extract-current', runner)


if __name__ == "__main__":
    unittest.main()
