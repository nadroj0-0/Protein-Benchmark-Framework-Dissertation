import csv
from dataclasses import replace
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from pfp_embedding_inventory.benchmark import BenchmarkError, load_aliases, parse_benchmark  # noqa: E402
from pfp_embedding_inventory.inventory import InventoryError, build_inventory  # noqa: E402
from pfp_embedding_inventory.models import MODALITIES  # noqa: E402
from pfp_embedding_inventory.reports import write_reports  # noqa: E402

from helpers import DIMS, DIRS, make_cache, make_config, write_nine_csvs  # noqa: E402


class InventoryTests(unittest.TestCase):
    def _benchmarks(self, root: Path, target_rows=None, source_rows=None):
        target_dir = root / "target"
        source_dir = root / "source"
        write_nine_csvs(target_dir, shared_rows=target_rows or [("P1", "ACDE", "1")])
        write_nine_csvs(source_dir, shared_rows=source_rows or target_rows or [("P1", "ACDE", "1")])
        config = make_config()
        return (
            parse_benchmark(target_dir, config.benchmark_contract),
            parse_benchmark(source_dir, config.benchmark_contract),
            config,
        )

    def _records(self, result):
        return {(record.protein_id, record.modality): record for record in result.records}

    def test_all_modalities_present_and_valid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            make_cache(cache, ["P1"])
            result = build_inventory(target, source, cache, config, "paper-faithful")
            records = self._records(result)
            prott5 = records[("P1", "prott5")]
            self.assertEqual(prott5.factual_status, "present-valid")
            self.assertEqual(prott5.requested_action, "reuse")
            for modality in ("text", "structure", "ppi"):
                record = records[("P1", modality)]
                self.assertTrue(record.valid)
                self.assertEqual(record.factual_status, "provenance-unknown")
                self.assertEqual(record.requested_action, "manual-review")

    def test_no_modalities_and_policy_difference(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            paper = self._records(
                build_inventory(target, source, cache, config, "paper-faithful")
            )
            maximum = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage")
            )
            for modality in MODALITIES:
                self.assertEqual(paper[("P1", modality)].factual_status, "missing")
                self.assertEqual(paper[("P1", modality)].requested_action, "leave-masked")
            self.assertEqual(maximum[("P1", "prott5")].requested_action, "generate")
            self.assertEqual(maximum[("P1", "text")].requested_action, "generate")
            self.assertEqual(maximum[("P1", "structure")].requested_action, "unavailable")
            self.assertEqual(maximum[("P1", "ppi")].requested_action, "unavailable")

    def test_some_modalities_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            make_cache(cache, ["P1"], modalities=["prott5", "ppi"])
            records = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage")
            )
            self.assertEqual(records[("P1", "prott5")].requested_action, "reuse")
            self.assertEqual(records[("P1", "ppi")].requested_action, "manual-review")
            self.assertEqual(records[("P1", "text")].requested_action, "generate")
            self.assertEqual(records[("P1", "structure")].requested_action, "unavailable")

    def test_wrong_dimensions_nan_inf_and_unreadable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            for directory in DIRS.values():
                (cache / directory).mkdir(parents=True, exist_ok=True)
            np.save(cache / DIRS["prott5"] / "P1.npy", np.zeros(10, dtype=np.float32))
            text = np.zeros(DIMS["text"], dtype=np.float32)
            text[0] = np.nan
            np.save(cache / DIRS["text"] / "P1.npy", text)
            structure = np.zeros(DIMS["structure"], dtype=np.float32)
            structure[1] = np.inf
            np.save(cache / DIRS["structure"] / "P1.npy", structure)
            (cache / DIRS["ppi"] / "P1.npy").write_bytes(b"not a numpy file")
            records = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage")
            )
            self.assertEqual(records[("P1", "prott5")].factual_status, "wrong-dimension")
            self.assertEqual(records[("P1", "text")].factual_status, "non-finite")
            self.assertEqual(records[("P1", "structure")].factual_status, "non-finite")
            self.assertEqual(records[("P1", "ppi")].factual_status, "unreadable")
            self.assertEqual(records[("P1", "prott5")].requested_action, "generate")
            self.assertEqual(records[("P1", "text")].requested_action, "generate")
            self.assertEqual(records[("P1", "structure")].requested_action, "unavailable")
            self.assertEqual(records[("P1", "ppi")].requested_action, "unavailable")

    def test_sequence_mismatch_invalidates_sequence_modalities(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(
                root,
                target_rows=[("P1", "AAAA", "1")],
                source_rows=[("P1", "CCCC", "1")],
            )
            cache = root / "cache"
            make_cache(cache, ["P1"])
            records = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage")
            )
            self.assertEqual(records[("P1", "prott5")].factual_status, "sequence-mismatch")
            self.assertEqual(records[("P1", "structure")].factual_status, "sequence-mismatch")
            self.assertEqual(records[("P1", "prott5")].requested_action, "generate")
            self.assertEqual(records[("P1", "structure")].requested_action, "unavailable")
            self.assertEqual(records[("P1", "ppi")].requested_action, "manual-review")

    def test_exact_sequence_cross_id_reuse_is_prott5_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(
                root,
                target_rows=[("NEW", "ACDE", "1")],
                source_rows=[("OLD", "ACDE", "1")],
            )
            cache = root / "cache"
            make_cache(cache, ["OLD"])
            records = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage")
            )
            prott5 = records[("NEW", "prott5")]
            self.assertEqual(prott5.requested_action, "reuse")
            self.assertEqual(prott5.source_protein_id, "OLD")
            self.assertEqual(prott5.match_route, "sequence-sha256")
            for modality in ("text", "structure", "ppi"):
                self.assertEqual(records[("NEW", modality)].factual_status, "missing")

    def test_explicit_alias_and_ambiguous_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(
                root,
                target_rows=[("NEW", "ACDE", "1")],
                source_rows=[("OLD", "ACDE", "1"), ("OLD2", "ACDE", "1")],
            )
            cache = root / "cache"
            make_cache(cache, ["OLD", "OLD2"])
            alias_file = root / "aliases.tsv"
            alias_file.write_text(
                "protein_id\tsource_protein_id\tmodality\tmapping_route\tsource_identity\tmapping_evidence\n"
                "NEW\tOLD\ttext\tcurated-description-alias\ttest-text-v1\tdescription-sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa;temporal-context:test-text-v1\n"
                "NEW\tOLD\tstructure\tcurated-uniprot-alias\ttest-structure|v1\tstructure-source:test-structure;structure-version:v1\n"
                "NEW\tOLD\tppi\tcurated-string-node\ttest-ppi-v1\tstring-id:9606.ENSP1;string-release:test-ppi-v1\n",
                encoding="utf-8",
            )
            aliases = load_aliases(alias_file)
            records = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage", aliases)
            )
            structure = records[("NEW", "structure")]
            self.assertEqual(structure.requested_action, "reuse")
            self.assertEqual(structure.match_route, "explicit-alias:curated-uniprot-alias")
            self.assertEqual(records[("NEW", "text")].requested_action, "reuse")
            self.assertEqual(records[("NEW", "ppi")].requested_action, "reuse")

            alias_file.write_text(
                "protein_id\tsource_protein_id\tmodality\tmapping_route\tsource_identity\tmapping_evidence\n"
                "NEW\tOLD\tppi\tmap-a\ttest-ppi-v1\tstring-id:node-a;string-release:test-ppi-v1\n"
                "NEW\tOLD2\tppi\tmap-b\ttest-ppi-v1\tstring-id:node-b;string-release:test-ppi-v1\n",
                encoding="utf-8",
            )
            ambiguous = self._records(
                build_inventory(
                    target,
                    source,
                    cache,
                    config,
                    "maximize-coverage",
                    load_aliases(alias_file),
                )
            )[("NEW", "ppi")]
            self.assertEqual(ambiguous.factual_status, "provenance-unknown")
            self.assertEqual(ambiguous.requested_action, "manual-review")
            self.assertEqual(ambiguous.match_route, "ambiguous-explicit-alias")

    def test_alias_source_cannot_escape_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aliases.tsv"
            path.write_text(
                "protein_id\tsource_protein_id\tmodality\tmapping_route\tsource_identity\tmapping_evidence\n"
                "P1\t../outside\tppi\tbad-route\ttest-ppi-v1\tunsafe\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(BenchmarkError, "unsafe or malformed"):
                load_aliases(path)

    def test_alias_release_evidence_requires_exact_identity_token(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            make_cache(cache, ["P1"], modalities=["ppi"])
            alias_file = root / "aliases.tsv"
            alias_file.write_text(
                "protein_id\tsource_protein_id\tmodality\tmapping_route\tsource_identity\tmapping_evidence\n"
                "P1\tP1\tppi\tcurated-string-node\ttest-ppi-v1\tstring-id:node;string-release:v1\n",
                encoding="utf-8",
            )
            record = self._records(
                build_inventory(
                    target,
                    source,
                    cache,
                    config,
                    "maximize-coverage",
                    load_aliases(alias_file),
                )
            )[("P1", "ppi")]
            self.assertEqual(record.factual_status, "provenance-unknown")
            self.assertEqual(record.requested_action, "manual-review")

    def test_unused_alias_target_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            alias_file = root / "aliases.tsv"
            alias_file.write_text(
                "protein_id\tsource_protein_id\tmodality\tmapping_route\tsource_identity\tmapping_evidence\n"
                "TYPO\tP1\tppi\tcurated\ttest-ppi-v1\tSTRING-v12:node\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(InventoryError, "absent from the benchmark"):
                build_inventory(
                    target,
                    source,
                    cache,
                    config,
                    "maximize-coverage",
                    load_aliases(alias_file),
                )

    def test_exact_id_precedes_prott5_alias_and_valid_hash_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(
                root,
                target_rows=[("P1", "ACDE", "1")],
                source_rows=[("P1", "ACDE", "1"), ("ALT", "ACDE", "1")],
            )
            cache = root / "cache"
            make_cache(cache, ["P1", "ALT"], modalities=["prott5"])
            alias_file = root / "aliases.tsv"
            alias_file.write_text(
                "protein_id\tsource_protein_id\tmodality\tmapping_route\tsource_identity\tmapping_evidence\n"
                "P1\tALT\tprott5\tcurated-alias\ttest-prott5-v1\tsequence-sha256:ignored-because-exact-id-is-primary\n",
                encoding="utf-8",
            )
            record = self._records(
                build_inventory(
                    target,
                    source,
                    cache,
                    config,
                    "maximize-coverage",
                    load_aliases(alias_file),
                )
            )[("P1", "prott5")]
            self.assertEqual(record.match_route, "exact-id")
            self.assertEqual(record.source_protein_id, "P1")

            np.save(cache / DIRS["prott5"] / "P1.npy", np.zeros(10, dtype=np.float32))
            record = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage")
            )[("P1", "prott5")]
            self.assertEqual(record.source_protein_id, "ALT")
            self.assertEqual(record.match_route, "sequence-sha256")

    def test_unsupported_integer_dtype_is_invalid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            directory = cache / DIRS["prott5"]
            directory.mkdir(parents=True)
            np.save(directory / "P1.npy", np.zeros(DIMS["prott5"], dtype=np.int32))
            record = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage")
            )[("P1", "prott5")]
            self.assertEqual(record.factual_status, "unsupported-dtype")
            self.assertFalse(record.valid)
            self.assertEqual(record.requested_action, "generate")

    def test_unknown_text_and_structure_provenance_routes_to_manual_review(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, _ = self._benchmarks(root)
            config = make_config({"text": "unknown", "structure": "unknown"})
            cache = root / "cache"
            make_cache(cache, ["P1"])
            records = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage")
            )
            for modality in ("text", "structure"):
                self.assertEqual(records[("P1", modality)].factual_status, "provenance-unknown")
                self.assertEqual(records[("P1", modality)].requested_action, "manual-review")

    def test_mixed_temporal_text_role_is_never_compatible(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            text_spec = config.modalities["text"]
            text_provenance = replace(
                text_spec.provenance,
                text_role_policy="cafa3-mixed-temporal",
            )
            modalities = dict(config.modalities)
            modalities["text"] = replace(text_spec, provenance=text_provenance)
            config = replace(config, modalities=modalities)
            cache = root / "cache"
            make_cache(cache, ["P1"], modalities=["text"])
            alias_file = root / "aliases.tsv"
            alias_file.write_text(
                "protein_id\tsource_protein_id\tmodality\tmapping_route\tsource_identity\tmapping_evidence\n"
                "P1\tP1\ttext\tcurated-description\ttest-text-v1\tdescription-sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa;temporal-context:test-text-v1\n",
                encoding="utf-8",
            )
            record = self._records(
                build_inventory(
                    target,
                    source,
                    cache,
                    config,
                    "maximize-coverage",
                    load_aliases(alias_file),
                )
            )[("P1", "text")]
            self.assertEqual(record.factual_status, "provenance-incompatible")
            self.assertEqual(record.requested_action, "manual-review")

    def test_reports_include_extras_and_all_required_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            make_cache(cache, ["P1", "EXTRA"])
            result = build_inventory(target, source, cache, config, "paper-faithful")
            output = root / "reports"
            summary = write_reports(result, output, cache)
            required = {
                "benchmark_proteins.tsv",
                "protein_embedding_summary.tsv",
                "embedding_inventory.tsv",
                "embedding_summary.json",
                "embedding_summary.md",
                "reuse_manifest.tsv",
                "generation_manifest.tsv",
                "manual_review.tsv",
                "reuse_prott5.txt",
                "generate_prott5.fasta",
                "reuse_text.txt",
                "generate_text.txt",
                "text_manual_review.txt",
                "reuse_structure.txt",
                "generate_structure.txt",
                "structure_unavailable.txt",
                "reuse_ppi.txt",
                "extract_ppi.txt",
                "ppi_unavailable.txt",
                "cache_extras.tsv",
            }
            self.assertTrue(required.issubset({path.name for path in output.iterdir()}))
            self.assertEqual(summary["cache_extras"], {modality: 1 for modality in MODALITIES})
            self.assertEqual(summary["coverage"]["global"]["at_least_one_modality"]["count"], 1)
            with (output / "protein_embedding_summary.tsv").open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["has_any_embedding"], "true")
            self.assertEqual(row["reusable_modalities"], "prott5")


if __name__ == "__main__":
    unittest.main()
