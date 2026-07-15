import csv
from dataclasses import replace
import json
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
from pfp_embedding_inventory.reports import ReportError, write_reports  # noqa: E402

from helpers import (  # noqa: E402
    DIMS, DIRS, bind_verified_cache, make_cache, make_config, write_nine_csvs,
)


class InventoryTests(unittest.TestCase):
    def _benchmarks(self, root: Path, target_rows=None, source_rows=None):
        target_dir = root / "target"
        source_dir = root / "source"
        write_nine_csvs(target_dir, shared_rows=target_rows or [("P1", "ACDE", "1")])
        write_nine_csvs(source_dir, shared_rows=source_rows or target_rows or [("P1", "ACDE", "1")])
        config = make_config()
        return (
            parse_benchmark(target_dir, config.target_benchmark_contract),
            parse_benchmark(source_dir, config.source_benchmark_contract),
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
            config, proof = bind_verified_cache(config, source, cache)
            result = build_inventory(
                target, source, cache, config, "paper-faithful",
                artifact_verification=proof,
            )
            records = self._records(result)
            prott5 = records[("P1", "prott5")]
            self.assertEqual(prott5.factual_status, "present-valid")
            self.assertEqual(prott5.requested_action, "reuse")
            for modality in ("text", "structure", "ppi"):
                record = records[("P1", modality)]
                self.assertTrue(record.valid)
                self.assertEqual(record.factual_status, "provenance-unknown")
                self.assertEqual(record.requested_action, "regenerate")

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
                self.assertEqual(paper[("P1", modality)].requested_action, "regenerate")
                self.assertEqual(maximum[("P1", modality)].requested_action, "regenerate")

    def test_some_modalities_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            make_cache(cache, ["P1"], modalities=["prott5", "ppi"])
            config, proof = bind_verified_cache(config, source, cache)
            records = self._records(
                build_inventory(
                    target, source, cache, config, "maximize-coverage",
                    artifact_verification=proof,
                )
            )
            self.assertEqual(records[("P1", "prott5")].requested_action, "reuse")
            self.assertEqual(records[("P1", "ppi")].requested_action, "regenerate")
            self.assertEqual(records[("P1", "text")].requested_action, "regenerate")
            self.assertEqual(records[("P1", "structure")].requested_action, "regenerate")

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
            self.assertEqual(records[("P1", "prott5")].requested_action, "regenerate")
            self.assertEqual(records[("P1", "text")].requested_action, "regenerate")
            self.assertEqual(records[("P1", "structure")].requested_action, "regenerate")
            self.assertEqual(records[("P1", "ppi")].requested_action, "regenerate")

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
            self.assertEqual(records[("P1", "prott5")].requested_action, "regenerate")
            self.assertEqual(records[("P1", "structure")].requested_action, "regenerate")
            self.assertEqual(records[("P1", "ppi")].requested_action, "regenerate")

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
            config, proof = bind_verified_cache(config, source, cache)
            records = self._records(
                build_inventory(
                    target, source, cache, config, "maximize-coverage",
                    artifact_verification=proof,
                )
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
            config, proof = bind_verified_cache(config, source, cache)
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
                build_inventory(
                    target, source, cache, config, "maximize-coverage", aliases,
                    artifact_verification=proof,
                )
            )
            structure = records[("NEW", "structure")]
            self.assertEqual(structure.requested_action, "regenerate")
            self.assertEqual(structure.match_route, "explicit-alias:curated-uniprot-alias")
            self.assertEqual(records[("NEW", "text")].requested_action, "regenerate")
            self.assertEqual(records[("NEW", "ppi")].requested_action, "regenerate")
            for modality in ("text", "structure", "ppi"):
                self.assertIn(
                    "reuse is unsupported",
                    records[("NEW", modality)].reason,
                )

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
                    artifact_verification=proof,
                )
            )[("NEW", "ppi")]
            self.assertEqual(ambiguous.factual_status, "provenance-unknown")
            self.assertEqual(ambiguous.requested_action, "regenerate")
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

    def test_embedding_file_symlink_cannot_escape_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            directory = cache / DIRS["prott5"]
            directory.mkdir(parents=True)
            outside = root / "outside.npy"
            np.save(outside, np.zeros(DIMS["prott5"], dtype=np.float32))
            (directory / "P1.npy").symlink_to(outside)
            record = self._records(
                build_inventory(target, source, cache, config, "maximize-coverage")
            )[("P1", "prott5")]
            self.assertEqual(record.factual_status, "unreadable")
            self.assertNotEqual(record.requested_action, "reuse")

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
            self.assertEqual(record.requested_action, "regenerate")

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
            self.assertEqual(record.requested_action, "regenerate")

    def test_unknown_text_and_structure_provenance_routes_to_regeneration(self):
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
                self.assertEqual(records[("P1", modality)].requested_action, "regenerate")

    def test_fixed_ppi_direct_id_is_reusable_but_unproven_modalities_regenerate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            modalities = dict(config.modalities)
            ppi = modalities["ppi"]
            modalities["ppi"] = replace(
                ppi,
                provenance=replace(ppi.provenance, allow_direct_id_reuse=True),
            )
            config = replace(config, modalities=modalities)
            cache = root / "cache"
            make_cache(cache, ["P1"])
            config, proof = bind_verified_cache(config, source, cache)
            records = self._records(
                build_inventory(
                    target, source, cache, config, "maximize-coverage",
                    artifact_verification=proof,
                )
            )
            self.assertEqual(records[("P1", "ppi")].requested_action, "reuse")
            self.assertTrue(records[("P1", "ppi")].scientifically_eligible)

            unknown_config = make_config({"text": "unknown", "structure": "unknown"})
            unknown_records = self._records(
                build_inventory(target, source, cache, unknown_config, "maximize-coverage")
            )
            for modality in ("text", "structure"):
                self.assertEqual(unknown_records[("P1", modality)].requested_action, "regenerate")
                self.assertIn(
                    "reuse is not positively proven",
                    unknown_records[("P1", modality)].reason,
                )

    def test_ppi_alias_cannot_displace_valid_direct_id_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(
                root,
                target_rows=[("P1", "ACDE", "1")],
                source_rows=[("P1", "ACDE", "1"), ("OLD", "AAAA", "1")],
            )
            modalities = dict(config.modalities)
            ppi = modalities["ppi"]
            modalities["ppi"] = replace(
                ppi,
                provenance=replace(ppi.provenance, allow_direct_id_reuse=True),
            )
            config = replace(config, modalities=modalities)
            cache = root / "cache"
            make_cache(cache, ["P1", "OLD"])
            config, proof = bind_verified_cache(config, source, cache)
            alias_file = root / "aliases.tsv"
            alias_file.write_text(
                "protein_id\tsource_protein_id\tmodality\tmapping_route\t"
                "source_identity\tmapping_evidence\n"
                "P1\tOLD\tppi\tlegacy-alias\ttest-ppi-v1\t"
                "string-id:node-old;string-release:test-ppi-v1\n",
                encoding="utf-8",
            )
            aliases = load_aliases(alias_file)
            for policy in ("paper-faithful", "maximize-coverage"):
                with self.subTest(policy=policy):
                    record = self._records(
                        build_inventory(
                            target, source, cache, config, policy, aliases,
                            artifact_verification=proof,
                        )
                    )[("P1", "ppi")]
                    self.assertEqual(record.requested_action, "reuse")
                    self.assertEqual(record.source_protein_id, "P1")
                    self.assertEqual(record.match_route, "exact-id")

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
            self.assertEqual(record.requested_action, "regenerate")

    def test_reports_include_extras_and_all_required_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            make_cache(cache, ["P1", "EXTRA"])
            config, proof = bind_verified_cache(config, source, cache)
            result = build_inventory(
                target, source, cache, config, "paper-faithful",
                artifact_verification=proof,
            )
            output = root / "reports"
            summary = write_reports(result, output, cache)
            required = {
                "benchmark_proteins.tsv.gz",
                "protein_embedding_summary.tsv.gz",
                "embedding_summary.json",
                "embedding_summary.md",
                "run_provenance.json",
                "run_provenance.md",
                "cache_extras.tsv.gz",
                "reuse.tsv",
                "regenerate.tsv",
                "regenerate_reasons.tsv",
            }
            self.assertTrue(required.issubset({path.name for path in output.iterdir()}))
            self.assertTrue((output / "reuse" / "prott5.txt").exists())
            self.assertTrue((output / "regenerate" / "prott5.fasta").exists())
            self.assertTrue((output / "regenerate" / "text.txt").exists())
            self.assertEqual(summary["cache_extras"], {modality: 1 for modality in MODALITIES})
            self.assertEqual(summary["coverage"]["global"]["at_least_one_modality"]["count"], 1)
            self.assertEqual(
                summary["coverage"]["global"]["requested_action_counts"],
                {"regenerate": 3, "reuse": 1},
            )
            import gzip
            with gzip.open(output / "protein_embedding_summary.tsv.gz", "rt") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["has_any_embedding"], "true")
            self.assertEqual(row["reuse_modalities"], "prott5")
            self.assertEqual(row["regenerate_modalities"], "text;structure;ppi")

            binary_rows = []
            rows_by_action = {}
            for name, action in (("reuse.tsv", "reuse"), ("regenerate.tsv", "regenerate")):
                with (output / name).open(newline="", encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertTrue(all(row["requested_action"] == action for row in rows))
                rows_by_action[action] = rows
                binary_rows.extend(rows)
            pairs = [(row["protein_id"], row["modality"]) for row in binary_rows]
            self.assertEqual(len(pairs), len(set(pairs)))
            self.assertEqual(set(pairs), {("P1", modality) for modality in MODALITIES})
            self.assertTrue(all(row["reason"] for row in binary_rows))
            for action in ("reuse", "regenerate"):
                for modality in MODALITIES:
                    expected_ids = sorted(
                        row["protein_id"]
                        for row in rows_by_action[action]
                        if row["modality"] == modality
                    )
                    actual_ids = (output / action / (modality + ".txt")).read_text(
                        encoding="utf-8"
                    ).splitlines()
                    self.assertEqual(actual_ids, expected_ids)
            fasta_ids = [
                line[1:]
                for line in (output / "regenerate" / "prott5.fasta").read_text(
                    encoding="utf-8"
                ).splitlines()
                if line.startswith(">")
            ]
            self.assertEqual(
                fasta_ids,
                sorted(
                    row["protein_id"]
                    for row in rows_by_action["regenerate"]
                    if row["modality"] == "prott5"
                ),
            )
            manifest = json.loads((output / "output_manifest.json").read_text())
            manifested_paths = {item["path"] for item in manifest["files"]}
            self.assertIn("reuse/prott5.txt", manifested_paths)
            self.assertIn("regenerate/prott5.fasta", manifested_paths)

    def test_reports_fail_closed_on_incomplete_binary_partition(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            result = build_inventory(target, source, cache, config, "maximize-coverage")
            result.records.pop()
            output = root / "reports"
            with self.assertRaisesRegex(ReportError, "binary reuse partition is invalid"):
                write_reports(result, output, cache)
            self.assertFalse(output.exists())

    def test_reports_reject_any_third_action(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, source, config = self._benchmarks(root)
            cache = root / "cache"
            result = build_inventory(target, source, cache, config, "maximize-coverage")
            with self.assertRaisesRegex(ValueError, "exactly 'reuse' or 'regenerate'"):
                replace(result.records[0], requested_action="third-state")


if __name__ == "__main__":
    unittest.main()
