import csv
import os
import sys
import unittest
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set, Tuple


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PACKAGE_ROOT.parent
SRC = PACKAGE_ROOT / "src"
sys.path.insert(0, str(SRC))

from pfp_embedding_inventory.benchmark import parse_benchmark  # noqa: E402
from pfp_embedding_inventory.config import load_config  # noqa: E402
from pfp_embedding_inventory.inventory import ArrayCache, build_inventory  # noqa: E402
from pfp_embedding_inventory.models import MODALITIES, ONTOLOGIES, SPLITS  # noqa: E402
from pfp_embedding_inventory.reports import write_reports  # noqa: E402


CANONICAL = Path(
    "/Users/jordansydney-darlin/CAFA3 Supplementary material/benchmark_inputs/"
    "cafa3_csv/cafa3_raw"
)
CACHE = Path(
    "/Users/jordansydney-darlin/CAFA3 Supplementary material/external_repos/"
    "PFP_reference_clone/data/embedding_cache"
)
HISTORICAL = (
    Path(
        "/Users/jordansydney-darlin/CAFA3 Supplementary material/benchmark_results/"
        "validation/cafa3_historical_validation/7061973_20260712_193237/generated"
    ),
    Path(
        "/Users/jordansydney-darlin/CAFA3 Supplementary material/benchmark_results/"
        "validation/cafa3_historical_validation/7061922_20260712_161031/generated"
    ),
)
CONFIG = REPO_ROOT / "configs" / "embedding_inventory.cafa3_published.json"
DIRS = {
    "prott5": "prott5",
    "text": "exp_text_embeddings_temporal",
    "structure": "IF1",
    "ppi": "ppi",
}


@unittest.skipUnless(
    os.environ.get("PFP_RUN_REAL_INTEGRATION") == "1",
    "set PFP_RUN_REAL_INTEGRATION=1 to run local CAFA3 cache integration",
)
class RealDataIntegrationTests(unittest.TestCase):
    def test_canonical_golden_and_historical_variants(self):
        for path in (CANONICAL, CACHE) + HISTORICAL:
            self.assertTrue(path.exists(), str(path))
        config = load_config(CONFIG)
        source = parse_benchmark(CANONICAL, config.benchmark_contract)
        independent_source = independent_parse(CANONICAL)
        self.assertEqual(set(source.proteins), set(independent_source["sequences"]))
        self.assertEqual(len(source.proteins), 69811)

        output_root = Path(
            os.environ.get(
                "PFP_INVENTORY_REAL_OUTPUT_ROOT",
                str(REPO_ROOT / "results" / "embedding_inventory" / "real_integration"),
            )
        )
        array_cache: ArrayCache = {}
        cache_ids = {
            modality: {path.stem for path in (CACHE / directory).glob("*.npy")}
            for modality, directory in DIRS.items()
        }

        canonical_result = build_inventory(
            source, source, CACHE, config, "paper-faithful", array_cache=array_cache
        )
        canonical_summary = write_reports(
            canonical_result, output_root / "canonical_paper_faithful", CACHE
        )
        self._assert_golden(canonical_result, canonical_summary, independent_source, cache_ids)

        for historical_dir in HISTORICAL:
            with self.subTest(historical=historical_dir.parent.name):
                target = parse_benchmark(historical_dir, config.benchmark_contract)
                independent_target = independent_parse(historical_dir)
                result = build_inventory(
                    target,
                    source,
                    CACHE,
                    config,
                    "maximize-coverage",
                    array_cache=array_cache,
                )
                summary = write_reports(
                    result,
                    output_root / (historical_dir.parent.name + "_maximize_coverage"),
                    CACHE,
                )
                self._assert_historical(
                    result,
                    summary,
                    independent_source,
                    independent_target,
                    cache_ids,
                )

    def _assert_golden(self, result, summary, independent, cache_ids):
        records = {(record.protein_id, record.modality): record for record in result.records}
        self.assertEqual(len(records), 4 * len(independent["sequences"]))
        self.assertFalse(any(record.requested_action == "generate" for record in result.records))
        for modality in MODALITIES:
            direct = sum(protein_id in cache_ids[modality] for protein_id in independent["sequences"])
            inventory_present = summary["coverage"]["global"]["by_modality"][modality]["present"]["count"]
            self.assertEqual(inventory_present, direct)
            self.assertEqual(len(cache_ids[modality]), direct)
            for protein_id in independent["sequences"]:
                record = records[(protein_id, modality)]
                if protein_id in cache_ids[modality]:
                    self.assertTrue(record.valid)
                    self.assertTrue(record.finite)
                    if modality == "prott5":
                        self.assertEqual(record.factual_status, "present-valid")
                        self.assertEqual(record.requested_action, "reuse")
                    else:
                        self.assertEqual(record.factual_status, "provenance-unknown")
                        self.assertEqual(record.requested_action, "manual-review")
                else:
                    self.assertEqual(record.factual_status, "missing")
                    self.assertEqual(record.requested_action, "leave-masked")

        for split, ids in independent["by_split"].items():
            for modality in MODALITIES:
                expected = sum(protein_id in cache_ids[modality] for protein_id in ids)
                observed = summary["coverage"]["by_split"][split]["by_modality"][modality]["present"]["count"]
                self.assertEqual(observed, expected)
        for ontology, ids in independent["by_ontology"].items():
            for modality in MODALITIES:
                expected = sum(protein_id in cache_ids[modality] for protein_id in ids)
                observed = summary["coverage"]["by_ontology"][ontology]["by_modality"][modality]["present"]["count"]
                self.assertEqual(observed, expected)

        all_ids = set(independent["sequences"])
        complete = sum(all(protein_id in cache_ids[m] for m in MODALITIES) for protein_id in all_ids)
        at_least_one = sum(any(protein_id in cache_ids[m] for m in MODALITIES) for protein_id in all_ids)
        self.assertEqual(summary["coverage"]["global"]["complete_four_modalities_present"]["count"], complete)
        self.assertEqual(summary["coverage"]["global"]["at_least_one_modality"]["count"], at_least_one)
        self.assertEqual(
            summary["coverage"]["global"]["complete_four_modalities_reusable"]["count"],
            0,
        )
        self.assertIn("zero vectors with mask 0.0", summary["pfp_missing_behavior"])

    def _assert_historical(self, result, summary, source, target, cache_ids):
        source_ids = set(source["sequences"])
        target_ids = set(target["sequences"])
        shared = source_ids & target_ids
        changed = {
            protein_id
            for protein_id in shared
            if source["sequences"][protein_id] != target["sequences"][protein_id]
        }
        new = target_ids - source_ids
        removed = source_ids - target_ids
        records = {(record.protein_id, record.modality): record for record in result.records}
        source_by_sequence = defaultdict(set)
        for protein_id, sequence in source["sequences"].items():
            source_by_sequence[sequence].add(protein_id)

        expected_prott5_reuse = {
            protein_id
            for protein_id, sequence in target["sequences"].items()
            if any(source_id in cache_ids["prott5"] for source_id in source_by_sequence.get(sequence, set()))
        }
        actual_prott5_reuse = {
            record.protein_id
            for record in result.records
            if record.modality == "prott5" and record.requested_action == "reuse"
        }
        self.assertEqual(actual_prott5_reuse, expected_prott5_reuse)
        for protein_id in new:
            record = records[(protein_id, "prott5")]
            if target["sequences"][protein_id] in source_by_sequence:
                self.assertEqual(record.match_route, "sequence-sha256")
            else:
                self.assertEqual(record.factual_status, "missing")
                self.assertEqual(record.requested_action, "generate")

        for protein_id in changed:
            self.assertNotEqual(records[(protein_id, "prott5")].match_route, "exact-id")
            self.assertNotEqual(records[(protein_id, "structure")].requested_action, "reuse")

        for modality in ("text", "structure", "ppi"):
            modality_records = [
                record for record in result.records if record.modality == modality
            ]
            self.assertFalse(
                any(record.requested_action == "reuse" for record in modality_records)
            )
            self.assertFalse(
                any(record.requested_action == "generate" for record in modality_records)
            )
            self.assertTrue(
                all(record.requested_action == "manual-review" for record in modality_records)
            )
            expected_present = len(target_ids & cache_ids[modality])
            observed_present = summary["coverage"]["global"]["by_modality"][modality][
                "present"
            ]["count"]
            self.assertEqual(observed_present, expected_present)

        for modality in MODALITIES:
            expected_extras = cache_ids[modality] - target_ids
            self.assertEqual(summary["cache_extras"][modality], len(expected_extras))
            self.assertTrue((removed & cache_ids[modality]).issubset(expected_extras))

        for record in result.records:
            if record.requested_action == "reuse":
                self.assertIn(
                    record.match_route,
                    {"exact-id", "sequence-sha256"},
                )


def independent_parse(directory: Path) -> Dict[str, object]:
    sequences: Dict[str, str] = {}
    memberships: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    by_split: Dict[str, Set[str]] = {split: set() for split in SPLITS}
    by_ontology: Dict[str, Set[str]] = {ontology: set() for ontology in ONTOLOGIES}
    for ontology in ONTOLOGIES:
        for split in SPLITS:
            path = directory / ("%s-%s.csv" % (ontology.lower(), split))
            with path.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.reader(handle)
                next(reader)
                for row in reader:
                    protein_id, sequence = row[:2]
                    if protein_id in sequences:
                        if sequences[protein_id] != sequence:
                            raise AssertionError("independent parser found sequence conflict")
                    else:
                        sequences[protein_id] = sequence
                    memberships[protein_id].add((ontology, split))
                    by_split[split].add(protein_id)
                    by_ontology[ontology].add(protein_id)
    return {
        "sequences": sequences,
        "memberships": memberships,
        "by_split": by_split,
        "by_ontology": by_ontology,
    }


if __name__ == "__main__":
    unittest.main()
