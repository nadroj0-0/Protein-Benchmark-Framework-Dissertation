from __future__ import annotations

import csv
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


DIAGNOSTICS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DIAGNOSTICS))

import analyse_homology_plateau_mechanism as module  # noqa: E402


class HomologyPlateauMechanismTests(unittest.TestCase):
    def test_size_bins_have_explicit_boundaries(self) -> None:
        self.assertEqual([module.size_bin(value) for value in (0, 1, 2, 3, 6, 11)], [
            "0", "1", "2", "3-5", "6-10", ">10"
        ])
        self.assertEqual(module.size_bin(6, exposure=True), "6+")

    def test_assignment_aliases_do_not_inflate_canonical_counts(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "assignments.tsv"
            fields = [
                "raw_uniprot_accession", "uniprot_accession", "uniref50_id",
                "mapping_status", "mmseqs_cluster_id", "split",
            ]
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow({
                    "raw_uniprot_accession": "OLD1", "uniprot_accession": "P1",
                    "uniref50_id": "UniRef50_P1",
                    "mapping_status": "mapped",
                    "mmseqs_cluster_id": "C1", "split": "training",
                })
            aliases, canonical, uniref50 = module.read_assignments(
                path, {"C1": {"split": "training", "full": 2, "supervised": 1}}
            )
            self.assertEqual(aliases["OLD1"], aliases["P1"])
            self.assertEqual(canonical, {"C1": 1})
            self.assertEqual(uniref50, {"C1": 1})

    def test_macro_metric_matches_cafa_definition(self) -> None:
        truth = np.asarray([[1, 0], [1, 1]], dtype=np.uint8)
        scores = np.asarray([[0.9, 0.8], [0.9, 0.2]], dtype=np.float32)
        components = module.metric_components(truth, scores, 0.5)
        result = module.aggregate_metric(components, np.asarray([0, 1]))
        self.assertAlmostEqual(result["precision"], 0.75)
        self.assertAlmostEqual(result["recall"], 0.75)
        self.assertAlmostEqual(result["f"], 0.75)

    def test_cluster_bootstrap_is_deterministic(self) -> None:
        truth = np.asarray([[1, 0], [1, 1], [0, 1]], dtype=np.uint8)
        left = module.metric_components(truth, np.asarray([[.9, .2], [.8, .7], [.2, .8]]), 0.5)
        right = module.metric_components(truth, np.asarray([[.7, .6], [.8, .2], [.4, .6]]), 0.5)
        indices = np.arange(3)
        first = module.bootstrap_delta(left, right, ["A", "A", "B"], indices, 50, 42)
        second = module.bootstrap_delta(left, right, ["A", "A", "B"], indices, 50, 42)
        self.assertEqual(first, second)

    def test_chunked_components_apply_writer_rounding_and_parent_propagation(self) -> None:
        class Graph:
            def ancestor_closure(self, aspect: str) -> dict[str, set[str]]:
                self.assert_aspect = aspect
                return {
                    "GO:ROOT": {"GO:ROOT"},
                    "GO:CHILD": {"GO:ROOT", "GO:CHILD"},
                }

        graph = Graph()
        terms = ("GO:ROOT", "GO:CHILD")
        plan = module.propagation_plan(terms, graph, "BPO")
        truth = np.asarray([[1, 1], [1, 0]], dtype=np.uint8)
        scores = np.asarray([[0.1, 0.9000004], [0.4999996, 0.2]], dtype=np.float32)
        result = module.chunked_metric_components(
            truth, scores, np.asarray([0, 1]), (0.5,), plan, chunk_size=1
        )[0.5]
        self.assertEqual(graph.assert_aspect, "BPO")
        np.testing.assert_allclose(result["precision"], [1.0, 1.0])
        np.testing.assert_allclose(result["recall"], [1.0, 1.0])
        np.testing.assert_array_equal(result["covered"], [True, True])


if __name__ == "__main__":
    unittest.main()
