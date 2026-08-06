from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "compare_go_relationship_policy.py"
SPEC = importlib.util.spec_from_file_location("go_relationship_policy", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _write_benchmark(root: Path, *, broad: bool) -> None:
    root.mkdir()
    terms = ["GO:0008150", "GO:0000001"] if broad else ["GO:0008150"]
    for aspect in MODULE.ASPECTS:
        root_term = MODULE.ROOTS[aspect]
        aspect_terms = [root_term, "GO:0000001"] if broad else [root_term]
        for split_name in MODULE.SPLITS:
            path = root / f"{aspect}-{split_name}.csv"
            path.write_text(
                "proteins,sequences," + ",".join(aspect_terms) + "\n"
                + "P1,MA," + ("1,1\n" if broad else "1\n"),
                encoding="ascii",
            )


class GoRelationshipPolicyTests(unittest.TestCase):
    def test_reports_final_label_and_term_differences(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            _write_benchmark(root / "broad", broad=True)
            _write_benchmark(root / "narrow", broad=False)
            obo = root / "go.obo"
            obo.write_text(
                "[Term]\nid: GO:0008150\nnamespace: biological_process\n"
                "[Term]\nid: GO:0000001\nnamespace: biological_process\n"
                "relationship: regulates GO:0008150\n",
                encoding="utf-8",
            )
            args = MODULE.argparse.Namespace(
                broad_label="broad", broad_dir=root / "broad",
                narrow_label="narrow", narrow_dir=root / "narrow",
                obo_file=obo, output_dir=root / "out",
            )
            payload = MODULE.run(args)
            row = payload["comparisons"]["bp"]["test"]
            self.assertEqual(row["common_proteins_with_changed_labels"], 1)
            self.assertEqual(row["broad_only_positive_labels_on_common_proteins"], 1)
            self.assertTrue((root / "out" / "RUN_COMPLETE.json").is_file())


if __name__ == "__main__":
    unittest.main()
