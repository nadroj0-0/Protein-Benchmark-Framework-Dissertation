from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cafa_benchmark_builder.ontology import Ontology


OBO = """format-version: 1.2

[Term]
id: GO:0008150
name: root
namespace: biological_process

[Term]
id: GO:0000001
name: parent
namespace: biological_process
is_a: GO:0008150 ! root

[Term]
id: GO:0000002
name: child
namespace: biological_process
is_a: GO:0000001 ! parent
relationship: regulates GO:0000003 ! regulated

[Term]
id: GO:0000003
name: regulated
namespace: biological_process
is_a: GO:0008150 ! root
"""


class RelationshipPolicyTests(unittest.TestCase):
    def test_default_remains_all_relationships_and_filter_is_opt_in(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "go.obo"
            path.write_text(OBO, encoding="utf-8")
            broad = Ontology(path, with_rels=True)
            narrow = Ontology(
                path, with_rels=True, relationship_types=frozenset({"part_of"})
            )
            self.assertIn("GO:0000003", broad.get_ancestors("GO:0000002"))
            self.assertNotIn("GO:0000003", narrow.get_ancestors("GO:0000002"))
            self.assertEqual(
                narrow.get_ancestors("GO:0000002"),
                {"GO:0000002", "GO:0000001", "GO:0008150"},
            )


if __name__ == "__main__":
    unittest.main()
