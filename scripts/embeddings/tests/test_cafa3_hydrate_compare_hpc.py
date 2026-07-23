from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "hpc_jobs" / "active" / "hpc_cafa3_hydrate_compare.sh"


class Cafa3HydrateCompareHpcTests(unittest.TestCase):
    def test_workflow_is_non_destructive_and_archived(self) -> None:
        text = SCRIPT.read_text(encoding="utf-8")

        self.assertIn("--preserve-evidence", text)
        self.assertIn("manage_embedding_archive.py create", text)
        self.assertIn("manage_embedding_archive.py extract", text)
        self.assertIn("compare_embeddings.py", text)
        self.assertIn("source_state_unchanged", text)
        self.assertIn("$1 != $3 || $2 != $4", text)
        self.assertIn("output_manifest.tsv", text)
        self.assertNotIn("--retire-source-embeddings", text)
        self.assertNotIn("SOURCE_EMBEDDINGS_RETIRED", text)

    def test_workflow_authenticates_all_published_archives(self) -> None:
        text = SCRIPT.read_text(encoding="utf-8")

        self.assertIn("mmfp_embeddings_prott5.tar.gz", text)
        self.assertIn("mmfp_embeddings_struct_ppi.tar.gz", text)
        self.assertIn("mmfp_embeddings_text_temporal.tar.gz", text)
        self.assertIn("Published archive checksum mismatch", text)


if __name__ == "__main__":
    unittest.main()
