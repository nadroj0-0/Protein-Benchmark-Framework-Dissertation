import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUBMISSION_SCRIPT = ROOT / "hpc_jobs" / "submit_multilinkage30_campaign.sh"


class Multilinkage30SubmissionTests(unittest.TestCase):
    def test_calibration_receives_campaign_source_run_and_ontology(self) -> None:
        source = SUBMISSION_SCRIPT.read_text(encoding="utf-8")
        match = re.search(
            r"submit_job calibration ml30_cal.*?CALIBRATION_JOB=",
            source,
            flags=re.DOTALL,
        )
        self.assertIsNotNone(match)
        block = match.group(0)
        self.assertIn('--source-run "$MODEL_RUN_DIR"', block)
        self.assertIn('--obo "$HOMOLOGY_OBO"', block)


if __name__ == "__main__":
    unittest.main()
