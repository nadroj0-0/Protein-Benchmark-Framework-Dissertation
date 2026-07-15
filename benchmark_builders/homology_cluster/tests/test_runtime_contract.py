from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

from homology_cluster_benchmark.attrition import load_attrition_policy
from homology_cluster_benchmark.frozen_inputs import load_frozen_input_manifest
from homology_cluster_benchmark.runtime_contract import (
    RuntimeInput,
    main,
    write_runtime_contract,
)


FIXTURES = Path(__file__).resolve().parent / "fixtures"
COMMIT = "a" * 40
URLS = {
    "uniref90_fasta": (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "uniref/uniref90/uniref90.fasta.gz"
    ),
    "idmapping": (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "knowledgebase/idmapping/idmapping_selected.tab.gz"
    ),
    "uniprot_sprot_sequences": (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "knowledgebase/complete/uniprot_sprot.dat.gz"
    ),
    "goa": "https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz",
    "go_obo": "https://release.geneontology.org/2026-06-19/ontology/go-basic.obo",
}


def _runtime_files(root: Path) -> dict[str, Path]:
    paths = {
        "uniref90_fasta": root / "uniref90.fasta.gz",
        "idmapping": root / "idmapping_selected.tab.gz",
        "uniprot_sprot_sequences": root / "uniprot_sprot.dat.gz",
        "goa": root / "goa_uniprot_all.gaf.234.gz",
        "go_obo": root / "go-basic.obo",
    }
    for name in ("uniref90_fasta", "idmapping", "uniprot_sprot_sequences"):
        paths[name].write_bytes(f"fixture-{name}\n".encode())
    with gzip.open(paths["goa"], "wb") as handle:
        handle.write((FIXTURES / "goa.gaf").read_bytes())
    paths["go_obo"].write_bytes((FIXTURES / "go-mini.obo").read_bytes())
    return paths


class RuntimeContractTests(unittest.TestCase):
    def test_runtime_contract_is_accepted_by_existing_strict_loaders(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = _runtime_files(root)
            manifest_path = root / "manifest.json"
            policy_path = root / "policy.json"
            result = write_runtime_contract(
                manifest_path,
                policy_path,
                "sprot-only",
                COMMIT,
                [
                    RuntimeInput(
                        name,
                        paths[name],
                        URLS[name],
                        "provided-path-staged-to-scratch",
                    )
                    for name in paths
                ],
            )
            manifest = load_frozen_input_manifest(
                manifest_path, uniprot_source_scope="sprot-only"
            )
            self.assertEqual(manifest.sha256, result["manifest_sha256"])
            policy, digest = load_attrition_policy(
                policy_path,
                source_scope="sprot-only",
                expected_releases={
                    "uniprot_uniref": "2026_02",
                    "goa": "234",
                    "ontology": "releases/2026-06-15",
                },
                framework_commit=COMMIT,
                frozen_input_manifest_sha256=manifest.sha256,
            )
            self.assertEqual(digest, result["policy_sha256"])
            self.assertIn("without making a prior pilot mandatory", policy["rationale"])

    def test_prepare_cli_does_not_require_the_unselected_trembl_role(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = _runtime_files(root)
            arguments = [
                "prepare",
                "--manifest-out", str(root / "manifest.json"),
                "--policy-out", str(root / "policy.json"),
                "--source-scope", "sprot-only",
                "--framework-revision", COMMIT,
            ]
            for name, path in paths.items():
                option = name.replace("_", "-")
                arguments.extend([
                    f"--{option}", str(path),
                    f"--{option}-url", URLS[name],
                    f"--{option}-acquisition", "downloaded-to-scratch",
                ])
            self.assertEqual(main(arguments), 0)
            self.assertTrue((root / "manifest.json").is_file())

    def test_wrong_embedded_goa_release_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = _runtime_files(root)
            wrong = (FIXTURES / "goa.gaf").read_text().replace(
                "!date-generated: 2026-06-17", "!date-generated: 2026-06-18"
            )
            with gzip.open(paths["goa"], "wt") as handle:
                handle.write(wrong)
            with self.assertRaisesRegex(ValueError, "not release 234"):
                write_runtime_contract(
                    root / "manifest.json",
                    root / "policy.json",
                    "sprot-only",
                    COMMIT,
                    [
                        RuntimeInput(
                            name,
                            paths[name],
                            URLS[name],
                            "downloaded-to-scratch",
                        )
                        for name in paths
                    ],
                )


if __name__ == "__main__":
    unittest.main()
