from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest

from homology_cluster_benchmark.frozen_inputs import (
    bind_frozen_inputs,
    expected_input_names,
    write_synthetic_fixture_manifest,
    load_frozen_input_manifest,
)
from homology_cluster_benchmark.clustering import (
    connect_proteins_to_clusters,
    mapping_counters_by_source,
    retained_cluster_info,
)
from homology_cluster_benchmark.idmapping import load_uniref90_mappings
from homology_cluster_benchmark.inputs import resolve_input
from homology_cluster_benchmark.mapping import load_requested_proteins_from_sources
from homology_cluster_benchmark.mmseqs import ClusterIndex
from homology_cluster_benchmark.models import InputSpec, MappingDecision, ProteinCatalog
from homology_cluster_benchmark.uniref import UniRefIndex

from tests.helpers import FIXTURES, fixture_config


def _dat(
    path: Path,
    records: list[tuple[str, tuple[str, ...], str]],
    review_status: str = "Reviewed",
) -> Path:
    chunks = []
    for primary, aliases, sequence in records:
        accessions = "; ".join((primary, *aliases)) + ";"
        chunks.append(
            f"ID   {primary}_ENTRY {review_status}; {len(sequence)} AA.\n"
            f"AC   {accessions}\n"
            "OX   NCBI_TaxID=9606;\n"
            f"SQ   SEQUENCE   {len(sequence)} AA;\n"
            f"     {sequence}\n"
            "//\n"
        )
    path.write_text("".join(chunks), encoding="utf-8")
    return path


class UniProtScopeTests(unittest.TestCase):
    def test_each_scope_requires_exactly_its_declared_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = fixture_config(root / "out", root / "temp")
            base.validate()
            trembl = replace(
                base,
                uniprot_source_scope="trembl-only",
                uniprot_sprot_sequences=None,
                uniprot_trembl_sequences=InputSpec(
                    "uniprot_trembl_sequences", FIXTURES / "uniprot.fasta",
                    release="2026_02", source_population="trembl",
                ),
            )
            trembl.validate()
            combined = replace(
                base,
                uniprot_source_scope="sprot-and-trembl",
                uniprot_trembl_sequences=trembl.uniprot_trembl_sequences,
            )
            combined.validate()
            with self.assertRaisesRegex(ValueError, "requires uniprot_trembl_sequences"):
                replace(base, uniprot_source_scope="sprot-and-trembl").validate()
            with self.assertRaisesRegex(ValueError, "forbids irrelevant source"):
                replace(
                    base,
                    uniprot_trembl_sequences=trembl.uniprot_trembl_sequences,
                ).validate()
            with self.assertRaisesRegex(ValueError, "explicitly set"):
                replace(base, uniprot_source_scope="").validate()

    def test_scope_specific_manifest_cardinality_roles_and_determinism(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = fixture_config(root / "out", root / "temp")
            cases = {
                "sprot-only": {"uniprot_sprot_sequences": base.uniprot_sprot_sequences},
                "trembl-only": {"uniprot_trembl_sequences": InputSpec(
                    "uniprot_trembl_sequences", FIXTURES / "uniprot.fasta",
                    release="2026_02", source_population="trembl",
                )},
                "sprot-and-trembl": {
                    "uniprot_sprot_sequences": base.uniprot_sprot_sequences,
                    "uniprot_trembl_sequences": InputSpec(
                        "uniprot_trembl_sequences", FIXTURES / "uniprot.fasta",
                        release="2026_02", source_population="trembl",
                    ),
                },
            }
            for scope, source_specs in cases.items():
                with self.subTest(scope=scope):
                    shared = {
                        name: getattr(base, name)
                        for name in ("uniref90_fasta", "idmapping", "goa", "go_obo")
                    }
                    specs = {**shared, **source_specs}
                    resolved = {
                        name: resolve_input(spec, root / "downloads", allow_downloads=False)
                        for name, spec in specs.items()
                    }
                    first = write_synthetic_fixture_manifest(
                        root / f"{scope}-a.json", specs, resolved, scope
                    )
                    second = write_synthetic_fixture_manifest(
                        root / f"{scope}-b.json", specs, resolved, scope
                    )
                    self.assertEqual(len(first.entries), 5 if scope != "sprot-and-trembl" else 6)
                    self.assertEqual(set(first.entries), set(expected_input_names(scope)))
                    self.assertEqual(first.source_fingerprint, second.source_fingerprint)
                    self.assertEqual(
                        bind_frozen_inputs(first, specs, resolved),
                        bind_frozen_inputs(second, specs, resolved),
                    )
                    payload = json.loads((root / f"{scope}-a.json").read_text())
                    payload["inputs"][0]["logical_role"] = "wrong-role"
                    wrong_role = root / f"{scope}-wrong-role.json"
                    wrong_role.write_text(json.dumps(payload))
                    with self.assertRaisesRegex(ValueError, "logical_role mismatch"):
                        load_frozen_input_manifest(
                            wrong_role, uniprot_source_scope=scope, fixture_mode=True
                        )

            manifest = root / "sprot-only-a.json"
            with self.assertRaisesRegex(ValueError, "source scope mismatch"):
                load_frozen_input_manifest(
                    manifest, uniprot_source_scope="trembl-only", fixture_mode=True
                )
            payload = json.loads(manifest.read_text())
            payload["inputs"][0]["release"] = "different-release"
            release_mismatch = root / "release-mismatch.json"
            release_mismatch.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "releases disagree"):
                load_frozen_input_manifest(
                    release_mismatch,
                    uniprot_source_scope="sprot-only",
                    fixture_mode=True,
                )
            payload = json.loads(manifest.read_text())
            payload["inputs"][0]["sha256"] = "not-a-hash"
            invalid_hash = root / "invalid-hash.json"
            invalid_hash.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                load_frozen_input_manifest(
                    invalid_hash,
                    uniprot_source_scope="sprot-only",
                    fixture_mode=True,
                )

    def test_combined_dat_is_fixed_order_and_reports_primary_secondary_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sprot = _dat(root / "sprot.dat", [("P1", ("S1",), "MAAA")])
            trembl = _dat(
                root / "trembl.dat", [("T1", ("TALT",), "MCCC")], "Unreviewed"
            )
            first = load_requested_proteins_from_sources(
                {"trembl": trembl, "sprot": sprot},
                {"S1", "T1"}, strict_collisions=True,
                collision_database=root / "first.sqlite",
            )
            second = load_requested_proteins_from_sources(
                {"sprot": sprot, "trembl": trembl},
                {"S1", "T1"}, strict_collisions=True,
                collision_database=root / "second.sqlite",
            )
            self.assertEqual(sorted(first.records), sorted(second.records))
            self.assertEqual(first.alias_to_primary["S1"], "P1")
            self.assertEqual(first.primary_source["P1"], "sprot")
            self.assertEqual(first.primary_source["T1"], "trembl")
            self.assertEqual(first.source_counts["sprot"]["primary_accessions_read"], 1)
            self.assertEqual(first.source_counts["trembl"]["secondary_aliases_read"], 1)

    def test_identical_primary_and_conflicting_sequence_collisions_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = _dat(root / "sprot.dat", [("P1", ("SHARED",), "MAAA")])
            cases = {
                "duplicate-primary-identical": [("P1", (), "MAAA")],
                "conflicting-sequence": [("P1", (), "MCCC")],
            }
            for expected, records in cases.items():
                with self.subTest(expected=expected):
                    other = _dat(root / f"{expected}.dat", records, "Unreviewed")
                    with self.assertRaisesRegex(ValueError, expected):
                        load_requested_proteins_from_sources(
                            {"sprot": base, "trembl": other},
                            {"P1", "SHARED", "T1"},
                            strict_collisions=True,
                            collision_database=root / f"{expected}.sqlite",
                        )

    def test_identical_ambiguous_secondary_is_reported_and_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sprot = _dat(root / "sprot.dat", [
                ("P68744", ("P18556",), "MAAA"),
                ("P68745", ("P18556",), "MAAA"),
            ])

            catalog = load_requested_proteins_from_sources(
                {"sprot": sprot},
                {"P18556"},
                strict_collisions=True,
                collision_database=root / "collisions.sqlite",
            )

            self.assertEqual(set(catalog.records), {"P68744", "P68745"})
            self.assertIn("P18556", catalog.ambiguous_aliases)
            self.assertNotIn("P18556", catalog.alias_to_primary)
            self.assertEqual(
                catalog.collision_counts["ambiguous-secondary-identical"], 1
            )
            self.assertEqual(
                catalog.source_counts["sprot"]["ambiguous_secondary_aliases"], 1
            )

    def test_dat_review_status_must_match_declared_source_role(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reviewed = _dat(root / "reviewed.dat", [("P1", (), "MAAA")])
            unreviewed = _dat(
                root / "unreviewed.dat", [("T1", (), "MCCC")], "Unreviewed"
            )
            for source, path in (("sprot", unreviewed), ("trembl", reviewed)):
                with self.subTest(source=source), self.assertRaisesRegex(
                    ValueError, "source-role mismatch"
                ):
                    load_requested_proteins_from_sources(
                        {source: path},
                        {"P1", "T1"},
                        strict_collisions=True,
                        collision_database=root / f"{source}.sqlite",
                    )

            missing_id = root / "missing-id.dat"
            missing_id.write_text(
                "AC   P1;\nSQ   SEQUENCE   4 AA;\n     MAAA\n//\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "has no ID line"):
                load_requested_proteins_from_sources(
                    {"sprot": missing_id},
                    {"P1"},
                    strict_collisions=True,
                    collision_database=root / "missing-id.sqlite",
                )

    def test_source_conflict_counts_are_measured_not_fabricated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sprot = _dat(root / "sprot.dat", [("P1", (), "MAAA")])
            trembl = _dat(
                root / "trembl.dat", [("P1", (), "MCCC")], "Unreviewed"
            )
            catalog = load_requested_proteins_from_sources(
                {"sprot": sprot, "trembl": trembl},
                {"UNRELATED"},
                strict_collisions=False,
                collision_database=root / "collisions.sqlite",
            )
            self.assertEqual(catalog.source_counts["sprot"]["conflicting_sequences"], 1)
            self.assertEqual(catalog.source_counts["trembl"]["conflicting_sequences"], 1)
            counts = mapping_counters_by_source(
                [], {"sprot", "trembl"}, catalog.source_counts
            )
            self.assertEqual(counts["sprot"]["conflicting_sequences"], 1)
            self.assertEqual(counts["trembl"]["conflicting_sequences"], 1)

    def test_out_of_scope_idmapping_hit_cannot_retain_a_uniref_cluster(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            uniref = UniRefIndex.build(
                FIXTURES / "uniref90.fasta", root / "uniref.sqlite"
            )
            clusters = ClusterIndex.build(
                FIXTURES / "clusters.tsv", uniref, root / "clusters.sqlite"
            )
            decisions = load_uniref90_mappings(
                FIXTURES / "idmapping_selected.tab",
                {"P1BP"},
                ProteinCatalog(),
                uniref,
            )
            self.assertEqual(decisions[0].status, "not-in-selected-uniprot")
            self.assertFalse(decisions[0].canonical_sequence_available)
            connected = connect_proteins_to_clusters(decisions, clusters)
            self.assertFalse(connected[0].mmseqs_cluster_id)
            self.assertEqual(retained_cluster_info(connected, clusters), {})

    def test_source_mapping_counts_separate_accession_uniref_and_cluster_stages(self):
        decisions = [
            MappingDecision(
                raw_accession="S1",
                protein_id="P1",
                accession_action="secondary-to-primary",
                uniref90_id="UniRef90_P1",
                status="missing-mmseqs-assignment",
                exists_in_fasta=True,
                canonical_sequence_available=True,
                source_population="sprot",
            )
        ]
        counts = mapping_counters_by_source(
            decisions,
            {"sprot", "trembl"},
            {"sprot": {"conflicting_sequences": 0}, "trembl": {}},
        )
        self.assertEqual(counts["sprot"]["mapped_by_secondary_alias"], 1)
        self.assertEqual(counts["sprot"]["mapped_to_uniref90"], 1)
        self.assertEqual(counts["sprot"]["mapped_to_mmseqs_cluster"], 0)
        self.assertEqual(counts["trembl"]["goa_accessions"], 0)


if __name__ == "__main__":
    unittest.main()
