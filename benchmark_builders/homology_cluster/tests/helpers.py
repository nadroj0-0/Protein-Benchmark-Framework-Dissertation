from __future__ import annotations

from pathlib import Path

from homology_cluster_benchmark.config import BuildConfig
from homology_cluster_benchmark.models import InputSpec


FIXTURES = Path(__file__).resolve().parent / "fixtures"


def fixture_config(output_root: Path, temp_root: Path, **overrides) -> BuildConfig:
    values = {
        "identity": 0.30,
        "output_dir": output_root,
        "temp_dir": temp_root,
        "uniref90_fasta": InputSpec(
            "uniref90_fasta", FIXTURES / "uniref90.fasta", release="2026_02",
            source_population="uniref90-clustering-scaffold",
        ),
        "idmapping": InputSpec(
            "idmapping", FIXTURES / "idmapping_selected.tab", release="2026_02",
            source_population="uniprotkb-shared-mapping",
        ),
        "uniprot_source_scope": "sprot-only",
        "uniprot_sprot_sequences": InputSpec(
            "uniprot_sprot_sequences", FIXTURES / "uniprot.fasta", release="2026_02",
            source_population="sprot",
        ),
        "uniprot_trembl_sequences": None,
        "goa": InputSpec(
            "goa", FIXTURES / "goa.gaf", release="234",
            source_population="uniprotkb-goa",
        ),
        "go_obo": InputSpec(
            "go_obo", FIXTURES / "go-mini.obo", release="releases/2026-06-15",
            source_population="gene-ontology",
        ),
        "cluster_assignments": FIXTURES / "clusters.tsv",
        "fixture_mode": True,
        "split_policy": "sequence-balanced",
        "min_count": 1,
        "threads": 2,
        "strict_qc": True,
    }
    values.update(overrides)
    return BuildConfig(**values)
