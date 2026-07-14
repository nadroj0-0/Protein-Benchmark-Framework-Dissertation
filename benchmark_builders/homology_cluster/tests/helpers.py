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
        "uniref90_fasta": InputSpec("uniref90_fasta", FIXTURES / "uniref90.fasta", release="2026_02"),
        "idmapping": InputSpec("idmapping", FIXTURES / "idmapping_selected.tab", release="2026_02"),
        "uniprot_sequences": InputSpec("uniprot_sequences", FIXTURES / "uniprot.fasta", release="2026_02"),
        "goa": InputSpec("goa", FIXTURES / "goa.gaf", release="234"),
        "go_obo": InputSpec("go_obo", FIXTURES / "go-mini.obo", release="releases/2026-06-15"),
        "cluster_assignments": FIXTURES / "clusters.tsv",
        "fixture_mode": True,
        "split_policy": "sequence-balanced",
        "min_count": 1,
        "threads": 2,
        "strict_qc": True,
    }
    values.update(overrides)
    return BuildConfig(**values)
