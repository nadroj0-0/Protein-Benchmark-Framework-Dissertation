from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


# Final CAFA3 benchmark policy from benchmark20171115/00README.txt.
CAFA3_FINAL_EXP_CODES = frozenset({
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC",
})

# Public CAFA_benchmark/create_benchmark.py policy. The final CAFA3 release
# extended this six-code set with TAS and IC.
CAFA3_PUBLIC_PYTHON_EXP_CODES = frozenset({
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
})

# Dissertation-supervisor policy. It is deliberately a named alternative, not
# folded into the CAFA3 profile, so results can state exactly which was used.
SUPERVISOR_EXP_CODES = frozenset({
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP",
    "HGI", "HEP", "TAS", "NAS", "IGC", "RCA", "ND", "IC",
})

EVIDENCE_POLICIES = {
    "cafa3-final": CAFA3_FINAL_EXP_CODES,
    "cafa3-public-python": CAFA3_PUBLIC_PYTHON_EXP_CODES,
    "supervisor": SUPERVISOR_EXP_CODES,
}


@dataclass(frozen=True)
class BenchmarkProfile:
    name: str
    training_taxon_policy: str
    target_taxon_policy: str
    evidence_policy: str
    training_reviewed_only: bool
    target_reviewed_only: bool
    t0_cutoff: str
    exclude_t1_backfill: bool = True
    require_t0_presence: bool = True
    sequence_change_policy: str = "exclude"
    protein_binding_policy: str = "drop-mf-protein-binding-only"


# Profiles encode policy only. Paths remain required CLI arguments.
BENCHMARK_PROFILES = {
    "cafa3-reconstructed": BenchmarkProfile(
        name="cafa3-reconstructed",
        training_taxon_policy="all",
        target_taxon_policy="cafa3-targets",
        evidence_policy="cafa3-final",
        training_reviewed_only=True,
        target_reviewed_only=False,
        t0_cutoff="20170213",
    ),
    "contemporary-cafa3-style": BenchmarkProfile(
        name="contemporary-cafa3-style",
        training_taxon_policy="all",
        target_taxon_policy="cafa3-targets",
        evidence_policy="cafa3-final",
        training_reviewed_only=True,
        target_reviewed_only=False,
        t0_cutoff="20250308",
    ),
    "supervisor": BenchmarkProfile(
        name="supervisor",
        training_taxon_policy="cafa3-targets",
        target_taxon_policy="cafa3-targets",
        evidence_policy="supervisor",
        training_reviewed_only=False,
        target_reviewed_only=False,
        t0_cutoff="20250308",
    ),
}

ASPECT_TO_PREFIX = {"P": "bp", "C": "cc", "F": "mf"}

PREFIX_TO_NAMESPACE = {
    "bp": "biological_process",
    "cc": "cellular_component",
    "mf": "molecular_function",
}

NAMESPACE_TO_PREFIX = {v: k for k, v in PREFIX_TO_NAMESPACE.items()}


@dataclass(frozen=True)
class BuildConfig:
    uniprot_t0: tuple[Path, ...]
    uniprot_t1: tuple[Path, ...]
    goa_t0: Path
    goa_t1: Path
    go_obo: Path
    output_dir: Path
    go_obo_t0: Path | None = None
    go_obo_t1: Path | None = None
    report_dir: Path | None = None
    profile_name: str = "contemporary-cafa3-style"
    training_taxa: frozenset[str] = field(default_factory=frozenset)
    target_taxa: frozenset[str] = field(default_factory=frozenset)
    evidence_codes: frozenset[str] = CAFA3_FINAL_EXP_CODES
    t0_cutoff: str | None = "20250308"
    exclude_t1_backfill: bool = True
    require_t0_presence: bool = True
    sequence_change_policy: str = "exclude"
    protein_binding_policy: str = "drop-mf-protein-binding-only"
    min_count: int = 50
    split: float = 0.9
    seed: int = 0
    reviewed_only: bool = True
    target_reviewed_only: bool = False
    include_rels: bool = True
    write_intermediates: bool = True
    write_checksums: bool = True
    strict_qc: bool = True
    max_gaf_records: int | None = None

    @property
    def ontology_t0(self) -> Path:
        return self.go_obo_t0 or self.go_obo

    @property
    def ontology_t1(self) -> Path:
        return self.go_obo_t1 or self.go_obo

    @property
    def reports(self) -> Path:
        return self.report_dir or self.output_dir / "reports"


def normalise_taxa(values: Iterable[str]) -> frozenset[str]:
    taxa = set()
    for value in values:
        value = value.strip()
        if not value:
            continue
        if value.startswith("taxon:"):
            value = value.split(":", 1)[1]
        taxa.add(value)
    return frozenset(taxa)


def normalise_gaf_date(value: str | None) -> str | None:
    if value is None:
        return None
    compact = value.replace("-", "").strip()
    if len(compact) != 8 or not compact.isdigit():
        raise ValueError(f"Expected a GAF date in YYYYMMDD or YYYY-MM-DD form, got {value!r}")
    return compact
