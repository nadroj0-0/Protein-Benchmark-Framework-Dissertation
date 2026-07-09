from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


# Final CAFA3 benchmark README policy:
# benchmark20171115/00README.txt says EXP, IDA, IPI, IMP, IGI, IEP, TAS, IC.
CAFA3_FINAL_EXP_CODES = frozenset({
    "EXP",
    "IDA",
    "IPI",
    "IMP",
    "IGI",
    "IEP",
    "TAS",
    "IC",
})

# Public CAFA_benchmark/create_benchmark.py policy. This is useful for auditing
# the public Python script, but the final CAFA3 README extends it with TAS and IC.
CAFA3_PUBLIC_PYTHON_EXP_CODES = frozenset({
    "EXP",
    "IDA",
    "IPI",
    "IMP",
    "IGI",
    "IEP",
})

# Supervisor-specified policy for dissertation benchmark variants. Keep this as
# a named preset so the benchmark can switch policy without code changes.
SUPERVISOR_EXP_CODES = frozenset({
    "EXP",
    "IDA",
    "IPI",
    "IMP",
    "IGI",
    "IEP",
    "HTP",
    "HDA",
    "HMP",
    "HGI",
    "HEP",
    "TAS",
    "NAS",
    "IGC",
    "RCA",
    "ND",
    "IC",
})

EVIDENCE_POLICIES = {
    "cafa3-final": CAFA3_FINAL_EXP_CODES,
    "cafa3-public-python": CAFA3_PUBLIC_PYTHON_EXP_CODES,
    "supervisor": SUPERVISOR_EXP_CODES,
}

ASPECT_TO_PREFIX = {
    "P": "bp",
    "C": "cc",
    "F": "mf",
}

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
    target_taxa: frozenset[str] = field(default_factory=frozenset)
    evidence_codes: frozenset[str] = CAFA3_FINAL_EXP_CODES
    min_count: int = 50
    split: float = 0.9
    seed: int = 0
    reviewed_only: bool = False
    include_rels: bool = True
    write_intermediates: bool = True
    max_gaf_records: int | None = None


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
