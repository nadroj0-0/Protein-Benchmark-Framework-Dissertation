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
