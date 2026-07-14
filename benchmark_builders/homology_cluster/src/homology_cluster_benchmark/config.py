from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .models import InputSpec


SUPPORTED_IDENTITIES = (0.30, 0.25, 0.20, 0.15, 0.10, 0.05)
FROZEN_UNIPROT_RELEASE = "2026_02"
FROZEN_GOA_RELEASE = "234"
FROZEN_ONTOLOGY_RELEASE = "releases/2026-06-15"
IDENTITY_DIRECTORY = {
    0.30: "identity_30",
    0.25: "identity_25",
    0.20: "identity_20",
    0.15: "identity_15",
    0.10: "identity_10",
    0.05: "identity_05",
}

SUPERVISOR_EVIDENCE_CODES = frozenset({
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP",
    "HGI", "HEP", "TAS", "NAS", "IGC", "RCA", "ND", "IC",
})

ROOT_TERMS = frozenset({"GO:0008150", "GO:0005575", "GO:0003674"})
ASPECT_TO_NAMESPACE = {
    "P": "biological_process",
    "C": "cellular_component",
    "F": "molecular_function",
}
PREFIX_TO_NAMESPACE = {
    "bp": "biological_process",
    "cc": "cellular_component",
    "mf": "molecular_function",
}
SPLITS = ("training", "validation", "test")
SPLIT_POLICIES = ("cluster-count-random", "sequence-balanced")
TRAINING_POPULATIONS = ("annotated-only", "all-cluster-members")


def parse_identity(value: str | int | float) -> float:
    if isinstance(value, str):
        stripped = value.strip().lower()
        if stripped == "all":
            raise ValueError("'all' represents multiple identities, not one identity")
        numeric = float(stripped)
    else:
        numeric = float(value)
    if numeric > 1:
        numeric /= 100.0
    for supported in SUPPORTED_IDENTITIES:
        if abs(numeric - supported) < 1e-12:
            return supported
    allowed = ", ".join(str(int(item * 100)) for item in SUPPORTED_IDENTITIES)
    raise ValueError(f"Unsupported identity {value!r}; Daniel's allowed percentages are {allowed}")


@dataclass(frozen=True)
class BuildConfig:
    identity: float
    output_dir: Path
    temp_dir: Path
    uniref90_fasta: InputSpec
    idmapping: InputSpec
    uniprot_sequences: InputSpec
    goa: InputSpec
    go_obo: InputSpec
    split_policy: str = "sequence-balanced"
    training_population: str = "annotated-only"
    mmseqs_bin: str = "mmseqs"
    expected_mmseqs_version: str | None = None
    cluster_assignments: Path | None = None
    frozen_input_manifest: Path | None = None
    fixture_mode: bool = False
    threads: int = 1
    seed: int = 0
    min_count: int = 50
    development_fraction: float = 0.80
    training_fraction_within_development: float = 0.90
    coverage: float = 0.80
    cov_mode: int = 0
    cluster_mode: int = 0
    alignment_mode: int = 3
    cluster_reassign: int = 1
    sensitivity: float = 7.5
    evalue: float = 1e-4
    include_relationships: bool = True
    evidence_codes: frozenset[str] = field(default_factory=lambda: SUPERVISOR_EVIDENCE_CODES)
    allow_downloads: bool = True
    strict_qc: bool = True
    allow_empty_fixture_outputs: bool = False
    keep_temp: bool = False
    release_uniprot: str = "2026_02"
    release_goa: str = "234"
    release_ontology: str = "releases/2026-06-15"
    giant_cluster_threshold: int = 10000
    scratch_safety_multiplier: float = 8.0
    minimum_free_disk_bytes: int = 0
    persistent_results_root: Path | None = None
    mmseqs_work_multiplier: float = 8.0
    publication_safety_multiplier: float = 2.0
    excluded_sample_per_reason: int = 1000

    def validate(self, require_pinned_inputs: bool = True) -> None:
        canonical_identity = parse_identity(self.identity)
        if canonical_identity != self.identity:
            raise ValueError("identity must be stored as a canonical fraction")
        if self.split_policy not in SPLIT_POLICIES:
            raise ValueError(f"Unsupported split policy: {self.split_policy}")
        if self.training_population not in TRAINING_POPULATIONS:
            raise ValueError(f"Unsupported training population: {self.training_population}")
        if self.training_population == "all-cluster-members":
            raise ValueError(
                "all-cluster-members is intentionally unsupported: it cannot emit supervised PFP "
                "rows without an approved label and embedding-cost policy; unannotated members "
                "remain visible only in retained-member manifests"
            )
        if self.coverage != 0.80:
            raise ValueError("Coverage is methodologically locked to exactly 0.80")
        if self.cov_mode != 0:
            raise ValueError("The frozen UniRef longest-sequence overlap policy requires --cov-mode 0")
        if self.cluster_mode != 0 or self.alignment_mode != 3 or self.cluster_reassign != 1:
            raise ValueError("MMseqs2 clustering/alignment/reassignment modes are fixed by this implementation")
        if self.development_fraction != 0.80:
            raise ValueError("Daniel's development/test policy is locked to exactly 0.80/0.20")
        if self.training_fraction_within_development != 0.90:
            raise ValueError(
                "The established development split is locked to exactly 0.90/0.10 "
                "training/validation"
            )
        if self.sensitivity != 7.5:
            raise ValueError("MMseqs2 sensitivity is methodologically locked to exactly 7.5")
        if self.evalue != 1e-4:
            raise ValueError("MMseqs2 E-value is methodologically locked to exactly 1e-4")
        if self.threads < 1:
            raise ValueError("threads must be positive")
        if self.min_count < 1:
            raise ValueError("min_count must be positive")
        if not self.fixture_mode and self.min_count < 50:
            raise ValueError("Production min_count is locked to at least 50; lower values require fixture mode")
        if self.evidence_codes != SUPERVISOR_EVIDENCE_CODES:
            raise ValueError("The evidence-code policy is locked to Daniel's exact supplied set")
        if not self.include_relationships and not self.fixture_mode:
            raise ValueError(
                "Production GO propagation is locked to the established builder's relationship policy"
            )
        if not self.strict_qc and not self.fixture_mode:
            raise ValueError("Production builds cannot disable strict contract validation")
        if self.cluster_assignments is not None and not self.fixture_mode:
            raise ValueError(
                "Precomputed --cluster-assignments are fixture-only because their generating "
                "MMseqs2 identity, coverage, command, and version cannot be proven by this run"
            )
        if self.allow_empty_fixture_outputs and not self.fixture_mode:
            raise ValueError("Empty ontology/split outputs may be allowed only in fixture mode")
        if require_pinned_inputs and not self.fixture_mode:
            if self.frozen_input_manifest is None:
                raise ValueError(
                    "Production requires --frozen-input-manifest; caller-supplied self-hashes alone "
                    "do not establish recorded frozen-source provenance"
                )
            expected_version = (self.expected_mmseqs_version or "").strip()
            placeholders = {"", "unknown", "latest", "replace_me", "placeholder"}
            lowered_version = expected_version.lower()
            if (
                lowered_version in placeholders
                or any(token in lowered_version for token in ("replace", "placeholder", "unknown"))
                or "<" in expected_version
                or ">" in expected_version
            ):
                raise ValueError(
                    "Production requires an exact non-placeholder --expected-mmseqs-version"
                )
        frozen = {
            "UniProt/UniRef": (self.release_uniprot, FROZEN_UNIPROT_RELEASE),
            "GOA": (self.release_goa, FROZEN_GOA_RELEASE),
            "ontology": (self.release_ontology, FROZEN_ONTOLOGY_RELEASE),
        }
        for label, (observed, expected) in frozen.items():
            if observed != expected:
                raise ValueError(
                    f"{label} release is frozen to {expected!r}; observed configured label {observed!r}"
                )
        expected_spec_releases = {
            "uniref90_fasta": self.release_uniprot,
            "idmapping": self.release_uniprot,
            "uniprot_sequences": self.release_uniprot,
            "goa": self.release_goa,
            "go_obo": self.release_ontology,
        }
        for name, expected_release in expected_spec_releases.items():
            spec = getattr(self, name)
            if spec.release != expected_release:
                raise ValueError(
                    f"{name} release metadata {spec.release!r} does not match {expected_release!r}"
                )
            if require_pinned_inputs and not self.fixture_mode and not spec.expected_sha256:
                raise ValueError(
                    f"Production input {name} requires an expected SHA-256 so a file cannot merely "
                    "be relabelled as the frozen release"
                )
        if self.scratch_safety_multiplier < 1:
            raise ValueError("scratch_safety_multiplier must be at least 1")
        if self.mmseqs_work_multiplier < 1:
            raise ValueError("mmseqs_work_multiplier must be at least 1")
        if self.publication_safety_multiplier < 1:
            raise ValueError("publication_safety_multiplier must be at least 1")
        if self.minimum_free_disk_bytes < 0:
            raise ValueError("minimum_free_disk_bytes cannot be negative")
        if self.excluded_sample_per_reason < 0:
            raise ValueError("excluded_sample_per_reason cannot be negative")

    @property
    def identity_directory(self) -> str:
        return IDENTITY_DIRECTORY[self.identity]

    @property
    def publication_relative_path(self) -> Path:
        return (
            Path(self.identity_directory)
            / self.split_policy
            / self.training_population
            / f"seed_{self.seed}"
            / f"min_count_{self.min_count}"
        )
