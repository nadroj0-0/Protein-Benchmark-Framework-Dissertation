import csv
from dataclasses import replace
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

import numpy as np

from pfp_embedding_inventory.models import (
    ArtifactScopeSpec,
    ArtifactVerification,
    BenchmarkData,
    BenchmarkContract,
    MODALITIES,
    ModalitySpec,
    PlannerConfig,
    ProvenanceSpec,
)
from pfp_embedding_inventory.provenance import compute_cache_catalog


DIMS = {"prott5": 1024, "text": 768, "structure": 512, "ppi": 512}
DIRS = {
    "prott5": "prott5",
    "text": "exp_text_embeddings_temporal",
    "structure": "IF1",
    "ppi": "ppi",
}
ONTOLOGIES = ("bp", "cc", "mf")
SPLITS = ("training", "validation", "test")


def make_contract(
    id_overlap: str = "allow", sequence_overlap: str = "allow"
) -> BenchmarkContract:
    return BenchmarkContract(
        id_overlap=id_overlap,
        sequence_overlap=sequence_overlap,
        protein_id_pattern=r"^[^\s/\\]+$",
        sequence_pattern=r"^[A-Za-z*.-]+$",
    )


def make_config(
    compatibility: Optional[Mapping[str, str]] = None,
    target_contract: Optional[BenchmarkContract] = None,
    source_contract: Optional[BenchmarkContract] = None,
) -> PlannerConfig:
    compatibility = compatibility or {}
    modalities: Dict[str, ModalitySpec] = {}
    for modality in MODALITIES:
        state = compatibility.get(modality, "compatible")
        identity = "test-%s|v1" % modality if modality == "structure" else "test-%s-v1" % modality
        modalities[modality] = ModalitySpec(
            name=modality,
            directory=DIRS[modality],
            expected_dim=DIMS[modality],
            sequence_dependent=modality in {"prott5", "structure"},
            allow_sequence_hash_reuse=modality == "prott5",
            provenance=ProvenanceSpec(
                compatibility=state,
                label="test-%s" % modality,
                source_identity=identity,
                target_identity=identity,
                evidence="synthetic test evidence",
                text_role_policy="none",
                requires_mapping_evidence=modality != "prott5",
            ),
        )
    return PlannerConfig(
        schema_version=3,
        name="synthetic",
        target_benchmark_contract=target_contract or make_contract(),
        source_benchmark_contract=source_contract or make_contract(),
        modalities=modalities,
        artifact_scope=ArtifactScopeSpec(
            mode="none",
            artifact_id="",
            metadata_url="",
            expected_benchmark_fingerprint="",
            expected_cache_catalog_fingerprint="",
            expected_modality_counts={},
            expected_total_files=0,
            expected_total_bytes=0,
            archives=(),
            expected_reference_commit="",
            reference_files=(),
        ),
    )


def bind_verified_cache(
    config: PlannerConfig, source: BenchmarkData, cache: Path
) -> Tuple[PlannerConfig, ArtifactVerification]:
    catalog = compute_cache_catalog(cache, config)
    artifact_id = "synthetic-verified-cache"
    scope = replace(
        config.artifact_scope,
        mode="verified-published-cache",
        artifact_id=artifact_id,
        expected_benchmark_fingerprint=source.fingerprint,
        expected_cache_catalog_fingerprint=catalog.fingerprint,
        expected_modality_counts=catalog.modality_counts,
        expected_total_files=catalog.total_files,
        expected_total_bytes=catalog.total_bytes,
    )
    bound_config = replace(config, artifact_scope=scope)
    proof = ArtifactVerification(
        configured=True,
        verified=True,
        artifact_id=artifact_id,
        checks={"synthetic-cache-proof": True},
        reasons=[],
        expected={},
        observed={
            "source_benchmark_fingerprint": source.fingerprint,
            "embedding_cache_root": str(cache.resolve()),
            "cache_catalog": catalog.as_dict(),
        },
    )
    return bound_config, proof


def write_nine_csvs(
    directory: Path,
    rows_by_file: Optional[
        Mapping[Tuple[str, str], Iterable[Tuple[str, str, str]]]
    ] = None,
    shared_rows: Optional[Iterable[Tuple[str, str, str]]] = None,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    if rows_by_file is None:
        shared = list(shared_rows or [("P1", "ACDE", "1")])
        rows_by_file = {
            (ontology, split): shared for ontology in ONTOLOGIES for split in SPLITS
        }
    for ontology in ONTOLOGIES:
        for split in SPLITS:
            path = directory / ("%s-%s.csv" % (ontology, split))
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["proteins", "sequences", "GO:0000001"])
                writer.writerows(rows_by_file[(ontology, split)])


def unique_rows_by_file() -> Dict[Tuple[str, str], List[Tuple[str, str, str]]]:
    rows: Dict[Tuple[str, str], List[Tuple[str, str, str]]] = {}
    counter = 0
    for ontology in ONTOLOGIES:
        for split in SPLITS:
            counter += 1
            rows[(ontology, split)] = [("P%d" % counter, "ACD%s" % ("E" * counter), "1")]
    return rows


def make_cache(cache: Path, protein_ids: Iterable[str], modalities: Iterable[str] = MODALITIES) -> None:
    for modality in modalities:
        directory = cache / DIRS[modality]
        directory.mkdir(parents=True, exist_ok=True)
        for protein_id in protein_ids:
            np.save(directory / (protein_id + ".npy"), np.arange(DIMS[modality], dtype=np.float32))
