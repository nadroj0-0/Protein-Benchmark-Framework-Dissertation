from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Iterable

from .config import CAFA3_FINAL_EXP_CODES, normalise_gaf_date
from .goa import load_normalized_annotation_map
from .models import AnnotationLoadResult, ProteinCatalog
from .official_targets import load_official_target_catalog
from .ontology import Ontology
from .parsers import load_protein_catalog


ASPECTS = {
    "bpo": ("BPO", "biological_process"),
    "cco": ("CCO", "cellular_component"),
    "mfo": ("MFO", "molecular_function"),
}
PROTEIN_BINDING = "GO:0005515"
SCHEMA_VERSION = 1


@dataclass(frozen=True)
class OrganizerReplayConfig:
    uniprot_t0: tuple[Path, ...]
    uniprot_t1: tuple[Path, ...]
    goa_t0: Path
    goa_t1: Path
    go_obo: Path
    official_target_fastas: tuple[Path, ...]
    official_target_mapping_dir: Path
    output_dir: Path
    go_obo_t0: Path | None = None
    go_obo_t1: Path | None = None
    official_benchmark_dir: Path | None = None
    evidence_codes: frozenset[str] = CAFA3_FINAL_EXP_CODES
    t0_cutoff: str = "20170213"
    t1_cutoff: str = "20171115"
    t1_endpoint_policy: str = "assigned-date-proxy"
    relationship_policy: str = "is-a-part-of"
    strict_qc: bool = False
    allow_frozen_source_fallback: bool = True
    max_gaf_records: int | None = None
    write_input_checksums: bool = True

    @property
    def ontology_t0(self) -> Path:
        return self.go_obo_t0 or self.go_obo

    @property
    def ontology_t1(self) -> Path:
        return self.go_obo_t1 or self.go_obo


def _relationship_types(policy: str) -> tuple[bool, frozenset[str] | None]:
    if policy == "is-a":
        return False, frozenset()
    if policy == "is-a-part-of":
        return True, frozenset({"part_of"})
    if policy == "all":
        return True, None
    raise ValueError(f"Unknown relationship policy: {policy}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _write_lines(path: Path, values: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{value}\n" for value in values))


def _write_tsv(path: Path, fields: tuple[str, ...], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(field, "")) for field in fields) + "\n")


def _validate(config: OrganizerReplayConfig) -> None:
    if config.output_dir.exists():
        raise FileExistsError(f"Replay output already exists: {config.output_dir}")
    if config.t1_endpoint_policy not in {"assigned-date-proxy", "snapshot-membership"}:
        raise ValueError("t1 endpoint policy must be assigned-date-proxy or snapshot-membership")
    _relationship_types(config.relationship_policy)
    if config.t1_cutoff <= config.t0_cutoff:
        raise ValueError("t1 cutoff must be later than t0 cutoff")
    if config.evidence_codes != CAFA3_FINAL_EXP_CODES:
        raise ValueError(
            "The organiser replay requires the final CAFA3 eight-code evidence policy"
        )
    paths = (
        list(config.uniprot_t0)
        + list(config.uniprot_t1)
        + list(config.official_target_fastas)
        + [config.goa_t0, config.goa_t1, config.go_obo, config.ontology_t0, config.ontology_t1]
    )
    if not config.uniprot_t0 or not config.uniprot_t1:
        raise ValueError("Both t0 and t1 UniProt inputs are required")
    if not config.official_target_fastas:
        raise ValueError("At least one official CAFA3 target FASTA is required")
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    if not config.official_target_mapping_dir.is_dir():
        raise FileNotFoundError(config.official_target_mapping_dir)
    if config.official_benchmark_dir is not None:
        _resolve_official_root(config.official_benchmark_dir)


def _resolve_official_root(path: Path) -> Path:
    candidates = (path, path / "benchmark20171115")
    for candidate in candidates:
        if (candidate / "lists").is_dir() and (candidate / "groundtruth").is_dir():
            return candidate
    raise FileNotFoundError(
        f"Official benchmark root must contain lists/ and groundtruth/: {path}"
    )


def _load_targets(
    config: OrganizerReplayConfig,
) -> tuple[ProteinCatalog, ProteinCatalog, list[dict[str, object]]]:
    t0_reference = load_protein_catalog(config.uniprot_t0, frozenset(), False)
    t1_reference = load_protein_catalog(config.uniprot_t1, frozenset(), False)
    t0 = load_official_target_catalog(
        config.official_target_fastas,
        config.official_target_mapping_dir,
        t0_reference,
        frozenset(),
        "t0",
    )
    t1 = load_official_target_catalog(
        config.official_target_fastas,
        config.official_target_mapping_dir,
        t1_reference,
        frozenset(),
        "t1",
    )
    if set(t0.catalog.records) != set(t1.catalog.records):
        raise ValueError("The fixed CAFA3 target universe changed between mapping passes")
    return t0.catalog, t1.catalog, t0.rows + t1.rows


def _load_annotations(
    config: OrganizerReplayConfig,
    benchmark_go: Ontology,
    t0_go: Ontology,
    t1_go: Ontology,
    t0_targets: ProteinCatalog,
    t1_targets: ProteinCatalog,
) -> tuple[AnnotationLoadResult, AnnotationLoadResult]:
    t0 = load_normalized_annotation_map(
        config.goa_t0,
        alias_to_primary=t0_targets.alias_to_primary,
        source_ontology=t0_go,
        benchmark_ontology=benchmark_go,
        other_ontology=t1_go,
        snapshot="t0",
        allow_frozen_source_fallback=config.allow_frozen_source_fallback,
        evidence_codes=config.evidence_codes,
        max_records=config.max_gaf_records,
    )
    t1 = load_normalized_annotation_map(
        config.goa_t1,
        alias_to_primary=t1_targets.alias_to_primary,
        source_ontology=t1_go,
        benchmark_ontology=benchmark_go,
        other_ontology=t0_go,
        snapshot="t1",
        allow_frozen_source_fallback=config.allow_frozen_source_fallback,
        evidence_codes=config.evidence_codes,
        exclude_on_or_before=config.t0_cutoff,
        include_on_or_before=(
            config.t1_cutoff
            if config.t1_endpoint_policy == "assigned-date-proxy"
            else None
        ),
        max_records=config.max_gaf_records,
    )
    if config.strict_qc and (t0.unmapped_terms or t1.unmapped_terms):
        raise ValueError("Strict QC found GO IDs unresolved in their source ontology")
    return t0, t1


def _direct_by_aspect(
    annotations: dict[str, set[str]], go: Ontology
) -> dict[str, dict[str, set[str]]]:
    by_aspect: dict[str, dict[str, set[str]]] = {
        prefix: defaultdict(set) for prefix in ASPECTS
    }
    namespace_to_prefix = {namespace: prefix for prefix, (_, namespace) in ASPECTS.items()}
    for protein_id, terms in annotations.items():
        for term in terms:
            prefix = namespace_to_prefix[go.get_namespace(term)]
            by_aspect[prefix][protein_id].add(term)
    return {prefix: dict(values) for prefix, values in by_aspect.items()}


def _propagate(go: Ontology, terms: Iterable[str]) -> set[str]:
    propagated: set[str] = set()
    for term in terms:
        propagated.update(go.get_ancestors(term))
    return propagated


def _remove_protein_binding_only(
    annotations: dict[str, set[str]], go: Ontology
) -> tuple[dict[str, set[str]], set[str]]:
    protein_binding_closure = go.get_ancestors(PROTEIN_BINDING)
    retained: dict[str, set[str]] = {}
    removed: set[str] = set()
    for protein_id, direct_terms in annotations.items():
        propagated = _propagate(go, direct_terms)
        if propagated and not (propagated - protein_binding_closure):
            removed.add(protein_id)
        else:
            retained[protein_id] = set(direct_terms)
    return retained, removed


def _leaf_terms(go: Ontology, propagated: set[str]) -> set[str]:
    return {
        term
        for term in propagated
        if not (set(go.get_term(term).get("children", set())) & propagated)
    }


def _build_membership(
    go: Ontology,
    t0_annotations: dict[str, set[str]],
    t1_annotations: dict[str, set[str]],
) -> tuple[
    dict[str, dict[str, set[str]]],
    dict[str, dict[str, set[str]]],
    dict[str, dict[str, set[str]]],
    dict[str, dict[str, set[str]]],
    dict[str, dict[str, set[str]]],
    list[dict[str, object]],
]:
    t0_direct_by_aspect = _direct_by_aspect(t0_annotations, go)
    t1_direct_by_aspect = _direct_by_aspect(t1_annotations, go)
    t0_by_aspect = {
        prefix: {protein: set(terms) for protein, terms in values.items()}
        for prefix, values in t0_direct_by_aspect.items()
    }
    t1_by_aspect = {
        prefix: {protein: set(terms) for protein, terms in values.items()}
        for prefix, values in t1_direct_by_aspect.items()
    }
    t0_by_aspect["mfo"], removed_t0 = _remove_protein_binding_only(
        t0_by_aspect["mfo"], go
    )
    t1_by_aspect["mfo"], removed_t1 = _remove_protein_binding_only(
        t1_by_aspect["mfo"], go
    )

    t0_all = set().union(*(set(values) for values in t0_by_aspect.values()))
    memberships: dict[str, dict[str, set[str]]] = {}
    rows: list[dict[str, object]] = []

    all_proteins = set(t0_annotations) | set(t1_annotations)
    for prefix, (_, namespace) in ASPECTS.items():
        t0_aspect = set(t0_by_aspect[prefix])
        t1_aspect = set(t1_by_aspect[prefix])
        type1 = t1_aspect - t0_all
        type2 = (t1_aspect - t0_aspect) - type1
        typex = type1 | type2
        memberships[prefix] = {"type1": type1, "type2": type2, "typex": typex}
        for protein_id in sorted(all_proteins):
            category = (
                "type1" if protein_id in type1 else
                "type2" if protein_id in type2 else
                "excluded"
            )
            rows.append({
                "protein": protein_id,
                "aspect": prefix,
                "namespace": namespace,
                "t0_any": int(protein_id in t0_all),
                "t0_aspect": int(protein_id in t0_aspect),
                "t1_gained_aspect": int(protein_id in t1_aspect),
                "protein_binding_only_removed_t0": int(
                    prefix == "mfo" and protein_id in removed_t0
                ),
                "protein_binding_only_removed_t1": int(
                    prefix == "mfo" and protein_id in removed_t1
                ),
                "category": category,
            })
    return (
        t0_direct_by_aspect,
        t1_direct_by_aspect,
        t0_by_aspect,
        t1_by_aspect,
        memberships,
        rows,
    )


def _write_annotation_stage(
    path: Path, annotations: dict[str, dict[str, set[str]]]
) -> None:
    rows = []
    for prefix in ASPECTS:
        for protein_id, terms in sorted(annotations[prefix].items()):
            rows.extend(
                {"protein": protein_id, "aspect": prefix, "go_id": term}
                for term in sorted(terms)
            )
    _write_tsv(path, ("protein", "aspect", "go_id"), rows)


def _write_release(
    root: Path,
    go: Ontology,
    t1_by_aspect: dict[str, dict[str, set[str]]],
    memberships: dict[str, dict[str, set[str]]],
) -> dict[str, Path]:
    written: dict[str, Path] = {}
    combined_truth: list[str] = []
    combined_typex: list[str] = []
    for prefix, (official_aspect, _) in ASPECTS.items():
        for category in ("type1", "type2", "typex"):
            path = root / "lists" / f"{prefix}_all_{category}.txt"
            values = sorted(memberships[prefix][category])
            _write_lines(path, values)
            written[f"{prefix}_{category}"] = path
        combined_typex.extend(sorted(memberships[prefix]["typex"]))

        truth_rows = []
        for protein_id in sorted(memberships[prefix]["typex"]):
            propagated = _propagate(go, t1_by_aspect[prefix][protein_id])
            for term in sorted(_leaf_terms(go, propagated)):
                truth_rows.append(f"{protein_id}\t{term}")
        truth_path = root / "groundtruth" / f"leafonly_{official_aspect}.txt"
        _write_lines(truth_path, truth_rows)
        written[f"leafonly_{official_aspect}"] = truth_path
        combined_truth.extend(truth_rows)

    xxo = root / "lists" / "xxo_all_typex.txt"
    _write_lines(xxo, combined_typex)
    written["xxo_typex"] = xxo
    aggregate = root / "groundtruth" / "leafonly_all.txt"
    _write_lines(aggregate, sorted(set(combined_truth)))
    written["leafonly_all"] = aggregate
    return written


def _read_set(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def _read_counter(path: Path) -> Counter[str]:
    return Counter(line.strip() for line in path.read_text().splitlines() if line.strip())


def _set_metrics(generated: set[str], official: set[str]) -> dict[str, object]:
    overlap = generated & official
    union = generated | official
    return {
        "generated": len(generated),
        "official": len(official),
        "overlap": len(overlap),
        "generated_only": len(generated - official),
        "official_only": len(official - generated),
        "precision": len(overlap) / len(generated) if generated else None,
        "recall": len(overlap) / len(official) if official else None,
        "jaccard": len(overlap) / len(union) if union else 1.0,
    }


def _multiset_metrics(
    generated: Counter[str], official: Counter[str]
) -> dict[str, object]:
    overlap = generated & official
    union = generated | official
    generated_count = sum(generated.values())
    official_count = sum(official.values())
    overlap_count = sum(overlap.values())
    return {
        "generated": generated_count,
        "official": official_count,
        "overlap": overlap_count,
        "generated_only": sum((generated - official).values()),
        "official_only": sum((official - generated).values()),
        "precision": overlap_count / generated_count if generated_count else None,
        "recall": overlap_count / official_count if official_count else None,
        "jaccard": overlap_count / sum(union.values()) if union else 1.0,
    }


def _expanded_counter(values: Counter[str]) -> list[str]:
    return [value for value in sorted(values) for _ in range(values[value])]


def _truth_by_protein(rows: set[str]) -> dict[str, set[str]]:
    truth: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        fields = row.split("\t")
        if len(fields) != 2:
            raise ValueError(f"Expected two-column leaf truth, got: {row!r}")
        truth[fields[0]].add(fields[1])
    return dict(truth)


def _compare_release(
    generated_root: Path,
    official_root: Path,
    difference_root: Path,
) -> dict[str, dict]:
    comparisons: dict[str, dict] = {}
    difference_root.mkdir(parents=True, exist_ok=False)
    for prefix, (official_aspect, _) in ASPECTS.items():
        for category in ("type1", "type2", "typex"):
            name = f"{prefix}_{category}"
            generated = _read_set(
                generated_root / "lists" / f"{prefix}_all_{category}.txt"
            )
            official = _read_set(
                official_root / "lists" / f"{prefix}_all_{category}.txt"
            )
            comparisons[name] = _set_metrics(generated, official)
            _write_lines(difference_root / f"{name}_generated_only.txt", sorted(generated - official))
            _write_lines(difference_root / f"{name}_official_only.txt", sorted(official - generated))
        name = f"leafonly_{official_aspect}"
        generated = _read_set(
            generated_root / "groundtruth" / f"leafonly_{official_aspect}.txt"
        )
        official = _read_set(
            official_root / "groundtruth" / f"leafonly_{official_aspect}.txt"
        )
        metrics = _set_metrics(generated, official)
        generated_truth = _truth_by_protein(generated)
        official_truth = _truth_by_protein(official)
        protein_union = set(generated_truth) | set(official_truth)
        exact = sum(
            generated_truth.get(protein_id, set()) == official_truth.get(protein_id, set())
            for protein_id in protein_union
        )
        metrics.update({
            "generated_proteins": len(generated_truth),
            "official_proteins": len(official_truth),
            "exact_protein_truth": exact,
            "exact_protein_truth_fraction": exact / len(protein_union) if protein_union else 1.0,
        })
        comparisons[name] = metrics
        _write_lines(difference_root / f"{name}_generated_only.txt", sorted(generated - official))
        _write_lines(difference_root / f"{name}_official_only.txt", sorted(official - generated))
    generated_xxo = _read_counter(generated_root / "lists" / "xxo_all_typex.txt")
    official_xxo = _read_counter(official_root / "lists" / "xxo_all_typex.txt")
    comparisons["xxo_typex"] = _multiset_metrics(generated_xxo, official_xxo)
    _write_lines(
        difference_root / "xxo_typex_generated_only.txt",
        _expanded_counter(generated_xxo - official_xxo),
    )
    _write_lines(
        difference_root / "xxo_typex_official_only.txt",
        _expanded_counter(official_xxo - generated_xxo),
    )
    return comparisons


def _comparison_markdown(comparisons: dict[str, dict]) -> str:
    lines = [
        "# CAFA3 organiser-replay comparison",
        "",
        "Aspect-list and leaf-truth comparisons are order-independent sets. The",
        "combined xxo list is a multiset because cross-aspect proteins legitimately",
        "repeat. Precision uses generated output as the denominator; recall uses the",
        "released CAFA3 output as the denominator.",
        "",
        "| Artefact | Generated | Official | Overlap | Generated only | Official only | Precision | Recall | Jaccard |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, values in sorted(comparisons.items()):
        def fmt(value: object) -> str:
            return "n/a" if value is None else f"{value:.6f}" if isinstance(value, float) else str(value)

        lines.append(
            f"| {name} | {values['generated']} | {values['official']} | "
            f"{values['overlap']} | {values['generated_only']} | "
            f"{values['official_only']} | {fmt(values['precision'])} | "
            f"{fmt(values['recall'])} | {fmt(values['jaccard'])} |"
        )
    return "\n".join(lines) + "\n"


def _counter_payload(result: AnnotationLoadResult) -> dict[str, object]:
    return {
        "counters": dict(sorted(result.counters.items())),
        "evidence_counts": dict(sorted(result.evidence_counts.items())),
        "taxon_counts": dict(sorted(result.taxon_counts.items())),
        "unmapped_terms": dict(sorted(result.unmapped_terms.items())),
        "out_of_benchmark_terms": dict(sorted(result.out_of_benchmark_terms.items())),
    }


def _input_evidence(config: OrganizerReplayConfig) -> list[dict[str, object]]:
    named_paths: list[tuple[str, Path]] = []
    named_paths.extend(("uniprot_t0", path) for path in config.uniprot_t0)
    named_paths.extend(("uniprot_t1", path) for path in config.uniprot_t1)
    named_paths.extend(("official_target_fasta", path) for path in config.official_target_fastas)
    named_paths.extend(
        ("official_target_mapping", path)
        for path in sorted(config.official_target_mapping_dir.iterdir())
        if path.is_file()
    )
    named_paths.extend([
        ("goa_t0", config.goa_t0),
        ("goa_t1", config.goa_t1),
        ("go_obo", config.go_obo),
        ("go_obo_t0", config.ontology_t0),
        ("go_obo_t1", config.ontology_t1),
    ])
    if config.official_benchmark_dir is not None:
        official_root = _resolve_official_root(config.official_benchmark_dir)
        named_paths.extend(
            ("official_output_oracle", path)
            for path in sorted((official_root / "lists").glob("*_all_type[12x].txt"))
        )
        named_paths.append(
            ("official_output_oracle", official_root / "lists" / "xxo_all_typex.txt")
        )
        named_paths.extend(
            ("official_output_oracle", path)
            for path in sorted((official_root / "groundtruth").glob("leafonly_*.txt"))
        )
    rows = []
    for role, path in named_paths:
        rows.append({
            "role": role,
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path) if config.write_input_checksums else None,
        })
    return rows


def _framework_provenance() -> dict[str, object]:
    root = Path(__file__).resolve().parents[4]
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
        ).strip()
        branch = subprocess.check_output(
            ["git", "-C", str(root), "branch", "--show-current"], text=True
        ).strip()
        dirty = bool(subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain"], text=True
        ).strip())
    except (OSError, subprocess.CalledProcessError):
        return {"available": False, "commit": None, "branch": None, "dirty": None}
    return {"available": True, "commit": commit, "branch": branch, "dirty": dirty}


def run_organizer_replay(config: OrganizerReplayConfig) -> dict[str, Path]:
    _validate(config)
    with_rels, relationship_types = _relationship_types(config.relationship_policy)
    benchmark_go = Ontology(
        config.go_obo,
        with_rels=with_rels,
        relationship_types=relationship_types,
    )
    t0_go = benchmark_go if config.ontology_t0 == config.go_obo else Ontology(
        config.ontology_t0,
        with_rels=with_rels,
        relationship_types=relationship_types,
    )
    t1_go = benchmark_go if config.ontology_t1 == config.go_obo else Ontology(
        config.ontology_t1,
        with_rels=with_rels,
        relationship_types=relationship_types,
    )
    t0_targets, t1_targets, target_rows = _load_targets(config)
    t0_result, t1_result = _load_annotations(
        config, benchmark_go, t0_go, t1_go, t0_targets, t1_targets
    )
    (
        t0_direct_by_aspect,
        t1_direct_by_aspect,
        t0_by_aspect,
        t1_by_aspect,
        memberships,
        flow_rows,
    ) = _build_membership(benchmark_go, t0_result.annotations, t1_result.annotations)

    config.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(
        prefix=f".{config.output_dir.name}.staging-",
        dir=config.output_dir.parent,
    ))
    try:
        release_root = staging / "benchmark20171115"
        written = _write_release(release_root, benchmark_go, t1_by_aspect, memberships)
        reports = staging / "reports"
        reports.mkdir()
        _write_tsv(
            reports / "stage_00_official_target_mapping.tsv",
            (
                "snapshot", "target_id", "source_identifiers", "mapping_files",
                "taxon_id", "sequence_length", "status", "reason", "mapping_method",
                "uniprot_accession", "present_in_snapshot",
            ),
            target_rows,
        )
        _write_annotation_stage(
            reports / "stage_01_t0_filtered_direct_annotations.tsv",
            t0_direct_by_aspect,
        )
        _write_annotation_stage(
            reports / "stage_02_t1_gained_direct_annotations.tsv",
            t1_direct_by_aspect,
        )
        _write_annotation_stage(
            reports / "stage_03_t0_after_protein_binding_filter.tsv",
            t0_by_aspect,
        )
        _write_annotation_stage(
            reports / "stage_04_t1_after_protein_binding_filter.tsv",
            t1_by_aspect,
        )
        _write_tsv(
            reports / "stage_05_benchmark_membership_flow.tsv",
            (
                "protein", "aspect", "namespace", "t0_any", "t0_aspect",
                "t1_gained_aspect", "protein_binding_only_removed_t0",
                "protein_binding_only_removed_t1", "category",
            ),
            flow_rows,
        )

        comparisons: dict[str, dict] = {}
        if config.official_benchmark_dir is not None:
            official_root = _resolve_official_root(config.official_benchmark_dir)
            comparisons = _compare_release(
                release_root,
                official_root,
                reports / "official_comparison_differences",
            )
            _atomic_json(reports / "official_comparison.json", comparisons)
            (reports / "official_comparison.md").write_text(
                _comparison_markdown(comparisons)
            )

        membership_counts = {
            prefix: {
                category: len(memberships[prefix][category])
                for category in ("type1", "type2", "typex")
            }
            for prefix in ASPECTS
        }
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "experiment": "cafa3-organizer-replay",
            "boundary": "raw historical databases to CAFA3 benchmark release artefacts",
            "explicitly_excluded": [
                "DeepGOPlus", "TEMPROT", "nine PFP CSVs", "embeddings", "PFP training",
            ],
            "policy": {
                "evidence_codes": sorted(config.evidence_codes),
                "t0": config.t0_cutoff,
                "t1": config.t1_cutoff,
                "t1_endpoint_policy": config.t1_endpoint_policy,
                "exclude_t1_rows_assigned_on_or_before_t0": True,
                "knowledge_classes": "CAFA2 type1 global-NK and type2 ontology-LK",
                "protein_binding_policy": "CAFA2 pre-classification removal at t0 and t1",
                "relationship_policy": config.relationship_policy,
                "not_qualified_rows": "excluded by framework GAF parser",
            },
            "evidence_status": {
                "eight_code_policy": "official CAFA3 final README",
                "assigned_date_backfill_removal": "official CAFA3 final README",
                "protein_binding_removal": "official CAFA3 final README and CAFA2 code",
                "type1_type2_algebra": "CAFA2 code",
                "is_a_part_of_default": "CAFA2 pfp_ontbuild default",
                "not_handling": "framework policy; exact final CAFA3 implementation not recovered",
                "exact_historical_goa_bytes": "caller supplied; not claimed to be organizer-private bytes",
                "final_cafa3_glue_and_manual_curation": "not recovered",
            },
            "target_count": len(t0_targets.records),
            "membership_counts": membership_counts,
            "t0_load": _counter_payload(t0_result),
            "t1_load": _counter_payload(t1_result),
            "official_comparison_present": bool(comparisons),
            "inputs": _input_evidence(config),
            "framework_provenance": _framework_provenance(),
            "reference_sources": {
                "cafa2_repository": {
                    "url": "https://github.com/yuxjiang/CAFA2",
                    "inspected_commit": "c56b9e7247f3c065c6ac71e14ec758b3839d835b",
                },
                "public_cafa_benchmark_repository": {
                    "url": "https://github.com/nguyenngochuy91/CAFA_benchmark",
                    "inspected_commit": "a6b6a09c84b002e6932e7f175c7185127a5d3ccd",
                },
                "cafa3_supplementary_dataset": {
                    "doi": "10.6084/m9.figshare.8135393.v3",
                    "archive_md5": "2eae900f6f3b60228fb99771576c3a12",
                },
            },
        }
        _atomic_json(reports / "replay_manifest.json", manifest)

        output_hashes = []
        for path in sorted(staging.rglob("*")):
            if path.is_file() and path.name != "REPLAY_COMPLETE.json":
                output_hashes.append({
                    "path": str(path.relative_to(staging)),
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                })
        completion = {
            "schema_version": SCHEMA_VERSION,
            "status": "complete",
            "output_files": output_hashes,
            "manifest_sha256": _sha256(reports / "replay_manifest.json"),
        }
        _atomic_json(staging / "REPLAY_COMPLETE.json", completion)
        os.replace(staging, config.output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    return {
        "output": config.output_dir,
        "manifest": config.output_dir / "reports" / "replay_manifest.json",
        "completion": config.output_dir / "REPLAY_COMPLETE.json",
        **{
            name: config.output_dir / path.relative_to(staging)
            for name, path in written.items()
        },
    }


def _discover_targets(root: Path) -> tuple[tuple[Path, ...], Path]:
    fasta_dir = root / "Target files"
    mapping_dir = root / "Mapping files"
    fastas = tuple(sorted(fasta_dir.glob("*.fasta")))
    if not fastas:
        raise FileNotFoundError(f"No target FASTAs found under {fasta_dir}")
    if not mapping_dir.is_dir():
        raise FileNotFoundError(mapping_dir)
    return fastas, mapping_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Replay only the raw-database to CAFA3 benchmark-release boundary; "
            "do not run DeepGOPlus, TEMPROT, embeddings or PFP."
        )
    )
    parser.add_argument("--uniprot-t0", action="append", type=Path, required=True)
    parser.add_argument("--uniprot-t1", action="append", type=Path, required=True)
    parser.add_argument("--goa-t0", type=Path, required=True)
    parser.add_argument("--goa-t1", type=Path, required=True)
    parser.add_argument("--go-obo", type=Path, required=True)
    parser.add_argument("--go-obo-t0", type=Path)
    parser.add_argument("--go-obo-t1", type=Path)
    targets = parser.add_mutually_exclusive_group(required=True)
    targets.add_argument(
        "--official-target-root",
        type=Path,
        help="Extracted CAFA3_targets root containing Target files/ and Mapping files/.",
    )
    targets.add_argument("--official-target-fasta", action="append", type=Path)
    parser.add_argument("--official-target-mapping-dir", type=Path)
    parser.add_argument("--official-benchmark-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--t0-cutoff", default="20170213")
    parser.add_argument("--t1-cutoff", default="20171115")
    parser.add_argument(
        "--t1-endpoint-policy",
        choices=("assigned-date-proxy", "snapshot-membership"),
        default="assigned-date-proxy",
    )
    parser.add_argument(
        "--relationship-policy",
        choices=("is-a-part-of", "is-a", "all"),
        default="is-a-part-of",
    )
    parser.add_argument("--strict-qc", action="store_true")
    parser.add_argument("--no-frozen-source-fallback", action="store_true")
    parser.add_argument("--max-gaf-records", type=int)
    parser.add_argument("--skip-input-checksums", action="store_true")
    return parser


def config_from_args(args: argparse.Namespace) -> OrganizerReplayConfig:
    if args.official_target_root:
        fastas, mapping_dir = _discover_targets(args.official_target_root)
    else:
        fastas = tuple(args.official_target_fasta or ())
        if args.official_target_mapping_dir is None:
            raise SystemExit(
                "--official-target-mapping-dir is required with --official-target-fasta"
            )
        mapping_dir = args.official_target_mapping_dir
    return OrganizerReplayConfig(
        uniprot_t0=tuple(args.uniprot_t0),
        uniprot_t1=tuple(args.uniprot_t1),
        goa_t0=args.goa_t0,
        goa_t1=args.goa_t1,
        go_obo=args.go_obo,
        go_obo_t0=args.go_obo_t0,
        go_obo_t1=args.go_obo_t1,
        official_target_fastas=fastas,
        official_target_mapping_dir=mapping_dir,
        official_benchmark_dir=args.official_benchmark_dir,
        output_dir=args.output_dir,
        t0_cutoff=normalise_gaf_date(args.t0_cutoff) or "20170213",
        t1_cutoff=normalise_gaf_date(args.t1_cutoff) or "20171115",
        t1_endpoint_policy=args.t1_endpoint_policy,
        relationship_policy=args.relationship_policy,
        strict_qc=args.strict_qc,
        allow_frozen_source_fallback=not args.no_frozen_source_fallback,
        max_gaf_records=args.max_gaf_records,
        write_input_checksums=not args.skip_input_checksums,
    )


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    written = run_organizer_replay(config_from_args(args))
    print(f"CAFA3 organiser replay complete: {written['output']}")


if __name__ == "__main__":
    main()
