from __future__ import annotations

from collections import Counter
from pathlib import Path
import logging

import numpy as np
import pandas as pd

from .config import BuildConfig, NAMESPACE_TO_PREFIX, PREFIX_TO_NAMESPACE
from .goa import load_annotation_map
from .ontology import Ontology
from .parsers import iter_uniprot

LOGGER = logging.getLogger(__name__)


def load_proteins(paths: tuple[Path, ...], target_taxa: frozenset[str], reviewed_only: bool) -> dict[str, str]:
    proteins: dict[str, str] = {}
    for path in paths:
        LOGGER.info("Loading UniProt sequences from %s", path)
        for rec in iter_uniprot(path):
            if target_taxa and rec.taxon_id not in target_taxa:
                continue
            if reviewed_only and rec.reviewed is False:
                continue
            proteins[rec.protein_id] = rec.sequence
            for accession in rec.accessions:
                proteins.setdefault(accession, rec.sequence)
    LOGGER.info("Loaded %d unique sequence IDs", len(proteins))
    return proteins


def propagate_annotations(go: Ontology, annotations: set[str]) -> set[str]:
    annots_set = set()
    for go_id in annotations:
        annots_set |= go.get_anchestors(go_id)
    return annots_set


def build_training_dataframe(
    go: Ontology,
    sequences: dict[str, str],
    train_annots: dict[str, set[str]],
) -> tuple[pd.DataFrame, Counter]:
    proteins = []
    seqs = []
    prop_annotations = []
    cnt = Counter()
    skipped_missing_sequence = 0
    skipped_empty_propagation = 0

    for prot_id, annots in train_annots.items():
        if prot_id not in sequences:
            skipped_missing_sequence += 1
            continue
        annots_set = propagate_annotations(go, annots)
        if not annots_set:
            skipped_empty_propagation += 1
            continue
        proteins.append(prot_id)
        seqs.append(sequences[prot_id])
        prop_annotations.append(annots_set)
        for go_id in annots_set:
            cnt[go_id] += 1

    df = pd.DataFrame({
        "proteins": proteins,
        "sequences": seqs,
        "annotations": prop_annotations,
    })
    if skipped_missing_sequence:
        LOGGER.warning(
            "Dropped %d t0-annotated proteins from training because no t0 sequence was loaded",
            skipped_missing_sequence,
        )
    if skipped_empty_propagation:
        LOGGER.warning(
            "Dropped %d t0-annotated proteins from training because annotation propagation returned no GO terms",
            skipped_empty_propagation,
        )
    return df, cnt


def build_test_dataframe(
    go: Ontology,
    sequences: dict[str, str],
    t0_annots: dict[str, set[str]],
    t1_annots: dict[str, set[str]],
    train_proteins: set[str],
) -> pd.DataFrame:
    proteins = []
    seqs = []
    prop_annotations = []
    skipped_train_id = 0
    skipped_missing_sequence = 0
    skipped_no_gain = 0
    skipped_empty_propagation = 0

    for prot_id, annots_t1 in t1_annots.items():
        if prot_id in train_proteins:
            skipped_train_id += 1
            continue
        if prot_id not in sequences:
            skipped_missing_sequence += 1
            continue
        gained = annots_t1 - t0_annots.get(prot_id, set())
        if not gained:
            skipped_no_gain += 1
            continue
        annots_set = propagate_annotations(go, gained)
        if not annots_set:
            skipped_empty_propagation += 1
            continue
        proteins.append(prot_id)
        seqs.append(sequences[prot_id])
        prop_annotations.append(annots_set)

    df = pd.DataFrame({
        "proteins": proteins,
        "sequences": seqs,
        "annotations": prop_annotations,
    })
    LOGGER.info(
        "Test candidate filtering: kept=%d skipped_train_or_valid_id=%d "
        "skipped_missing_t1_sequence=%d skipped_no_gained_terms=%d skipped_empty_propagation=%d",
        len(df),
        skipped_train_id,
        skipped_missing_sequence,
        skipped_no_gain,
        skipped_empty_propagation,
    )
    return df


def make_terms_dataframe(cnt: Counter, min_count: int) -> pd.DataFrame:
    # Mirrors DeepGOPlus cafa3_data.py: filter propagated training terms by
    # min_count and store a single-column DataFrame named "terms".
    terms = []
    for key, val in cnt.items():
        if val >= min_count:
            terms.append(key)
    return pd.DataFrame({"terms": terms})


def split_train_valid(df: pd.DataFrame, split: float = 0.9, seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Mirrors DeepGOPlus deepgoplus.py load_data split block.
    n = len(df)
    index = np.arange(n)
    train_n = int(n * split)
    np.random.seed(seed=seed)
    np.random.shuffle(index)
    train_df = df.iloc[index[:train_n]]
    valid_df = df.iloc[index[train_n:]]
    return train_df, valid_df


def terms_for_namespace(go: Ontology, terms_df: pd.DataFrame, namespace: str) -> list[str]:
    terms = []
    for term in terms_df.terms.tolist():
        if term in go.ont and go.get_namespace(term) == namespace:
            terms.append(term)
    return terms


def preprocess_for_terms(go: Ontology, df: pd.DataFrame, list_terms: list[str]) -> pd.DataFrame:
    # Closely follows TEMPROT src/<ontology>/dataset.py preprocess().
    f_df = {}
    f_df["proteins"] = []
    f_df["sequences"] = []
    term_index = {term: idx for idx, term in enumerate(list_terms)}

    for i in list_terms:
        f_df[i] = []

    for i in range(len(df)):
        is_in_dataset = False
        actual_terms = [0 for _ in range(len(list_terms))]
        protein, sequence, annotation = df.iloc[i, :3].values

        for term in list(annotation):
            if term in term_index:
                is_in_dataset = True
                actual_terms[term_index[term]] = 1
                for ant in go.get_anchestors(term):
                    if ant in term_index:
                        actual_terms[term_index[ant]] = 1

        if is_in_dataset:
            f_df["proteins"].append(protein)
            f_df["sequences"].append(sequence)

            for j in range(len(list_terms)):
                f_df[list_terms[j]].append(actual_terms[j])

    return pd.DataFrame(f_df)


def export_pfp_csvs(
    go: Ontology,
    train_df: pd.DataFrame,
    valid_df: pd.DataFrame,
    test_df: pd.DataFrame,
    terms_df: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    for prefix, namespace in PREFIX_TO_NAMESPACE.items():
        list_terms = terms_for_namespace(go, terms_df, namespace)
        LOGGER.info("Exporting %s with %d GO terms", prefix, len(list_terms))

        n_train = preprocess_for_terms(go, train_df, list_terms)
        n_val = preprocess_for_terms(go, valid_df, list_terms)
        n_test = preprocess_for_terms(go, test_df, list_terms)

        # Same sequence-duplicate removal policy as TEMPROT dataset.py.
        if not n_test.empty:
            n_train = n_train[~n_train.sequences.isin(n_test.sequences.values)]
            n_val = n_val[~n_val.sequences.isin(n_test.sequences.values)]
        if not n_train.empty:
            n_val = n_val[~n_val.sequences.isin(n_train.sequences.values)]

        outputs = {
            "training": n_train,
            "validation": n_val,
            "test": n_test,
        }
        for split, split_df in outputs.items():
            path = output_dir / f"{prefix}-{split}.csv"
            split_df.to_csv(path, index=False)
            written[f"{prefix}-{split}"] = path
    return written


def build_benchmark(config: BuildConfig) -> dict[str, Path]:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")

    LOGGER.info("Loading GO ontology from %s", config.go_obo)
    go = Ontology(config.go_obo, with_rels=config.include_rels)

    seq_t0 = load_proteins(config.uniprot_t0, config.target_taxa, config.reviewed_only)
    seq_t1 = load_proteins(config.uniprot_t1, config.target_taxa, config.reviewed_only)

    t0_annotation_universe = set(seq_t0) | set(seq_t1)
    t1_annotation_universe = set(seq_t1)

    LOGGER.info("Loading t0 GOA annotations from %s", config.goa_t0)
    t0_annots = load_annotation_map(
        config.goa_t0,
        evidence_codes=config.evidence_codes,
        target_taxa=config.target_taxa,
        max_records=config.max_gaf_records,
        allowed_proteins=t0_annotation_universe,
    )
    LOGGER.info("Loading t1 GOA annotations from %s", config.goa_t1)
    t1_annots = load_annotation_map(
        config.goa_t1,
        evidence_codes=config.evidence_codes,
        target_taxa=config.target_taxa,
        max_records=config.max_gaf_records,
        allowed_proteins=t1_annotation_universe,
    )

    LOGGER.info("Building DeepGOPlus-style training dataframe")
    train_all_df, cnt = build_training_dataframe(go, seq_t0, t0_annots)
    LOGGER.info("Training proteins before split: %d", len(train_all_df))

    terms_df = make_terms_dataframe(cnt, config.min_count)
    LOGGER.info("Terms with min_count >= %d: %d", config.min_count, len(terms_df))

    train_df, valid_df = split_train_valid(train_all_df, split=config.split, seed=config.seed)
    LOGGER.info("Train/valid split: %d/%d", len(train_df), len(valid_df))

    LOGGER.info("Building temporal gained-annotation test dataframe")
    test_df = build_test_dataframe(go, seq_t1, t0_annots, t1_annots, set(train_all_df["proteins"]))
    LOGGER.info("Test proteins: %d", len(test_df))

    written: dict[str, Path] = {}
    if config.write_intermediates:
        paths = {
            "train_data": config.output_dir / "train_data.pkl",
            "train_data_train": config.output_dir / "train_data_train.pkl",
            "train_data_valid": config.output_dir / "train_data_valid.pkl",
            "test_data": config.output_dir / "test_data.pkl",
            "terms": config.output_dir / "terms.pkl",
        }
        train_all_df.to_pickle(paths["train_data"])
        train_df.to_pickle(paths["train_data_train"])
        valid_df.to_pickle(paths["train_data_valid"])
        test_df.to_pickle(paths["test_data"])
        terms_df.to_pickle(paths["terms"])
        written.update(paths)

    csv_paths = export_pfp_csvs(go, train_df, valid_df, test_df, terms_df, config.output_dir)
    written.update(csv_paths)
    return written
