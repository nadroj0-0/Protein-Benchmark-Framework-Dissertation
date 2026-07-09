from __future__ import annotations

from collections import Counter
from pathlib import Path
import logging
import shutil

import numpy as np
import pandas as pd

from .config import BuildConfig, NAMESPACE_TO_PREFIX, PREFIX_TO_NAMESPACE
from .goa import load_annotation_map
from .ontology import Ontology
from .parsers import iter_fasta, iter_uniprot

LOGGER = logging.getLogger(__name__)

DEEPGOPLUS_PICKLE_FILES = {
    "train_data": "train_data.pkl",
    "train_data_train": "train_data_train.pkl",
    "train_data_valid": "train_data_valid.pkl",
    "test_data": "test_data.pkl",
    "terms": "terms.pkl",
}


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


def make_deepgoplus_terms_dataframe(cnt: Counter, min_count: int) -> pd.DataFrame:
    """Mirror DeepGOPlus cafa3_data.py term filtering/order.

    The original groups by the prefix before ":" even though all GO terms share
    the same "GO" prefix. Keeping the shape here makes the historical pickle
    generation path easier to compare directly against cafa3_data.py.
    """
    res: dict[str, list[str]] = {}
    for key, val in cnt.items():
        if val >= min_count:
            ont = key.split(":")[0]
            if ont not in res:
                res[ont] = []
            res[ont].append(key)

    terms = []
    for _, val in res.items():
        terms += val
    return pd.DataFrame({"terms": terms})


def load_deepgoplus_annotation_file(path: Path) -> dict[str, set[str]]:
    """Read DeepGOPlus/CAFA two-column tab annotation files."""
    annots: dict[str, set[str]] = {}
    with open(path, "r") as handle:
        for line in handle:
            it = line.strip().split("\t")
            if len(it) < 2:
                continue
            prot_id = it[0]
            if prot_id not in annots:
                annots[prot_id] = set()
            go_id = it[1]
            annots[prot_id].add(go_id)
    return annots


def build_deepgoplus_dataframe_from_fasta(
    go: Ontology,
    sequences_file: Path,
    annots: dict[str, set[str]],
    count_terms: bool,
) -> tuple[pd.DataFrame, Counter]:
    """Mirror the dataframe construction in DeepGOPlus cafa3_data.py."""
    proteins = []
    sequences = []
    annotations = []
    cnt = Counter()

    for prot_info, sequence in iter_fasta(sequences_file):
        prot_id = prot_info.split()[0]
        if prot_id in annots:
            proteins.append(prot_id)
            sequences.append(sequence)
            annotations.append(annots[prot_id])

    prop_annotations = []
    for direct_annots in annotations:
        annots_set = set()
        for go_id in direct_annots:
            annots_set |= go.get_anchestors(go_id)
        prop_annotations.append(annots_set)
        if count_terms:
            for go_id in annots_set:
                cnt[go_id] += 1

    df = pd.DataFrame({
        "proteins": proteins,
        "sequences": sequences,
        "annotations": prop_annotations,
    })
    return df, cnt


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


def require_deepgoplus_file(deepgoplus_dir: Path, filename: str) -> Path:
    path = deepgoplus_dir / filename
    if not path.is_file():
        raise FileNotFoundError(f"Missing required DeepGOPlus file: {path}")
    return path


def export_from_deepgoplus_pickles(
    deepgoplus_dir: Path,
    go_obo: Path,
    output_dir: Path,
    include_rels: bool = True,
    write_intermediates: bool = True,
) -> dict[str, Path]:
    """Export PFP-compatible CSVs from released DeepGOPlus/TEMPROT pickles.

    This is the historical validation path. It intentionally starts from the
    same intermediate artefacts used by TEMPROT instead of rebuilding them from
    raw GOA snapshots.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")

    LOGGER.info("Loading GO ontology from %s", go_obo)
    go = Ontology(go_obo, with_rels=include_rels)

    paths = {
        key: require_deepgoplus_file(deepgoplus_dir, filename)
        for key, filename in DEEPGOPLUS_PICKLE_FILES.items()
    }

    LOGGER.info("Loading DeepGOPlus train split from %s", paths["train_data_train"])
    train_df = pd.read_pickle(paths["train_data_train"])
    LOGGER.info("Loading DeepGOPlus validation split from %s", paths["train_data_valid"])
    valid_df = pd.read_pickle(paths["train_data_valid"])
    LOGGER.info("Loading DeepGOPlus test data from %s", paths["test_data"])
    test_df = pd.read_pickle(paths["test_data"])
    LOGGER.info("Loading DeepGOPlus terms from %s", paths["terms"])
    terms_df = pd.read_pickle(paths["terms"])

    LOGGER.info(
        "Loaded DeepGOPlus dataframes: train=%s valid=%s test=%s terms=%s",
        train_df.shape,
        valid_df.shape,
        test_df.shape,
        terms_df.shape,
    )

    written: dict[str, Path] = {}
    if write_intermediates:
        for key, src in paths.items():
            dest = output_dir / src.name
            shutil.copy2(src, dest)
            written[key] = dest

    csv_paths = export_pfp_csvs(go, train_df, valid_df.iloc[:, :3], test_df, terms_df, output_dir)
    written.update(csv_paths)
    return written


def generate_deepgoplus_pickles_from_cafa_files(
    go_obo: Path,
    train_sequences_file: Path,
    train_annotations_file: Path,
    test_sequences_file: Path,
    test_annotations_file: Path,
    output_dir: Path,
    min_count: int = 50,
    include_rels: bool = True,
) -> dict[str, Path]:
    """Generate DeepGOPlus CAFA pickles from official CAFA-style files.

    This mirrors DeepGOPlus cafa3_data.py. It is intentionally narrower than
    the contemporary snapshot builder: it starts from the already-curated CAFA
    training/target/ground-truth files and writes only train_data.pkl,
    test_data.pkl and terms.pkl.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")

    LOGGER.info("Loading GO from %s", go_obo)
    go = Ontology(go_obo, with_rels=include_rels)

    LOGGER.info("Loading training annotations from %s", train_annotations_file)
    train_annots = load_deepgoplus_annotation_file(train_annotations_file)

    LOGGER.info("Loading training sequences from %s", train_sequences_file)
    train_df, cnt = build_deepgoplus_dataframe_from_fasta(
        go=go,
        sequences_file=train_sequences_file,
        annots=train_annots,
        count_terms=True,
    )
    LOGGER.info("Train proteins: %d", len(train_df))

    terms_df = make_deepgoplus_terms_dataframe(cnt, min_count)
    LOGGER.info("Number of terms %d", len(terms_df))

    LOGGER.info("Loading testing annotations from %s", test_annotations_file)
    test_annots = load_deepgoplus_annotation_file(test_annotations_file)

    LOGGER.info("Loading testing sequences from %s", test_sequences_file)
    test_df, _ = build_deepgoplus_dataframe_from_fasta(
        go=go,
        sequences_file=test_sequences_file,
        annots=test_annots,
        count_terms=False,
    )
    LOGGER.info("Test proteins: %d", len(test_df))

    written = {
        "train_data": output_dir / "train_data.pkl",
        "test_data": output_dir / "test_data.pkl",
        "terms": output_dir / "terms.pkl",
    }
    LOGGER.info("Saving training data to %s", written["train_data"])
    train_df.to_pickle(written["train_data"])
    LOGGER.info("Saving terms to %s", written["terms"])
    terms_df.to_pickle(written["terms"])
    LOGGER.info("Saving testing data to %s", written["test_data"])
    test_df.to_pickle(written["test_data"])
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
