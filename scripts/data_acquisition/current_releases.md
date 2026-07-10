# Current Release Resolution

**Benchmark generation date:** 9 July 2026

The evaluation (t1) component of this project intentionally downloads data from
the upstream **current** release endpoints rather than archived snapshots.

Because these endpoints change over time, this document records exactly what
"current" referred to on **9 July 2026**.

---

## UniProt

Current endpoint:

https://ftp.uniprot.org/pub/databases/uniprot/current_release/

Resolved release on **9 July 2026**:

- UniProtKB Release **2026_02**
- Official release date: **10 June 2026**

---

## GOA

Current endpoint:

https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/

Resolved release on **9 July 2026**:

- GOA UniProt Release **234**
- Official release date: **17 June 2026**

---

## Gene Ontology

Current endpoint:

https://current.geneontology.org/

Resolved release on **9 July 2026**:

- GOA-declared ontology version **2026-06-15** (used to normalise GOA 234)
- The independent current GO endpoint was release **2026-06-19** on 9 July;
  it is not substituted for the ontology version declared by the GAF.

---

## Reproducing this benchmark in the future

The download script checks the mutable UniProt and GOA endpoints before using
them and fails if they have advanced. The GO ontology URLs are fixed archived
release URLs.

This ensures that the exact evaluation snapshot used on **9 July 2026** can be
recreated in the future.
