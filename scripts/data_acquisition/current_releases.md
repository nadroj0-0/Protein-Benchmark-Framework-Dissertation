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

- UniProtKB Release **2026_01**
- Official release date: **28 January 2026**

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

- Gene Ontology Release **2026-06-19**

---

## Reproducing this benchmark in the future

If the upstream **current** endpoints have advanced beyond the versions listed
above, they should be replaced with the corresponding archived snapshot URLs
before rerunning the download script.

This ensures that the exact evaluation snapshot used on **9 July 2026** can be
recreated in the future.