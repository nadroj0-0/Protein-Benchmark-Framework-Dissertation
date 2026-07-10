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

- GOA-declared ontology version **2026-06-15**
- Public GO product directory **2026-06-19**, whose `go-basic.obo` declares
  `data-version: releases/2026-06-15`

---

## Reproducing this benchmark in the future

The download script checks the mutable UniProt and GOA endpoints before using
them and fails if they have advanced. The GO ontology URLs are fixed archived
release URLs.

GOA 225 declares GO version 2025-03-07. No standalone 2025-03-07 product
directory is retained in the public GO release bucket. Benchmark generation
therefore freezes predictions to 2025-02-06 and uses 2025-03-16 only for source
GO-ID resolution; the resulting exclusions are reported.

This ensures that the exact evaluation snapshot used on **9 July 2026** can be
recreated in the future.
