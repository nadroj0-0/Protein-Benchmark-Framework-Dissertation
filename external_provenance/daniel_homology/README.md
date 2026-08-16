# Daniel Buchan Homology Inputs

This directory preserves the small source-provenance files for the independent
supervisor-generated homology partitions used in the dissertation.

- `run_mmseqs_clustering.sh` is the exact script supplied by Daniel Buchan.
- `handoff_state.py` is the run-contract helper required by that script.
- `identity_*_SOURCE_PROVENANCE.json` records the received assignment hashes,
  source population, clustering settings, and validation facts for 30%, 25%,
  20%, and 15% identity.

The large assignment files and their source archives are external data and are
not stored in Git. Their SHA-256 digests are recorded in the provenance JSON.
These files document an external comparison stream; they are not the framework's
own MMseqs2 implementation.

The exact archive-to-gzip import used for the 25%, 20%, and 15% handoff is
preserved in `scripts/data_acquisition/import_daniel_clusters.sh`. It is a
historical machine-bound import record; the portable submission workflow may
replace its absolute storage bindings while retaining the checksums and
validation contract recorded here.
