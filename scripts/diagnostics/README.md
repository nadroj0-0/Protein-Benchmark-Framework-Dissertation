# Diagnostics

## Audit a PFP working copy

`audit_pfp_working_copy.py` compares an existing local PFP working directory
with a fresh temporary clone of `psipred/PFP`. It does not modify the supplied
working directory, upload anything, or retain the public clone.

Activate the Python environment that was used for PFP, then run:

```bash
python scripts/diagnostics/audit_pfp_working_copy.py /path/to/working/PFP \
  > pfp_working_copy_audit.md
```

The redirection is optional and is the only persistent output. Without it, the
Markdown report is printed to the terminal. A non-zero exit means that the
working copy differs materially from the current public release.

The report separates:

- public tracked files modified or missing locally;
- local Git-tracked files absent from the public release;
- names of local untracked entries, without reading or printing their contents;
- the active Python environment used to run the audit, including a complete
  `pip freeze --all`.

Local-only files are evidence of local development, not automatic evidence
that those files should have been included in the published repository.
