# Repository scope decision

## Included in the code-release scope

- Top-level analysis, workflow, validation, and figure-source-data scripts under scripts/.
- Root README.md, CITATION.cff, CODE_AVAILABILITY.md, REPRODUCIBILITY_NOTES.md, and .gitignore.
- The Stage22 archival closure package, including Git provenance, release, remote, archive, manuscript wording, QA, and Stage22 status records.
- The small local metadata-audit builder used to regenerate the closure records from Git/filesystem evidence.

## Excluded from Git tracking

- data/, figures/, metadata/, reports/, tmp/, large sequencing or matrix objects, and stage-generated output/input/quality-check folders.
- Restricted or patient-level data, local virtual environments, Python/R caches, .Rhistory, generated plots, Zotero local databases and temporary exports, OneDrive workspace artifacts, credentials, API keys, and tokens.
- Existing scientific results are preserved in place and are not deleted or rewritten by this closure.

## Staging rule

The release uses an explicit allowlist. No broad git add -A is used. The final manifest records the paths and SHA-256 values of the release scope, while protected files remain outside Git tracking.
