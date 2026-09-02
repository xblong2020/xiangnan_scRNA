# Code repository archival report

## Final local decision

| Field | Value |
|---|---|
| REPOSITORY_ARCHIVAL_STATUS | READY_FOR_MANUAL_ARCHIVAL |
| PENDING_REPOSITORY_ARCHIVAL | TRUE |
| STAGE22_REPOSITORY_BLOCKER | PENDING_EXTERNAL_ACTION |
| RELEASE_TAG | v1.0.0 |
| RELEASE_COMMIT | 1c0303049ba629b9b986cb9a7e088384f15f5a87 |
| CURRENT_HEAD_AT_AUDIT | 50acc2a18866e53d4fdd54deea92cd416ca1aaa0 |
| REMOTE_REPOSITORY_STATUS | NOT_CONFIGURED |
| ZENODO_STATUS | READY_FOR_MANUAL_ARCHIVAL |
| PERMANENT_IDENTIFIER | NULL |

## Decision rationale

The repository has a local Git provenance record and an annotated release tag when shown above. No public remote, archive response, DOI, SWHID, or other permanent identifier is verified locally.
Therefore the repository archival blocker remains pending external action. The local release package is ready for the project owner to publish once, then append the provider verification evidence.

## Scope protection

- No biological analysis was rerun.
- Figure 1–8, Results, and Discussion were not modified.
- Stage19 was not reopened.
- Stage23 was not entered automatically.
- GSE326201 remains Tier 1 exploratory; GSE282701 remains BLOCKED_PROVENANCE_UNRESOLVED; ICGC OS remains ESTIMABLE_BUT_NOT_VALIDATED and supplementary/Extended Data only; Figure 8 remains EXTENDED_DATA_ONLY.
- historical_exact_version_not_recoverable remains explicit.

## Minimum manual action

Create/configure the public remote, push the release branch and the immutable release tag, enable Zenodo or another approved archive, record the returned DOI or other permanent identifier, verify its URL and associated commit, and replace only the explicit pending fields in the closure records.
