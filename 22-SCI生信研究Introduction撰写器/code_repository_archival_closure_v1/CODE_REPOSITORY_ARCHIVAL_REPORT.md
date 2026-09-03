# Code repository archival report

## Final local decision

| Field | Value |
|---|---|
| REPOSITORY_ARCHIVAL_STATUS | OPEN_PENDING_EXTERNAL_ARCHIVAL |
| PENDING_REPOSITORY_ARCHIVAL | TRUE |
| STAGE22_REPOSITORY_BLOCKER | PENDING_ZENODO_ARCHIVAL |
| RELEASE_TAG | v1.0.0 |
| RELEASE_COMMIT | 1c0303049ba629b9b986cb9a7e088384f15f5a87 |
| CURRENT_HEAD_AT_AUDIT | d1b9d7e39a5e2d563b4ad286b770c8fef350c360 |
| REMOTE_REPOSITORY_STATUS | CONFIGURED_AND_VERIFIED |
| GITHUB_RELEASE_STATUS | PUBLISHED |
| GITHUB_RELEASE_URL | https://github.com/xblong2020/xiangnan_scRNA/releases/tag/v1.0.0 |
| ZENODO_STATUS | MANUAL_ACTION_REQUIRED |
| PERMANENT_IDENTIFIER | NULL |

## Decision rationale

The repository, GitHub Release, and immutable tag states are recorded from the supplied local and remote evidence. Zenodo permanent-identifier verification remains the independent archival gate.
The final status is driven by the explicit publication-state inputs and retains pending values whenever DOI evidence is absent.

## Scope protection

- No biological analysis was rerun.
- Figure 1–8, Results, and Discussion were not modified.
- Stage19 was not reopened.
- Stage23 was not entered automatically.
- GSE326201 remains Tier 1 exploratory; GSE282701 remains BLOCKED_PROVENANCE_UNRESOLVED; ICGC OS remains ESTIMABLE_BUT_NOT_VALIDATED and supplementary/Extended Data only; Figure 8 remains EXTENDED_DATA_ONLY.
- historical_exact_version_not_recoverable remains explicit.

## Minimum manual action

Enable the verified GitHub repository in Zenodo, obtain the version DOI for v1.0.0, verify DOI resolution and release provenance, and replace only the explicit pending fields in the closure records.
