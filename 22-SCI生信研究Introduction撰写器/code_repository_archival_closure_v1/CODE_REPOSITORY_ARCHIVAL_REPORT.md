# Code repository archival report

## Final local decision

| Field | Value |
|---|---|
| REPOSITORY_ARCHIVAL_STATUS | MANUAL_GITHUB_RELEASE_REQUIRED |
| PENDING_REPOSITORY_ARCHIVAL | TRUE |
| STAGE22_REPOSITORY_BLOCKER | PENDING_GITHUB_RELEASE |
| RELEASE_TAG | v1.0.0 |
| RELEASE_COMMIT | 1c0303049ba629b9b986cb9a7e088384f15f5a87 |
| CURRENT_HEAD_AT_AUDIT | a7ee1edb8c139e9dca94a3e3702c750701e019e0 |
| REMOTE_REPOSITORY_STATUS | CONFIGURED_AND_VERIFIED |
| GITHUB_REPOSITORY_URL | https://github.com/xblong2020/xiangnan_scRNA |
| GITHUB_MAIN_COMMIT | a7ee1edb8c139e9dca94a3e3702c750701e019e0 |
| GITHUB_RELEASE_STATUS | MANUAL_ACTION_REQUIRED |
| ZENODO_STATUS | WAITING_FOR_GITHUB_RELEASE |
| PERMANENT_IDENTIFIER | NULL |

## Decision rationale

The public GitHub repository, main branch, and annotated v1.0.0 tag are verified. The GitHub Release object remains a manual action because the GitHub CLI is unavailable. No DOI, SWHID, or other permanent identifier is verified locally.
Therefore the repository archival blocker remains pending GitHub Release publication and subsequent Zenodo archival.

## Scope protection

- No biological analysis was rerun.
- Figure 1–8, Results, and Discussion were not modified.
- Stage19 was not reopened.
- Stage23 was not entered automatically.
- GSE326201 remains Tier 1 exploratory; GSE282701 remains BLOCKED_PROVENANCE_UNRESOLVED; ICGC OS remains ESTIMABLE_BUT_NOT_VALIDATED and supplementary/Extended Data only; Figure 8 remains EXTENDED_DATA_ONLY.
- historical_exact_version_not_recoverable remains explicit.

## Minimum manual action

Publish the normal GitHub Release v1.0.0 from the existing public repository, enable Zenodo or another approved archive, record the returned DOI or other permanent identifier, verify its URL and associated commit, and replace only the explicit pending fields in the closure records.
