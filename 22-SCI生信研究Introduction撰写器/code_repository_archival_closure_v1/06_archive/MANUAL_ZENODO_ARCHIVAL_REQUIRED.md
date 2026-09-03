# Manual Zenodo archival required

Current status:

- GITHUB_RELEASE_STATUS = PUBLISHED
- GITHUB_RELEASE_TAG = v1.0.0
- GITHUB_RELEASE_URL = https://github.com/xblong2020/xiangnan_scRNA/releases/tag/v1.0.0
- ZENODO_STATUS = MANUAL_ACTION_REQUIRED
- PERMANENT_IDENTIFIER = null
- PENDING_REPOSITORY_ARCHIVAL = TRUE
- STAGE22_REPOSITORY_BLOCKER = PENDING_ZENODO_ARCHIVAL
- FINAL_GATE = MANUAL_ZENODO_ACTION_REQUIRED

Zenodo public API queries for xiangnan_scRNA and xblong2020 returned zero records in this audit. No Zenodo account token or authenticated archive endpoint is available to this workspace.

Verification timestamp: 2026-09-03T09:16:18+08:00

## Required manual action

1. Sign in to Zenodo with the project account.
2. Open the GitHub integration or Linked accounts settings.
3. Connect GitHub and find xblong2020/xiangnan_scRNA.
4. Enable repository archival.
5. Confirm that the normal GitHub Release v1.0.0 is published.
6. Check that Zenodo starts processing the release.
7. Open the generated Zenodo record.
8. Retrieve the real version DOI and Zenodo record URL.
9. Verify that the record is public, version v1.0.0, and associated with commit 1c0303049ba629b9b986cb9a7e088384f15f5a87.
10. Provide the DOI and record URL back to the workspace for three-way DOI verification.

Do not manually construct a DOI. Keep PERMANENT_IDENTIFIER null until a real Zenodo record resolves through https://doi.org/ and its metadata/provenance are independently verified.
