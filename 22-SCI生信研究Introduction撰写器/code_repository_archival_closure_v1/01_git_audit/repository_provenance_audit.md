# Repository provenance audit

Generated locally: 2026-09-03T09:34:41+08:00
Audit phase: final

| Field | Value |
|---|---|
| REPOSITORY_EXISTS | TRUE |
| GIT_ROOT | REDACTED_LOCAL_PATH |
| CURRENT_BRANCH | codex/module7-sctenifoldknk |
| HEAD_COMMIT | b024f6aae586b326bd91b8b0c2fc0d51b47a3540 |
| HEAD_STATE | COMMIT |
| REMOTE_EXISTS | TRUE |
| REMOTE_URL | https://github.com/xblong2020/xiangnan_scRNA.git |
| TAG_EXISTS | TRUE |
| LATEST_TAG | v1.0.0 |
| WORKTREE_CLEAN | FALSE |
| ARCHIVAL_READY | FALSE |
| RELEASE_TAG | v1.0.0 |
| RELEASE_COMMIT | 1c0303049ba629b9b986cb9a7e088384f15f5a87 |
| REMOTE_REPOSITORY_STATUS | CONFIGURED_AND_VERIFIED |
| REMOTE_PUSH_STATUS | SUCCESS |

## Evidence interpretation

The release anchor is a local annotated tag only until a public remote and archive response are independently verified.
No environment variable values were read. Secret-bearing files are excluded by the repository policy.
The local repository began this closure with an unborn branch; any commit and tag listed after the release phase are real Git outputs from this closure.
The release scope excludes protected data, generated figures/results, local environments, caches, Zotero local data, and credentials.

## Required boundary checks

| Field | Value |
|---|---|
| large release-scope files >= 50 MiB | 0 |
| sensitive scan | PASS_NO_RELEASE_SCOPE_SECRET_PATTERN |
| fake/unverified repository DOI scan | PASS_NO_UNVERIFIED_REPOSITORY_DOI |
| Stage19 reopening scan | PASS_NO_STAGE19_REOPEN_COMMAND_OR_TRUE_FLAG |
| protected Figure/Results/Discussion scan | PASS_NO_FIGURE_RESULTS_DISCUSSION_MODIFICATION |
| historical_exact_version_not_recoverable | TRUE |
| biological rerun | FALSE |
