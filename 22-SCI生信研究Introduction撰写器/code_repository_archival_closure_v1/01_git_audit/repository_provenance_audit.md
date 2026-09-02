# Repository provenance audit

Generated locally: 2026-09-02T18:19:21+08:00
Audit phase: final

| Field | Value |
|---|---|
| REPOSITORY_EXISTS | TRUE |
| GIT_ROOT | C:/Users/Administrator/OneDrive/文档/湘南学院单细胞 |
| CURRENT_BRANCH | codex/module7-sctenifoldknk |
| HEAD_COMMIT | e1f8122c26847a32629df37da8c0de21ea11a657 |
| HEAD_STATE | COMMIT |
| REMOTE_EXISTS | FALSE |
| REMOTE_URL | NULL |
| TAG_EXISTS | TRUE |
| LATEST_TAG | v1.0.0 |
| WORKTREE_CLEAN | TRUE |
| ARCHIVAL_READY | FALSE |
| RELEASE_TAG | v1.0.0 |
| RELEASE_COMMIT | 1c0303049ba629b9b986cb9a7e088384f15f5a87 |
| REMOTE_REPOSITORY_STATUS | NOT_CONFIGURED |
| REMOTE_PUSH_STATUS | NOT_APPLICABLE_NO_REMOTE |

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
