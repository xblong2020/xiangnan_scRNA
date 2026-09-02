# Repository provenance audit

Generated locally: 2026-09-02T08:36:24+08:00
Audit phase: initial

| Field | Value |
|---|---|
| REPOSITORY_EXISTS | TRUE |
| GIT_ROOT | C:/Users/Administrator/OneDrive/文档/湘南学院单细胞 |
| CURRENT_BRANCH | codex/module7-sctenifoldknk |
| HEAD_COMMIT | NULL |
| HEAD_STATE | UNBORN_OR_UNRESOLVED |
| REMOTE_EXISTS | FALSE |
| REMOTE_URL | NULL |
| TAG_EXISTS | FALSE |
| LATEST_TAG | NULL |
| WORKTREE_CLEAN | FALSE |
| ARCHIVAL_READY | FALSE |
| RELEASE_TAG | NULL |
| RELEASE_COMMIT | NULL |
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
| protected Figure/Results/Discussion scan | BASELINE_CREATED |
| historical_exact_version_not_recoverable | TRUE |
| biological rerun | FALSE |
