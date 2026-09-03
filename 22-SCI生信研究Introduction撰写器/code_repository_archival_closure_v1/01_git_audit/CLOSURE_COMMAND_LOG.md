# Stage22 closure command log

Generated locally: 2026-09-03T09:34:41+08:00
Audit phase: final

- Read-only Git commands: rev-parse, status, branch, log, remote, tag, ls-files, cat-file, and for-each-ref.
- Read-only filesystem inventory: source counts, protected-scope inventory, large-file scan, sensitive-pattern scan, placeholder scan, DOI classification, and SHA-256 hashing.
- Local metadata builder: 22-SCI生信研究Introduction撰写器/code_repository_archival_closure_v1/01_git_audit/build_stage22_closure.py; it does not import or execute biological analysis modules.
- Git mutation scope: only the explicitly allowlisted release source/metadata files and Stage22 closure records.
- External push/upload: remote publication refs were verified for the configured remote.
- Local tag re-anchor: v1.0.0 was re-anchored before external publication from 302d29c7570b20c549cddc90580fc60a3a2ce4f9 to the corrected release commit; no remote tag existed.

No scVI/scanVI, CopyKAT, trajectory, SCENIC, CellOracle, scTenifoldKnk, external-validation, or figure-generation entry point was executed by this closure.

## Current external verification

- GitHub REST repository: public=true, default_branch=main.
- GitHub REST release v1.0.0: draft=false, prerelease=false, published=true.
- Zenodo API search for xiangnan_scRNA: total=0.
- Zenodo API search for xblong2020: total=0.
- No Zenodo token, OAuth session, or archive record was available to this workspace.
- Lightweight QA: python -B -m pytest tests; 274 passed in 28.42s at 2026-09-03T09:28:41+08:00.
