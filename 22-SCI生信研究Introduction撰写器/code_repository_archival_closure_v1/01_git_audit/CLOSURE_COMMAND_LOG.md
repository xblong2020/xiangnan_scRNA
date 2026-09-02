# Stage22 closure command log

Generated locally: 2026-09-02T18:19:21+08:00
Audit phase: final

- Read-only Git commands: rev-parse, status, branch, log, remote, tag, ls-files, cat-file, and for-each-ref.
- Read-only filesystem inventory: source counts, protected-scope inventory, large-file scan, sensitive-pattern scan, placeholder scan, DOI classification, and SHA-256 hashing.
- Local metadata builder: 22-SCI生信研究Introduction撰写器/code_repository_archival_closure_v1/01_git_audit/build_stage22_closure.py; it does not import or execute biological analysis modules.
- Git mutation scope: only the explicitly allowlisted release source/metadata files and Stage22 closure records.
- External push/upload: not attempted when no remote or archive authorization is configured.
- Local tag re-anchor: v1.0.0 was re-anchored before external publication from 302d29c7570b20c549cddc90580fc60a3a2ce4f9 to the corrected release commit; no remote tag existed.

No scVI/scanVI, CopyKAT, trajectory, SCENIC, CellOracle, scTenifoldKnk, external-validation, or figure-generation entry point was executed by this closure.

## GitHub publication events

- Target repository preflight: git ls-remote https://github.com/xblong2020/xiangnan_scRNA.git returned zero refs with exit code 0.
- Remote configured: origin = https://github.com/xblong2020/xiangnan_scRNA.git.
- Non-force branch publication: git push -u origin HEAD:main; remote main matched a7ee1edb8c139e9dca94a3e3702c750701e019e0.
- Non-force tag publication: git push origin v1.0.0; remote peeled tag target matched 1c0303049ba629b9b986cb9a7e088384f15f5a87.
- GitHub CLI audit: gh was not installed; no Release API call was attempted. Manual release instructions were generated.
- Zenodo DOI: no DOI was generated or inferred; PERMANENT_IDENTIFIER remains null.
