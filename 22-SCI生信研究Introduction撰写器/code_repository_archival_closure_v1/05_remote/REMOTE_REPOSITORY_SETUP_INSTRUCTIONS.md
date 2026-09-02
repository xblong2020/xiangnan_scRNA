# Remote repository setup instructions

Current status: REMOTE_REPOSITORY_STATUS = NOT_CONFIGURED

Minimal one-time external action:

1. Create or select a public GitHub/GitLab repository owned by the project team. Do not upload raw/restricted data, patient-level files, large sequencing objects, caches, environments, credentials, or unapproved derived outputs.
2. From the project root, add the verified remote URL: git remote add origin <PUBLIC_REPOSITORY_URL>.
3. Push the current release branch: git push -u origin codex/module7-sctenifoldknk.
4. Push the immutable release tag: git push origin refs/tags/v1.0.0.
5. Verify the public repository URL, branch mapping, tag target commit, and release visibility. Then enable Zenodo or another approved archive for that repository.

The placeholder URL above is an instruction placeholder, not a repository claim. Replace it only with the URL returned by the repository provider.
