# Independent checksum and integrity report

Status: PASS.

- Analysis branch and HEAD: `experiment/task-a-perclass-c2-sensitivity` at `83503bad47b60193598b2b9ebe819c22c83e8ac1`.
- Frozen researcher bytes: charter SHA-256 `0dbd83c4dd69cca9fc70eecb9a8012bf6e91bfb9b51fcc0abd10f36c84d9851a`; spec SHA-256 `044bf4c5cdb80df5825e91cf5c6c9a29b17d2e75dac4f8d67f002f1655534340`.
- Producer manifest SHA-256: `29c8f67e4feefb4b2d7ab79e55834eda05e970bb8c2e41e084b86bd0f6ad9ecc`.
- Producer checksum inventory: 152/152 entries independently verified; zero missing or mismatched.
- Expanded pinned-input inventory: 111/111 entries independently rehashed; zero missing or mismatched. This covers the 21 directly pinned paths plus 90 nested seal entries (TASK-028 design 23, TASK-028 results 28, TASK-029 design 19, TASK-029 intermediate 7, TASK-029 results 13).
- Exact cohort guards independently reconstructed: TCGA primary 506 = 49/148/146/163; TCGA no-purity 507 = 49/148/147/163; CPTAC Discovery 95 = 7/25/43/20; CPTAC Confirmatory 135 = 6/47/66/16, in POLE/MMRd/NSMP/p53abn order.
- Original dirty worktree and TASK-030 worktree retained identical start/end HEAD, branch, status hash, tracked binary-diff hash, and staged binary-diff hash during the reviewer run. All relevant tracked `src/**` trees remain identical to HEAD and contain no untracked additions.
- The 44-row inherited DOCX/PDF audit still rehashes without a mismatch. In particular, the final TASK-030 DOCX is `01928833e3c595a57d4f42623a489b2f8966ef9d85443048d009c4a79d514d2c`, the final PDF is `0e31a0077b70fb75c3fd5e972ba3a2a634138e91290693c04e9a1d59c60a5e56`, and the original npj DOCX is `d22731647ff983976bec8a094a67c20945d81509a8d64a777216911741df81a2`.
- Definitive producer run1/run2: 27/27 files byte-identical.
- Reporting-schema correction: only six structured result files changed between preserved original run1 and definitive run1. The authoritative CPTAC meta table added exactly the six missing Discovery/Confirmatory coefficient, coefficient-SE, and residual-SD fields. The long-form representations gained a 15-column repeated meta provenance block; their d, bootstrap-SE, weight, and direction values already existed unchanged in the original meta table. Every common scientific field is identical; analytical report, counts, identities, BH families, taxonomy, presentation tables, and all figures are byte-identical. No model, seed, family, threshold, taxonomy, figure, verdict, or claim changed.

Machine evidence is in `independent_reproduction/INPUT_RECHECK.tsv`, `DEFINITIVE_RUN_PAIR.tsv`, `FULL_MACHINE_COMPARISON.tsv`, and `SUMMARY.json`.
