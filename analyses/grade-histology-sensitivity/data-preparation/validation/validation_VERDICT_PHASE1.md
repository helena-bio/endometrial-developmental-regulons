# SURVIVES

I cannot break it.

Binding verdict: **SURVIVES**. No FATAL or MAJOR objection remains. TASK B Phase 1 supports only the mechanical state `RESTRICT_TCGA_ONLY`; it does not support any molecular, attenuation, confounding, causal, target, or frozen-verdict claim.

## Independent reproduction

- Branch and HEAD are exact: `experiment/task-b-grade-histology-sensitivity`, `83503bad47b60193598b2b9ebe819c22c83e8ac1`.
- The four frozen protocol/access hashes, producer manifest hash, inventory hash, input-lock hash, harmonization hash, and decision hash all match the assigned seals.
- All 12 pre-existing allowlisted hashes match. The permitted schemas/columns were enforced. `gdc_case_mapping.json` is indeed a broad schema mapping rather than the described case crosswalk; it was denied and unused.
- The producer inventory has 68 data rows and all 68 resolve to the stated size and SHA-256. The mixed path convention is reproducible: `inputs/*` is relative to the task root; every other row is relative to `phase1_execution/`. This is a documented path-clarity issue, not a reproducibility failure.
- All 43 independently enumerated manifest references match. The canonical run comparison is 12/12 byte-identical.
- I reissued all six exact official acquisitions twice in fresh `env -i` processes under `strace`, validating wrapper and record keys before values. All returned HTTP 200, both reviewer rounds were byte-identical, and all raw and canonical SHA-256 values exactly equal the producer lock. Record counts are 529, 529, 83, 95, 123, and 159. There is no live API drift and no evidence of fabricated preserved bytes.
- From the locked bytes and allowlisted roster/subtype/linkage files, without importing producer analysis code, I reconstructed 737 patient-ledger rows. I compared 12,398 row/field values spanning the linkage ledger, immutable raw counts, missingness, raw and proposed subtype cross-tabs, harmonization map, and reported design diagnostics: 0 mismatches. Maximum numerical difference was `4.83081796787e-10`.

## Reproduced facts and decisions

### TCGA

- Exact patient/sample/aliquot/subtype linkage: n=507.
- Raw grade: G1=94, G2=113, G3=291, High Grade=9. Frozen binary grade: low=207, high=291, missing=9.
- Histology: endometrioid=382, non-endometrioid=124, missing=1.
- Base: PASS, n=507, p=4.
- Binary histology addition: PASS, n=506, p=5.
- Endometrioid-only plus grade: PASS, n=381, p=5.
- All-histology grade and all-histology grade-plus-histology are FORBIDDEN because TCGA `GRADE` does not document one common FIGO construct across retained histologies. Non-endometrioid tumours cannot be assigned grade 3 by assumption.
- Optional endometrioid grade-stratified models are FORBIDDEN. Low grade has n=205 and 13.6585% of rows above `3p/n`, exceeding the 5% ceiling; high grade passes alone at n=176. The rule requires at least two supported strata.

### CPTAC Discovery

- Frozen roster n=95; one case has two distinct primary-tumour samples and no exact frozen selector.
- cBio/PDC grade disagrees in 16 ordinal records; 7 cross the binary boundary. Histology conflicts in 12 records. No preferred-source repair is admissible.
- FIGO grade is missing in 12 and high grade has n=8. Histology after conflict is 83 endometrioid and 12 missing, hence one retained level.
- p53abn has 12/20 missing for both variables. The exact-link base is n=94 and fails leverage at 7.4468% above `3p/n`.
- Grade and histology are BLOCKED, independently of the support failures.

### CPTAC Confirmatory

- Frozen roster n=135; two cases have two aliquots and no exact frozen selector, leaving n=133 exact links.
- PDC `tumor_grade` is generic and has no locked FIGO semantics; grade is unavailable for all 135.
- Histology is one retained endometrioid level among the 133 exact links. Grade and histology are BLOCKED.

### Global state

The frozen definition of `RESTRICT_TCGA_ONLY` explicitly applies when TCGA supports at least one requested adjustment and every corresponding CPTAC route is BLOCKED or NOT_HARMONIZABLE. That is the observed state. Global `BLOCKED` would erase a specifically defined TCGA-only route; `PARTIAL_MODEL_HIERARCHY` would fail to preserve the CPTAC prohibition; `NOT_HARMONIZABLE` would ignore the clean, supported TCGA nodes. The producer global state is correct.

## Source semantics

- cBioPortal clinical attributes are study-defined/free-form; the TCGA `GRADE` label alone does not certify all-histology FIGO meaning.
- The locked CPTAC Discovery dictionary explicitly calls `HISTOLOGIC_GRADE_FIGO` Histologic FIGO Grade. Confirmatory PDC supplies only generic `tumor_grade`, so it cannot be mapped to FIGO without inference.
- ISGYP supports binary FIGO grades 1-2 versus grade 3 for endometrioid carcinoma and distinguishes serous, clear-cell, undifferentiated carcinoma, and carcinosarcoma from grade-3 endometrioid carcinoma.
- The ICD-O mapping is consistent with official morphology semantics: 8380/3 and 8382/3 are endometrioid; serous/clear-cell/mixed/undifferentiated codes remain non-endometrioid; 8140/3 adenocarcinoma NOS remains missing.
- The frozen map SHA is `c4252fd1d17a07212852a894631afd7a8df3b331ccf9aed314869b0e5512a362`. It was formed without molecular outcomes and cannot be changed after seeing Phase-2 results.

## Outcome-firewall attack

- Producer static scan: PASS, no findings.
- I independently parsed all five producer raw trace stages and the 1,325-row open audit: zero forbidden content opens and zero external analysis-network events. The two analysis runs used `bwrap --unshare-net`; acquisition used only the six declared endpoints.
- My two acquisition traces and offline reconstruction trace likewise show zero forbidden content opens and zero external analysis-network events.
- The initial HTTP 503 saved no body. The retry retained the identical endpoint and fields; the resulting bytes now reproduce live.
- The early superseded preflight report is not preserved, so its exact former bytes cannot be proven. The disclosed content was linkage identifiers only; the current report is sanitized. Linkage identifiers are permitted and are not grade/histology or molecular outcomes.
- NumPy 2.3.5 came from the host through a fresh `--system-site-packages` environment; no post-acquisition package download occurred. Independent diagnostics on the same NumPy version reproduce to sub-nanometric numerical tolerance.
- No expression, score, target result, coefficient, d/CI/p/q/direction, taxonomy, figure, manuscript, Task A scientific artifact, or TASK-030 result was accessed. No molecular model ran.

## Adversarial checklist C1-C8

- **C1 circularity: holds.** Clinical mappings/support rules were frozen before molecular access; no molecular outcome existed to tune a category, threshold, or model.
- **C2 leakage: holds.** Exact inputs and open traces show clinical/subtype/linkage-only access; acquisition and analysis were separated and analysis was offline.
- **C3 ecological fallacy: holds.** Phase 1 makes only roster/linkage/design feasibility statements, not per-patient biology or prognosis.
- **C4 reproducibility: holds.** Live source bytes, 68/68 inventory, 12/12 run pair, 43 manifest references, 737 ledger rows, and 12,398 compared fields independently regenerate.
- **C5 manuscript/code drift: holds within scope.** No manuscript or upstream science result is used. Protocol, producer, `src/**`, and `docs/**` current hashes/diffs remain unchanged.
- **C6 hidden terms: holds.** Phase 1 explicitly excludes the molecular covariate arrays and states that exact full-model diagnostics are untested and mandatory in Phase 2.
- **C7 statistics/design: holds for the binding state.** The fixed support, rank, condition, VIF/GVIF, Cramer's V, near-zero-variance, and leverage rules regenerate. Supported and forbidden TCGA nodes are correctly separated; CPTAC is blocked before adjusted inference.
- **C8 overreach: holds.** The decision states no molecular conclusion and cannot alter frozen TCGA/CPTAC verdicts.

## Minor unresolved reporting defect

**MINOR M1 (non-load-bearing): blocked CPTAC node labels reuse base diagnostics.** The producer reports Discovery `endometrioid_only` and `optional_grade_stratified` as n=94, p=4 base-design rows, and Confirmatory versions as n=133, p=4 base-design rows. Literal frozen-node reconstruction gives Discovery endometrioid-only plus grade n=83, p=5; optional low n=75 and high n=8 (both fail, high is rank-deficient); Confirmatory grade-dependent nodes have n=0. This is a machine-output labeling/coverage defect. It cannot change the verdict because exact linkage already BLOCKS both strata, Discovery high-grade support is inadequate, Confirmatory grade semantics are unavailable, and no CPTAC Phase-2 route exists. It must not be copied as if those placeholder rows were literal node diagnostics.

## Preservation and Phase-2 boundary

Producer start snapshots for external forbidden dirty worktrees were not captured. Current hashes plus zero-access traces cannot prove their historical byte equality; the producer disclosed this and did not claim false equality. This is a provenance limitation, not a material defect in the clinical-only result. The reviewer wrote only in `phase1_critic/`; protocol hashes, all 68 producer inventory entries, and `src/**`/`docs/**` diffs remain unchanged.

Phase 2 remains forbidden until this verdict is accepted and Sophia separately authorizes it. The only candidate TCGA nodes are: base, base plus binary histology, and endometrioid-only plus grade. Forbidden: all-histology grade, all-histology grade plus histology, optional grade-stratified models, all CPTAC adjusted inference, any rescue node, and any new target/category. Every exact full molecular design diagnostic and outcome-dependent influence diagnostic must be rerun in Phase 2.

This Phase-1 state is not evidence that GATA2 or SOX9 attenuates, that grade or histology confounds any molecular estimate, or that any frozen TCGA/CPTAC verdict changes.
