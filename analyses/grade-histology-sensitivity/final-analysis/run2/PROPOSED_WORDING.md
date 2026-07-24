# Proposed wording (proposal only; no manuscript edited)

## Results

In the post-hoc TCGA-UCEC explanatory sensitivity, the frozen C2 contrasts were refit on exactly matched rows before binary histology adjustment and, within endometrioid carcinoma, binary FIGO grade adjustment. Results are descriptive and are provided in MANUSCRIPT_READY_SENSITIVITY_TABLE.tsv; raw p values and diagnostic permutation p values carry no confirmatory or multiplicity credit.

## Methods

We used frozen TASK-028 signed regulon outcomes for six prespecified targets and OLS models containing the four-level molecular subtype factor, M4 proliferation, ESTIMATE-derived composition, and, in the primary model, CPE purity. C2 was 0.5*POLE + 0.5*MMRd - NSMP with p53abn retained at weight zero. Histology and endometrioid-only grade additions were treatment coded under the frozen Phase-1 map. We used 2,000 attempted patient bootstraps and 2,000 diagnostic subtype-label permutations with deterministic SHA-256-derived sub-seeds.

## Limitations

This post-hoc analysis is not a unique causal estimand. Grade and histology may be confounders, mediators, consequences of subtype biology, pathology proxies, or composition markers. CPE and ESTIMATE are estimated proxies and do not establish purity independence or complete cell-type adjustment. TCGA-only internal bootstrap is not external validation. The analysis supports no individual-patient biomarker, treatment response, target category, q value, verdict revision, or manuscript change.
