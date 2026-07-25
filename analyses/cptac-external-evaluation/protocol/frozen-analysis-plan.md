TASK-029 -- FROZEN EXTERNAL-REPLICATION DESIGN (authoritative; approved by Vladimir 2026-07-14)
================================================================================
This is the SEALED frozen design for the CPTAC-UCEC external replication of the sealed TASK-028
finding (F34). It is FEASIBILITY-validated, NOT an achieved replication. Final replicated/not-
replicated is assigned ONLY after frozen execution + independent verification. No auto-launch
before the independent critic confirms this sealed design. F34 stays sealed/staged; git is not a
blocker. Sub-artifacts are SHA-pinned in SEAL_MANIFEST.sha256; on any conflict, this file governs
the design and FINAL_FREEZE_RULES_bm.md governs the b-m rule wording.

--------------------------------------------------------------------------------
0. STATUS (pre-execution; feasibility)
--------------------------------------------------------------------------------
  C2 (GATA2, HOXA9, SOX9, WT1): confirmatory-evaluable (simulated BH-of-6 = 0.916 on joinable pooled).
  C1 (LHX1, PAX8): underpowered for confirmatory replication (BH 0.766), sensitivity-only.
  M1 negative: equivalence-underpowered, sensitivity-only.
  Headline: "CPTAC permits a rigorous confirmatory test of C2, but NOT of C1 or M1."

--------------------------------------------------------------------------------
1. COHORT, INDEPENDENCE, JOIN (all verified pre-freeze, tumour-blind)
--------------------------------------------------------------------------------
  Cohort = pooled CPTAC-UCEC: Discovery (Dou 2020 Cell) + Confirmatory (Dou 2023 Cancer Cell),
  two DONOR-INDEPENDENT strata. Independence CLEAN: 0 overlap at case/sample/aliquot (CPTAC_
  INDEPENDENCE.md); dedup rule k (one patient = one primary-tumour analytic unit; keep study-
  primary aliquot for the 3 replicate-aliquot patients). CompRef 'Ref' pool = technical reference
  only, excluded from donor counts / subtype contrasts / donor-level uncertainty.
  Patient-level subtype-to-RNA join CONFIRMED both strata (PATIENT_LEVEL_SUBTYPE_GATE.md):
  Discovery 95/95, Confirmatory 135/135, zero attrition, 3 Confirmatory NA excluded (missing WGS-
  CNV, not outcome-derived). Joined per-subtype: Disc 7/25/43/20; Conf 6/47/66/16. Pooled joinable
  contrast sizes: C1 p53abn 36 vs 194; C2 POLE+MMRd 85 vs NSMP 109.

--------------------------------------------------------------------------------
2. SUBTYPE ROSTER -- FROZEN PROVENANCE (defining exposure input)
--------------------------------------------------------------------------------
  Confirmatory roster = Dou-2023 supplement mmc2.xlsx (open Elsevier CDN
  https://ars.els-cdn.com/content/image/1-s2.0-S1535610823002477-mmc2.xlsx), Table S1, retrieved
  2026-07-14T12:35:24Z, 277385 bytes, SHA-256 2ea92c4279918c3a6158b24ebfa1e0e36ffc2876d50fd2ed3f07265a389e0e31,
  worksheet Meta_table, join column Case_id, subtype column Genomic_subtype (POLE/MSI-H/CNV_L/CNV_H/NA).
  Immutable local copy: sources_task029gate/dou2023_mmc2.xlsx. RULE: every run VERIFIES sha256 AND
  size==277385 before analysis; a publisher-replaced file is NOT auto-accepted (new explicit
  decision + re-pin required). Discovery roster = cBioPortal ucec_cptac_2020 GENOMICS_SUBTYPE.

--------------------------------------------------------------------------------
3. SUBTYPE-DEFINITION CONCORDANCE (crosswalk) + MMRd/C2 SENSITIVITY (frozen interpretation)
--------------------------------------------------------------------------------
  Verdict = `concordant with disclosed implementation differences` (status 2 of 3;
  SUBTYPE_PROVENANCE_AND_CROSSWALK.md): SAME TCGA/Kandoth-2013 integrated-genomic four-class
  estimand + identical hierarchy (POLE>MSI>CNV-high>CNV-low residual); map POLE->POLE, MSI-H->MMRd,
  CNV_L->NSMP, CNV_H->p53abn. Technical implementation differs (MSIsensor vs TCGA MSI/MMR; WGS
  CNV-burden vs SCNA cluster; POLE hotspot list). DISCLOSED LIMITATION: the exact STAR-Methods
  cutoffs (MSIsensor, CNV-burden, POLE alleles) are not recoverable from open sources -> the
  classifiers are NOT presented as algorithmically identical.
  MMRd/C2 SENSITIVITY -- FROZEN INTERPRETATION (Vladimir 2026-07-14):
    - Dou-2023 patient-level labels are PRIMARY and define the confirmatory C2 analysis.
    - An alternative TCGA-style MMRd/p53abn call is used ONLY IF: available BEFORE the biological
      read; INDEPENDENT of RNA expression / regulon scores / observed effects; with pre-documented
      provenance + mapping; requiring NO outcome-informed relabeling.
    - The sensitivity: cannot turn a non-passing primary C2 into `replicated`; does not replace the
      primary subtype schema; does not automatically overturn a passed primary result; assesses
      robustness to classifier implementation only.
    - If primary C2 passes but the sensitivity direction is unstable, the final language is
      EXACTLY: "C2 replicated under the prespecified Dou-2023 genomic subtype implementation, but
      not robust to the available alternative MMRd classification."
    - If a reliable alternate call is NOT available, the sensitivity is NOT improvised; the
      technical classifier difference remains a disclosed limitation.

--------------------------------------------------------------------------------
4. FAMILIES, MODEL, CONTRASTS, MULTIPLICITY, TAXONOMY (from FINAL_FREEZE_RULES_bm.md + v3)
--------------------------------------------------------------------------------
  F1 = six positive targets (C1 LHX1/PAX8; C2 GATA2/HOXA9/SOX9/WT1), BH FDR q<=0.05 across the six.
  F2 = M1 negative (separate equivalence; M1 NOT in F1 multiplicity).
  Scoring = the sealed TASK-028 modules verbatim (M1 ssGSEA; M3 signed VIPER/aREA over the 761
  CollecTRI-signed edges), NO retuning. PRIMARY MODEL = VARIANT 2 FROZEN NO-PURITY model: subtype
  + M4 proliferation (covariate-only) + ESTIMATE composition; SAME contrasts C1=[-1,-1,-1,+3],
  C2=[+1,+1,-2,0], C3=[+1,-1,0,0] (equal-within-side, full-model); SAME full-model structure.
  ESTIMATE-purity = disclosed sensitivity only (expression-derived, non-identical, possibly
  collinear); its frozen sub-rules (missingness/complete-case, subtype counts, collinearity
  diagnostics, model-instability rule, cannot-override-primary) apply. WORDING: "The primary model
  avoids adjustment with an expression-derived, non-identical purity surrogate, while residual
  confounding by tumour purity remains a disclosed limitation." Primary CPTAC analysis = a
  MODEL-ADAPTED external replication, NOT an identical-covariate replication of TASK-028.
  META = inverse-variance FIXED-EFFECT primary over the two strata; stratum-specific direction+CI
  mandatory; OPPOSITE-DIRECTION VETO; random-effects sensitivity-only; NO raw-expression pooling.
  b `replicated` (per F1 target) requires ALL FIVE: pre-specified direction; |d|>=0.5 + bootstrap
  CI excl 0; BH q<=0.05; stratum-consistency; no critical QC/model-fit failure. Power = evaluability
  BEFORE read, NOT the criterion. TAXONOMY: not evaluable / underpowered / evaluable-not-replicated
  / replicated (no "FAIL" for underpowered). c-k carried (v3 sec 7): labels=genomic 4-class; all-
  zero-across-cohort gene rule; gene-ID map = GDC GENCODE v36 / 60660 identity (verify at
  acquisition); seed 20260714; etc.

--------------------------------------------------------------------------------
5. M1 (F2), CENTER, COMPREF -- disclosed limitations
--------------------------------------------------------------------------------
  M1: NOT a confirmatory external replication target. RECORD verbatim: "External equivalence
  confirmation of the M1 negative finding is not feasible in the available independent cohort."
  Sensitivity-only; margin 0.30 NOT widened; non-significance != confirmed negative; no post-hoc
  PASS. CENTER: true tissue_source_site unrecoverable (limitation); C3L/C3N accrual families in
  BOTH strata (not fully confounded); a TMT-batch center sensitivity is a DISCLOSED PROXY only
  (not a true recruiting-site correction); if no reliable site variable, it stays a disclosed
  limitation with NO improvised statistical correction. COMPREF: technical reference only.

--------------------------------------------------------------------------------
6. GATE ORDER (no auto-launch before the critic confirms this sealed design)
--------------------------------------------------------------------------------
  seal final freeze [THIS] -> independent critic verification of the sealed frozen design ->
  RNA acquisition + frozen execution -> independent result verification. Final replicated/not-
  replicated only after execution + independent verification. No SECONDARY-B, no ATAC, no scoring
  before the critic confirms + the RNA-acquisition gate is opened. SEED master 20260714.
