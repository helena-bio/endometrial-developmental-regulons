TASK-029 -- PATIENT-LEVEL SUBTYPE-JOIN CHECKPOINT (the last true execution gate)
================================================================================
Date 2026-07-14. METADATA-ONLY. Establishes per-PATIENT subtype-to-RNA joinability per
stratum for the CPTAC-UCEC external replication. STOP after; no freeze, no scoring, no
expression opened. Git HEAD unchanged: start == end == 83503bad47b60193598b2b9ebe819c22c83e8ac1.

CONTRASTS (frozen): C1 = p53abn vs {POLE,MMRd,NSMP}; C2 = {POLE,MMRd} vs NSMP.
F1 = six targets (LHX1,PAX8 on C1; GATA2,HOXA9,SOX9,WT1 on C2), BH-of-6 q<=0.05.
Evaluability = simulated BH-of-6 procedure-level power >= 0.80 (data-free, seed 20260714).
4-class map: POLE->POLE; MSI-H->MMRd; CNV_low/CNV_L->NSMP; CNV_high/CNV_H->p53abn.

================================================================================
HEADLINE VERDICT
================================================================================
BOTH strata are patient-level subtype-to-RNA JOINABLE from AUTHORITATIVE OPEN sources.
The Confirmatory per-case gap (previously PMC-gated) is now CLOSED: the Dou-2023 official
supplement mmc2.xlsx (Meta_table, column Genomic_subtype) is openly hosted on the Elsevier
CDN (ars.els-cdn.com, HTTP 200) and gives per-CASE 4-class labels. All classified cases in
BOTH strata join to a GDC CPTAC-3 primary-tumour RNA-seq file by case_submitter_id with
ZERO attrition. Pooled analytic set is unchanged from the frozen baseline; C2 stays
confirmatory-evaluable (BH-of-6 = 0.916 >= 0.80).

================================================================================
1-6. PER-STRATUM JOIN VERDICT + PER-SUBTYPE JOINED COUNTS
================================================================================

DISCOVERY (Dou 2020 Cell, PMID 32059776; PDC000125; cBioPortal ucec_cptac_2020)
  Subtype source: cBioPortal ucec_cptac_2020 GENOMICS_SUBTYPE, per-case, live API pull
    reproduced identically against prior cache (95 cases classified).
  4-class (raw): POLE 7 | MSI-H 25 | CNV_low 43 | CNV_high 20  (total 95).
  RNA join anchor: GDC CPTAC-3 uterine primary-tumour RNA-seq (STAR-Counts) case set.
  JOIN VERDICT: CONFIRMED. 95/95 classified Discovery cases join to a primary-tumour RNA
    file by case_submitter_id. NOT joinable = 0.
  PER-SUBTYPE JOINED: POLE 7 | MMRd 25 | NSMP 43 | p53abn 20  (total 95).
  Missing/conflicting: none classified-but-unjoinable; no label disagreement.

CONFIRMATORY (Dou 2023 Cancer Cell, PMID 37567170; PDC000439)
  Subtype source: Dou-2023 supplement mmc2.xlsx -> sheet 'Meta_table' -> column
    'Genomic_subtype' (per-case), OPEN via Elsevier CDN (the PMC copy is gated; the CDN
    media-object copy is not). Confirmed Confirmatory-cohort table: all 138 tumor cases
    flagged Discovery_study = No; 0 rows Discovery_study = Yes.
  Analytic-unit resolution: 140 Tumor rows = 138 distinct patients + 2 replicate aliquots
    (C3L-05571, C3N-02978), which the source itself marks Case_excluded = Yes (dedup rule k
    already applied at source). 138 distinct non-excluded tumor patients.
  4-class (raw, 138 patients): POLE 6 | MSI-H 47 | CNV_L 66 | CNV_H 16 | NA 3.
    -> exactly reproduces the frozen aggregate 6/47/66/16 + 3 unclassified.
  RNA join anchor: same GDC CPTAC-3 primary-tumour RNA-seq case set.
  JOIN VERDICT: CONFIRMED. 135/135 CLASSIFIED Confirmatory cases join to a primary-tumour
    RNA file by case_submitter_id. NOT joinable = 0.
  PER-SUBTYPE JOINED: POLE 6 | MMRd 47 | NSMP 66 | p53abn 16  (total 135).
  Missing/conflicting: 3 unclassified (Genomic_subtype = NA): C3L-00898, C3L-02802,
    C3N-02298 -- they HAVE RNA but no genomic 4-class label -> EXCLUDED from C1/C2 (not
    derived from outcome; not imputed). No label disagreement (no cross-source overlap).

INDEPENDENCE (per-case, re-confirmed at this granularity):
  Discovery case set (95) INTERSECT Confirmatory tumor case set (138) = 0 overlap.
  Donor-disjoint at case_submitter_id. CompRef contributes no Case_id row (technical only).

MULTI-SAMPLE / DEDUP (rule k = one patient = one primary-tumour analytic unit):
  On the RNA side, 199/240 CPTAC-3 uterine cases carry >1 primary-tumour RNA file -- these
  are replicate vial/portion aliquots of the SAME primary tumour (suffix pairs -0X/-1X, e.g.
  -01/-11, -02/-12; all sample_type = Primary Tumor). They collapse to ONE analytic unit per
  patient by rule k; NO patient is dropped (case count == patient count). Confirmatory
  within-study replicate aliquots (C3L-05571, C3N-02978) are already source-excluded.
  A per-patient primary aliquot selection (study-primary / lowest suffix) is applied at
  scoring time; it does not change any per-stratum patient count here.

UNAMBIGUOUS 4-CLASS MAPPING: genomic 4-class only, from the published Genomic_subtype /
  GENOMICS_SUBTYPE fields. NO histology fallback, NO surrogate, NO derived classifier.

================================================================================
7. RECOMPUTED CONTRAST SIZES + BH-of-6 POWER ON ACTUALLY-JOINABLE CASES
================================================================================
Joinable per-subtype (pooled): POLE 13 | MMRd 72 | NSMP 109 | p53abn 36  (total 230).

  stratum          C1 (p53abn vs rest)   C2 (POLE+MMRd vs NSMP)
  Discovery        20 vs 75              32 vs 43
  Confirmatory     16 vs 119            53 vs 66
  POOLED           36 vs 194            85 vs 109

Simulated BH-of-6 procedure-level power (seed 20260714, N_REP 200000, true_d 0.50; power
model byte-identical to simulate_power_f1_completecase.py; DATA-FREE -- only arm sizes):

  scenario                 C2 min BH power   C1 min BH power   C2 status (>= 0.80)
  JOINABLE_POOLED              0.916            0.766          EVALUABLE
  JOINABLE_DISCOVERY_only      0.449            0.400          COLLAPSES (single-stratum floor)
  JOINABLE_CONFIRM_only        0.695            0.403          COLLAPSES (single-stratum floor)

Cross-check: JOINABLE_POOLED bh_sim is IDENTICAL (byte-for-byte) to the prior frozen
BASELINE_no_purity -- because 100% of classified cases in BOTH strata are RNA-joinable, so
the join introduces ZERO attrition and the analytic sizes equal the frozen design.

DOES C2 STAY CONFIRMATORY-EVALUABLE? YES. Because BOTH strata are patient-level joinable,
the POOLED design holds: C2 min BH power = 0.916 >= 0.80. C2 does NOT collapse. (Had the
Confirmatory patient-level join failed -> Confirmatory 'not evaluable' -> Discovery-only,
C2 = 0.449 < 0.80 -> C2 evaluability WOULD have collapsed. It does not, because the join is
confirmed.) C1 remains 0.766 < 0.80 = confirmatory-underpowered / sensitivity-only, exactly
as frozen; C1 is not the evaluability gate.

================================================================================
8. CENTER / COMPREF STATUS
================================================================================
C3L/C3N ACCRUAL-FAMILY OVERLAP (on JOINABLE cases): Discovery C3L 51 / C3N 44;
  Confirmatory C3L 60 / C3N 75. Both families appear in BOTH strata -> accrual family is
  NOT fully confounded with the Discovery-vs-Confirmatory split.

SITE / BATCH RECOVERABILITY:
  - True tissue_source_site: NOT recoverable -- unpopulated in GDC (0/241 CPTAC-3 uterine
    cases) and absent from PDC open metadata. LIMITATION.
  - Recoverable per-case batch proxies WITH variation, NOT fully confounded with stratum:
      Confirmatory: mmc2 'Batch' b1/b2/b3/b4 (33/34/35/36, ~balanced), per case.
      Discovery: cBioPortal PROTEOMICS_TMT_BATCH / PROTEOMICS_TMT_PLEX / TUMOR_SITE, per case.
  VERDICT: a PRE-SPECIFIED center/batch SENSITIVITY IS POSSIBLE (accrual-family C3L/C3N and
    per-case processing-batch labels are recoverable and vary within stratum). CAVEAT: the
    recoverable batch labels are PROTEOMIC-TMT processing batches, not the RNA-seq batch and
    not a true recruiting site; they are a defensible batch/center PROXY, to be DISCLOSED as
    such. A true tissue_source_site correction remains a LIMITATION (no improvised site
    correction is invented). Center is NOT fully confounded with stratum -> no forced
    limitation-only outcome; the sensitivity is admissible if pre-specified.

COMPREF: CONFIRMED TECHNICAL REFERENCE ONLY. The CompRef pool contributes no Case_id /
  observational unit (0 Ref-like Case_id rows in the Confirmatory Meta_table; the only
  cross-study shared token 'Ref' is the internal TMT reference). EXCLUDED from donor counts
  and from biological uncertainty (uncertainty stays at the donor level).

================================================================================
EVALUABILITY DECISION
================================================================================
BOTH strata ENTER the confirmatory meta-analysis at the patient level:
  Discovery   = EVALUABLE (per-case subtype: cBioPortal GENOMICS_SUBTYPE; 95/95 RNA-joined).
  Confirmatory = EVALUABLE (per-case subtype: Dou-2023 mmc2 Meta_table Genomic_subtype, open
    Elsevier CDN; 135/135 classified RNA-joined; 3 NA excluded).
Neither stratum is 'not evaluable'. No label derived from outcome; no post-read classifier
built; no expression opened; no freeze authored. The pooled design is confirmed feasible with
C2 confirmatory-evaluable (0.916) and C1 sensitivity-only (0.766), matching the frozen design.

================================================================================
GAPS NOW CLOSED vs the prior CONFIRM-AT-ACQUISITION list
================================================================================
G1 (Dou-2023 per-CASE subtype roster): CLOSED. Resolved per-case from the open Elsevier-CDN
   supplement mmc2.xlsx (POLE 6/MMRd 47/NSMP 66/p53abn 16 + 3 NA), reproducing the aggregate.
G3 (3 Confirmatory unclassified): RESOLVED -- C3L-00898, C3L-02802, C3N-02298 (Genomic_subtype
   NA); excluded from C1/C2.
G4 (tissue_source_site): still a LIMITATION (unpopulated); batch proxy recoverable instead.
