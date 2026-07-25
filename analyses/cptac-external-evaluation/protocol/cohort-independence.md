TASK-029 (A): CPTAC-UCEC DISCOVERY-vs-CONFIRMATORY INDEPENDENCE + REAL PER-STRATUM COUNTS
======================================================================================
Recon date 2026-07-14. TUMOUR-BLIND: only metadata (case/sample/aliquot IDs, subtype
labels) read from OPEN sources. NO expression value read. NO matrix downloaded.

Git HEAD unchanged during this work: 83503bad47b60193598b2b9ebe819c22c83e8ac1.

------------------------------------------------------------------------------
SOURCES (all open metadata)
------------------------------------------------------------------------------
- GDC CPTAC-3 case/biospecimen API (https://api.gdc.cancer.gov): 241 uterine
  adenocarcinoma cases in project CPTAC-3.
- cBioPortal study ucec_cptac_2020 (Discovery): 95 tumour samples with the
  curated GENOMICS_SUBTYPE 4-class attribute. (Dou 2020 Cell, PMID 32059776,
  PMC7233456.)
- PDC (Proteomic Data Commons) GraphQL biospecimenPerStudy:
    Discovery   = PDC000125 "CPTAC UCEC Discovery Study - Proteome"    (154 aliquots).
    Confirmatory= PDC000439 "CPTAC UCEC Confirmatory Study - Proteome" (161 aliquots).
- Dou 2023 Cancer Cell (PMID 37567170, PMC10631452) Table S1 (paper text) for the
  Confirmatory aggregate subtype breakdown.

------------------------------------------------------------------------------
1. ROSTERS
------------------------------------------------------------------------------
Discovery (Dou 2020, PDC000125): 103 primary-tumour patients profiled by proteome;
  95 of them carry a genomic 4-class subtype in the Cell-2020 analysis (cBioPortal).
Confirmatory (Dou 2023, PDC000439): 138 primary-tumour patients (140 tumour
  aliquots; 2 patients have a replicate aliquot of the same tumour).
Full case/sample/aliquot ID lists are in the PDC JSON dumps under sources/ and the
  per-subtype Discovery roster in sources/discovery_roster_by_subtype.json.

------------------------------------------------------------------------------
2. OVERLAP CHECK  --  VERDICT: CLEAN (zero donor-level overlap)
------------------------------------------------------------------------------
Computed on PDC biospecimen IDs, PRIMARY TUMOUR samples only:
    case_submitter_id  overlap = 0
    sample_submitter_id overlap = 0
    aliquot_submitter_id overlap = 0
The ONLY ID shared across the two studies (over all sample types) is the literal
string "Ref" = the CompRef internal reference pool, NOT a patient.
No patient appears in both cohorts. No Discovery sample is reused/renamed in the
Confirmatory study. Within each study, no patient appears more than once as a
distinct tumour (2 Discovery + 2 Confirmatory cases carry a replicate ALIQUOT of
the SAME tumour: Disc C3N-01825; Conf C3N-02978, C3L-05571 -> aggregate to one
value per patient at scoring; a per-patient dedup, not a cross-cohort issue).
Cross-validation: all 95 cBioPortal Discovery case IDs are present in GDC CPTAC-3
and NONE appears in the Confirmatory tumour roster.

------------------------------------------------------------------------------
3. EXACT PER-STRATUM SUBTYPE COUNTS (4-class genomic)
------------------------------------------------------------------------------
DISCOVERY (Dou 2020) -- RESOLVED from cBioPortal GENOMICS_SUBTYPE
  (mapping: POLE->POLE, MSI-H->MMRd, CNV_low->NSMP, CNV_high->p53abn):
    POLE=7  MMRd=25  NSMP=43  p53abn=20   (total classified = 95)

CONFIRMATORY (Dou 2023) -- aggregate RESOLVED from paper text (Table S1),
  per-CASE labels = CONFIRM-AT-ACQUISITION (see gaps):
    POLE=6  MMRd=47  NSMP=66  p53abn=16   (+3 unclassified; total = 138)

POOLED (Discovery + Confirmatory, subtype-classified tumours):
    POLE=13 MMRd=72 NSMP=109 p53abn=36   (total = 230)

REAL CONTRAST GROUP SIZES
  C1 = p53abn vs pooled {POLE,MMRd,NSMP}
       Discovery      20 vs 75
       Confirmatory   16 vs 119
       POOLED         36 vs 194
  C2 = pooled {POLE,MMRd} vs NSMP
       Discovery      32 vs 43
       Confirmatory   53 vs 66
       POOLED         85 vs 109

NOTE / DEVIATION vs spec's EXPECTED sizes: the spec anticipated C1 ~50 vs ~183 and
C2 ~78 vs ~95. The RESOLVED real pooled sizes are C1 36 vs 194 and C2 85 vs 109.
The p53abn arm (36) is notably SMALLER than the ~50 anticipated. Power (part B) is
computed on the REAL resolved sizes; the spec's expected sizes are ALSO run as a
sensitivity.

------------------------------------------------------------------------------
4. DEDUPLICATION RULE (frozen proposal)
------------------------------------------------------------------------------
No cross-cohort duplicate exists, so no cross-cohort dedup is needed. For the
within-study replicate ALIQUOTS (same tumour, 3 patients total), freeze:
  "Where a patient has >1 tumour aliquot in a study, keep the single aliquot the
   source study designated primary (lowest aliquot suffix / study-primary flag);
   collapse replicate aliquots to one score per patient BEFORE contrast tests."
Power is already computed on PATIENT counts (138, not 140), so no recompute needed.

------------------------------------------------------------------------------
5. SHARED RECRUITING CENTER (disclosure)
------------------------------------------------------------------------------
Both cohorts draw from the same two CPTAC accrual-group families (case-ID prefixes
C3L and C3N): Discovery tumour C3L=58/C3N=45; Confirmatory tumour C3L=62/C3N=76.
This is SHARED recruiting infrastructure -> possible common-site batch effect,
DISCLOSED. It does NOT break donor-level independence (zero patient overlap).
The true tissue_source_site field is not populated in GDC/PDC open metadata ->
site-level confounding is CONFIRM-AT-ACQUISITION.

------------------------------------------------------------------------------
CONFIRM-AT-ACQUISITION GAPS (do not fabricate)
------------------------------------------------------------------------------
G1. Dou 2023 Table S1 per-CASE genomic subtype labels: aggregate resolved
    (6/47/66/16 + 3 unclassified); per-case assignment behind the PMC proof-of-work
    download gate at recon time. Confirm the per-case roster on acquisition.
G2. 8 Discovery PDC tumour cases beyond the 95 subtyped in Cell-2020
    (C3L-00084, C3L-00157, C3L-00356, C3L-00938, C3L-01247, C3L-01253, C3L-01284,
    C3N-01001): confirm they are QC/subtype-unassignable exclusions.
G3. The 3 Confirmatory "unclassified" tumours: confirm they are excluded from C1/C2.
G4. tissue_source_site / recruiting-center identity: not in open metadata.

VERDICT: INDEPENDENCE = CLEAN. Discovery and Confirmatory are donor-disjoint at
case/sample/aliquot level. Real pooled contrast sizes C1 36/194, C2 85/109.
