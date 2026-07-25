TASK-029 -- SUBTYPE PROVENANCE + DEFINITION CROSSWALK (two pre-freeze checks)
============================================================================
Date 2026-07-14 (retrieval verify 2026-07-14T12:35:24Z). METADATA + PUBLISHED-METHODS ONLY.
Subject: the Dou-2023 CPTAC-UCEC confirmatory subtype roster mmc2.xlsx -- the defining input
for subtype EXPOSURE in the TASK-029 external replication (Confirmatory stratum).

CONSTRAINT ATTESTATION (tumour-blind): only roster METADATA (Case_id + subtype label
vocabulary + the classifier ingredient LABELS) and published METHODS text were read. NO
expression value, NO regulon score, NO subtype effect, NO outcome, NO RNA matrix were opened
or downloaded. NO subtype label was derived or corrected from RNA / regulon / observed effect.
TASK-028 SECONDARY-B / ATAC not opened; sealed_v3 and the staged F34 not altered.
Git HEAD unchanged: start == end == 83503bad47b60193598b2b9ebe819c22c83e8ac1. No git write op;
no write to the analysis source tree. The acquisition record was written in an
isolated workspace.

============================================================================
CHECK 1 -- FROZEN PROVENANCE of mmc2.xlsx
============================================================================
Publisher / CDN URL (Elsevier CDN, Dou 2023 Cancer Cell PMID 37567170):
  https://ars.els-cdn.com/content/image/1-s2.0-S1535610823002477-mmc2.xlsx
  (article PII S1535610823002477; DOI 10.1016/j.ccell.2023.07.007; PMC10631452, whose PMC
   copy is behind the open-access download gate -- the CDN media-object copy returns HTTP 200
   and is the source used.)

Supplement file TITLE / VERSION:
  Table S1 (mmc2.xlsx) -- "Sample clinical and molecular meta information". The publisher does
  not stamp an explicit revision number on this xlsx; byte-identity is pinned by the SHA-256.

RETRIEVAL date+time (UTC): provenance verified 2026-07-14T12:35:24Z (date -u). The immutable
  copy was fetched by the TASK-029 join gate (curl -s -L, command in roster_provenance.json).

FILE SIZE: 277385 bytes.
SHA-256 (computed here with sha256sum on the local immutable copy):
  2ea92c4279918c3a6158b24ebfa1e0e36ffc2876d50fd2ed3f07265a389e0e31

WORKSHEET + COLUMNS (confirmed by reading xlsx metadata with openpyxl 3.1.5; header +
label vocabulary only, NOT any biological outcome):
  Worksheets present: README, Meta_table, Normal samples.
  Analysis worksheet:                 Meta_table  (213 columns; 190 non-empty data rows).
  (a) case/sample ID join column:     Case_id     (CPTAC C3L-/C3N- case_submitter_id).
  (b) subtype label column:           Genomic_subtype.
  Observed Genomic_subtype vocabulary: POLE / MSI-H / CNV_L / CNV_H / NA.
  Classifier-ingredient LABEL columns (metadata, read for the crosswalk only): POLE (Yes/No/NA),
  MSI_status (MSS/MSI-H/NA), CNV_status (CNV_L/CNV_H/NA); informational numeric CNV_ratio and
  MSIsensor_ratio exist but were NOT read as an outcome.

Local immutable copy path:
  external-inputs/dou2023_mmc2.xlsx

ROSTER LABEL COUNTS (metadata only; Group=Tumor, Case_excluded!=Yes; n=138 eligible tumors):
  POLE 6 | MSI-H 47 | CNV_L 66 | CNV_H 16 | NA(unclassified) 3.
  This reproduces the aggregate stated in the Dou-2023 main text verbatim: "Tumors were
  classified into four genomic subtypes: 6 POLE, 47 MSI-H, 66 CNV-L, 16 CNV-H tumors and 3
  unclassified because of missing WGS-based CNV data." The 3 NA cases (C3L-00898, C3L-02802,
  C3N-02298) are unclassified for a DATA-COMPLETENESS reason (missing WGS-based CNV), NOT a
  biological outcome; they carry RNA but no 4-class label and are EXCLUDED (not imputed, not
  derived from outcome).

RULE (explicit, binding on every future run):
  BEFORE any analysis, recompute and VERIFY sha256 == 2ea92c4279918c3a6158b24ebfa1e0e36ffc287
  6d50fd2ed3f07265a389e0e31 AND size == 277385 bytes on the local immutable copy. On MISMATCH:
  STOP. If the publisher replaces/updates the CDN file, the NEW version is NOT auto-accepted --
  re-pin (record the new URL/size/sha256) and re-verify worksheet + column names before any
  use. Silent substitution of a re-hosted file is forbidden (P4 provenance).

============================================================================
CHECK 2 -- SUBTYPE-DEFINITION CROSSWALK (estimand comparability, not classifier identity)
============================================================================
This is a crosswalk of DEFINITIONS between the Dou-2023 (and companion Dou-2020, same schema)
subtype scheme and the FROZEN TASK-028 subtype schema. It is NOT merely a confirmation that the
column is named "Genomic_subtype".

--- The FROZEN TASK-028 schema (target estimand) ---
  TCGA INTEGRATED GENOMIC 4-class (Kandoth et al., Nature 2013, PMID 23636398). Frozen subtype
  table = cBioPortal ucec_tcga_pan_can_atlas_2018 SUBTYPE. Raw vocabulary
  UCEC_POLE / MSI / CN_LOW / CN_HIGH, with the frozen equivalences CN_LOW = NSMP (the copy-
  number-low residual) and CN_HIGH = p53abn (copy-number-high / serous-like). Estimand =
  integration of somatic mutation burden + microsatellite instability + somatic copy-number
  alteration (SCNA). HIERARCHY: POLE-ultramutated first, then MSI-hypermutated, then
  copy-number-high (serous-like/p53abn), residual = copy-number-low (NSMP).
  Canonical Kandoth-2013 definitions (confirmed from the TCGA/pan-cancer literature): (1)
  'ultramutated' = POLE exonuclease-domain mutant, very high mutation rate; (2) 'hypermutated'
  = MSI (often via MLH1 promoter methylation), high mutation rate, few CNAs; (3) 'copy-number
  low' = MSS, low mutation rate, mostly grade 1-2 endometrioid; (4) 'copy-number high' =
  extensive SCNA, recurrent TP53 mutation, serous + ~25% grade-3 endometrioid.

--- The DOU-2023 scheme (as defined by the CPTAC-EC methods; Dou-2020 shares the schema) ---
  Both Dou-2020 (Cell, PMID 32059776, PMC7233456) and Dou-2023 (Cancer Cell, PMID 37567170,
  PMC10631452) state EXPLICITLY that tumors are "classified into the four genomic subtypes
  outlined in the TCGA EC landmark study" -- i.e. the SAME Kandoth-2013 integrated-genomic
  estimand -- into POLE (ultramutated), MSI-H (hypermutated), CNV-low (endometrioid-like), and
  CNV-high (serous-like). The 2023 unclassified cases are those with "missing WGS-based CNV
  data".
  APPLIED HIERARCHY, verified directly from the mmc2 Meta_table LABEL vocabulary (label-only
  cross-tab of POLE / MSI_status / CNV_status against Genomic_subtype; no outcome read):
    POLE=Yes -> Genomic_subtype=POLE regardless of MSI/CNV  (6 cases: 4 with CNV_L, 1 CNV_H,
        1 CNV_NA; all -> POLE).  => POLE assigned FIRST.
    then MSI_status=MSI-H (and POLE=No) -> MSI-H            (47: 45 CNV_L + 2 CNV_NA).
    then CNV_status=CNV_H (and POLE=No, MSS) -> CNV_H       (16).
    residual CNV_L / MSS / POLE=No -> CNV_L                 (66).
    CNV_status=NA with no POLE/MSI trigger -> NA            (3, unclassified).
  This applied order (POLE > MSI > CNV-high > CNV-low residual) is IDENTICAL to the frozen
  TASK-028 hierarchy. Ambiguous/unclassified handling: cases lacking the CNV determinant
  (missing WGS-based CNV) are left NA and excluded, not force-assigned.

--- Dou label VOCABULARY mapped onto the frozen POLE/MMRd/NSMP/p53abn axis ---
  (this is the mapping already frozen in the TASK-029 join gate; reproduced and justified here)
    Dou Genomic_subtype   Kandoth/TCGA raw     Frozen TASK-028 label
    -------------------   ------------------   ---------------------
    POLE              ->  UCEC_POLE        ->  POLE       (POLE exonuclease-domain ultramutated)
    MSI-H             ->  MSI              ->  MMRd        (MSI-hypermutated / mismatch-repair-deficient)
    CNV_L             ->  CN_LOW           ->  NSMP        (copy-number-low residual / no specific molecular profile)
    CNV_H             ->  CN_HIGH          ->  p53abn      (copy-number-high / serous-like / TP53-driven)
    NA                ->  (unclassified)   ->  EXCLUDED    (missing WGS-based CNV; not imputed)

--- Per-class definition comparison ---
  POLE:   Both = POLE exonuclease-domain hotspot mutation -> ultramutated. Same determinant
          (POLE EDM). Potential technical difference = the exact hotspot allele list applied;
          not resolvable from the open STAR-Methods text within the tumour-blind constraint.
          Effect is expected to be nil-to-negligible on class membership (POLE is the smallest,
          most sharply-defined arm; 6 cases here).
  MSI-H (=MMRd): Both = microsatellite-instability-high hypermutated. TCGA/Kandoth used a
          MSI/MMR determination; CPTAC-EC uses a sequencing-based MSI caller (the roster carries
          MSIsensor_ratio and MSI_status MSS/MSI-H). Same estimand (MSI hypermutation);
          DIFFERENT technical caller (MSIsensor-based vs the TCGA MSI/MMR-IHC/MLH1 call). The
          exact MSIsensor score cutoff is not resolvable from the open text.
  CNV-high (=p53abn): Both = copy-number-high / serous-like, TP53-driven. TCGA/Kandoth derived
          CN-high from an SCNA-cluster; Dou-2023 derives CNV_H from a WGS-based CNV metric
          (roster CNV_ratio / CNV_status). Same estimand (high SCNA burden / serous-like);
          DIFFERENT technical criterion (WGS CNV-burden threshold vs SCNA cluster). The exact
          CNV-high burden threshold is not resolvable from the open text. NOTE: the frozen
          TASK-028 CN_HIGH<->p53abn equivalence is itself already disclosed as TCGA-equivalent
          but NOT byte-identical at the p53abn boundary (a genomic CNV call vs a p53-IHC/TP53
          surrogate can disagree on the serous-like boundary) -- this crosswalk inherits that
          same, already-disclosed p53abn-boundary caveat.
  CNV-low (=NSMP): Both = the copy-number-low residual (MSS, non-POLE, non-CNV-high). Same
          residual definition; membership depends only on the other three calls, so it inherits
          their technical differences and adds none of its own.

--- Priority / ambiguity rule ---
  IDENTICAL priority order (POLE > MSI-H > CNV-high > CNV-low residual), verified from the label
  vocabulary of mmc2 itself. Ambiguous cases (missing the CNV determinant) are NA/excluded in
  Dou-2023, consistent with the frozen "no derived classifier, no histology fallback, no
  surrogate" rule of the join gate.

============================================================================
CONCORDANCE VERDICT
============================================================================
VERDICT: `concordant with disclosed implementation differences` (status 2 of 3).

REASON: The Dou-2023 Genomic_subtype scheme is the SAME four-class molecular estimand as the
frozen TASK-028 schema -- both are, by the authors' own statement, "the four genomic subtypes
outlined in the TCGA EC landmark study" (Kandoth 2013): integration of somatic mutation + MSI +
SCNA, with the IDENTICAL POLE > MSI > CNV-high > CNV-low-residual hierarchy (independently
verified here from the mmc2 label vocabulary). The vocabulary maps unambiguously and
biologically onto POLE/MMRd/NSMP/p53abn. It is therefore NOT status 3 (non-concordant): the
mapping is unambiguous and does not change any contrast's meaning.
It is NOT status 1 (definition-concordant, no caveat) because the TECHNICAL criteria differ
between the CPTAC-EC implementation and the TCGA implementation used to build the frozen
TASK-028 SUBTYPE table:
  - MSI-H call: MSIsensor-based (CPTAC) vs the TCGA MSI/MMR determination -- different caller,
    same MSI-hypermutation estimand.
  - CNV-high call: a WGS-based CNV-burden threshold (CPTAC) vs an SCNA-cluster (TCGA) -- different
    technical threshold, same high-SCNA/serous-like estimand; inherits the already-disclosed
    p53abn-boundary caveat.
  - POLE call: same POLE-EDM determinant; possible difference only in the exact hotspot allele
    list (not resolvable from open text; expected negligible on membership).
These are DISCLOSED IMPLEMENTATION DIFFERENCES on the SAME estimand, not an estimand change.

AFFECTED CASES (bounding, from labels only): the technical-criterion differences could only
matter at the MSI and CNV-high boundaries. In this roster those arms are MSI-H = 47 and
CNV_H = 16 (plus the 3 CNV-NA unclassified that are already excluded). POLE (6) and the CNV_L
residual (66) are not directly threshold-sensitive. No case can be RE-labelled here without
reading outcome/expression, which is forbidden; so the affected set is bounded (not enumerated
per-case) as the MSI-H and CNV_H boundary arms.

PROPOSED PRE-SPECIFIED SENSITIVITY (to be frozen by Vladimir, NOT authored here):
  Because the discordance is confined to the MSI (=MMRd) and CNV-high (=p53abn) boundary calls,
  pre-specify a labelled sensitivity that re-derives the C1 (p53abn vs rest) and C2
  ({POLE,MMRd} vs NSMP) contrasts under a boundary-perturbation: (i) hold the frozen Dou
  Genomic_subtype labels as PRIMARY; (ii) as a sensitivity, cross-check the MSI-H and CNV_H
  memberships against an independently-obtained TCGA-style call for the SAME CPTAC cases IF and
  only if such a call becomes available at acquisition WITHOUT reading the replication outcome,
  and report contrast-effect direction/stability across the two labelings. If no independent
  TCGA-style call is obtainable tumour-blind, report the implementation difference as a DISCLOSED
  LIMITATION (comparable-estimand, not byte-identical classifier) rather than inventing a
  correction. No label is to be re-derived from the observed replication effect (that would be
  P1/P2 circularity).

The confirmatory execution is NOT blocked by this verdict (status 2, not 3): the estimand is the
same and the contrasts keep their meaning; the technical differences are documented and carried
as a disclosed implementation caveat + a pre-specifiable sensitivity, per the anti-fabrication
estimand-change rule.

============================================================================
SOURCES
============================================================================
- Dou Y, et al. Proteogenomic insights suggest druggable pathways in endometrial carcinoma.
  Cancer Cell 2023;41(9):1586-1605. PMID 37567170; PMC10631452; DOI 10.1016/j.ccell.2023.07.007.
  Supplement mmc2.xlsx via Elsevier CDN
  https://ars.els-cdn.com/content/image/1-s2.0-S1535610823002477-mmc2.xlsx
  (main-text quote on the 6/47/66/16 + 3-unclassified counts and "missing WGS-based CNV data",
   and "four genomic subtypes outlined in the TCGA EC landmark study").
- Dou Y, et al. Proteogenomic characterization of endometrial carcinoma. Cell 2020;180:729-748.
  PMID 32059776; PMC7233456; DOI 10.1016/j.cell.2020.01.026 (companion, same subtype schema;
  7 POLE / 25 MSI / 43 CNV-low / 20 CNV-high; classification per the TCGA EC landmark study).
- Kandoth C, et al. Integrated genomic characterization of endometrial carcinoma. Nature
  2013;497:67-73. PMID 23636398 (the four-class integrated-genomic definition + hierarchy).
- Frozen TASK-028 subtype table: cBioPortal ucec_tcga_pan_can_atlas_2018 SUBTYPE
  (UCEC_POLE/MSI/CN_LOW/CN_HIGH); frozen equivalences CN_LOW=NSMP, CN_HIGH=p53abn.
- TASK-029 join gate (this work dir): PATIENT_LEVEL_SUBTYPE_GATE.md,
  MANIFEST_patient_level_subtype_gate.md (the pinned URL/sha256 and the frozen 4-class map).

NO expression opened; NO label derived from outcome; NO freeze authored. STOP.
Git HEAD start == end == 83503bad47b60193598b2b9ebe819c22c83e8ac1.
