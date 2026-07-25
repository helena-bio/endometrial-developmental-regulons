TASK-028 -- B2 TUMOUR MODEL SPECIFICATION (FROZEN CANDIDATE for seal)
================================================================================
Encodes Vladimir's pre-seal corrections (2026-07-13). FROZEN BEFORE any tumour read.
No tumour value informs any threshold. On seal this file is read-only; its SHA-256
is recorded in SEAL_MANIFEST.sha256. Authoritative power justification:
researcher_v2/JUSTIFICATION_B2_H.md; a-priori power: experimenter_final/POWER_CONTRASTS.md.

--------------------------------------------------------------------------------
B2.1 INPUT / NORMALIZATION  [settled at Freeze A]
--------------------------------------------------------------------------------
Cohort = the frozen 507 primary tumours, one-per-patient (SELECTION_RULE b92d7478,
subtype table 17a899a4). No sample added/dropped. Values = augmented_star_gene_counts
"tpm_unstranded" -> log2(TPM+1). Single GDC STAR / GENCODE-v36 pipeline; NO TOIL/RSEM
mixing; NO batch correction as a primary substitute. Feature space = 60660 GENCODE v36
genes; module genes mapped via the frozen B1 mapping rule.

--------------------------------------------------------------------------------
B2.2 LOW-EXPRESSION RULE (frozen; Vladimir point 7)
--------------------------------------------------------------------------------
NO subtype-aware, variance, quantile, or top-variable filtering -- EVER. Module
membership is frozen (B1). After the separately-authorized tumour read, the ONLY
admissible mechanical removal is a gene that is TECHNICALLY UNDETECTABLE across the
ENTIRE cohort (TPM == 0 in ALL 507 samples). Every such removal is logged (gene,
module, "all-zero-507"); there is NO gene replacement and NO module re-optimization.
This is the sole post-read change permitted to any frozen module.

--------------------------------------------------------------------------------
B2.3 SCORING (Vladimir point 3)
--------------------------------------------------------------------------------
M1 (and compact-M1 sensitivity): per-sample ssGSEA on log2(TPM+1) (params pinned).
M3 regulons -- PRIMARY = SIGNED VIPER/aREA over the CollecTRI-signed regulons, under the
STRICTLY CONSENSUS-SIGN-BASED policy (v2 amendment, DEC-M3-DUAL-EDGE; edge ledger
B1/m3_primary_edge_ledger.tsv, reconciliation m3_edge_reconciliation.md):
  - PRIMARY includes EXACTLY the 761 CONSENSUS-RESOLVED admitted edges, each ONCE with its
    CollecTRI consensus sign (+1 activation / -1 repression). This includes the 59 both-
    flagged edges that resolved to a single consensus sign.
  - ALL 41 NO-CONSENSUS admitted edges are EXCLUDED from primary (20 both-flagged/no-
    consensus + 21 single-flag/no-consensus). An unresolved edge is NOT split, NOT given
    a raw-flag sign, NOT given any inferred/manual sign. Consensus-only; no exceptions.
  - per-edge likelihood weight = CollecTRI consensus confidence (1.0 consensus / 0.5
    non-consensus).
  - ORPHAN TARGETS (the 21 unique non-CORE targets whose only incoming edges are no-
    consensus, listed in B1/orphan_target_ledger.tsv) receive NO zero / arbitrary / raw-
    flag contribution: they simply do not enter primary signed scoring, recorded by name
    with reason "no consensus-signed primary edge". (Source target membership 583; primary
    targets with >=1 included edge 562; orphan targets 21. Membership files 603/577 are the
    frozen provenance universe, byte-identical to v1 -- unchanged.)
TWO PRE-FROZEN SENSITIVITY analyses (reported SEPARATELY from primary; NOT chosen post hoc):
  - S1 flag-split: ALL 79 both-flagged edges (the 59 consensus-resolved AND the 20
    unresolved) are represented as TWO separate contributions (activation +1, repression
    -1), each carrying its consensus-confidence weight; non-both-flagged edges enter once.
  - S2 conservative exclusion: ALL 79 both-flagged edges EXCLUDED; single-flag consensus-
    resolved edges retained by consensus sign; single-flag no-consensus edges remain
    excluded as in primary.
Signed-ssGSEA = a PRE-FROZEN SENSITIVITY method. The scoring method is fixed a priori;
it is NOT chosen by the result.
Scoring-method sensitivity: mean-z / singscore concordance (reported, not primary).

--------------------------------------------------------------------------------
B2.4 PRIMARY MODEL + CONTRAST MATRIX (Vladimir points 1 + full-model weights)
--------------------------------------------------------------------------------
One covariate-adjusted model per scorable module m (M1_analysis_ready[236]; each M3
CORE_TF signed regulon):
    score_{m,i} = b0 + SUBTYPE_i(4 levels, 3 df) + b_prolif*M4_i + b_purity*purity_i
                  + b_comp*composition_i + eps_i
PRIMARY (purity-adjusted) is fit on the CPE-COMPLETE-CASE n=506 (B2.6; TCGA-BS-A0TG
excluded ONLY from CPE-including models; subtype counts POLE49/MMRd148/NSMP146/p53abn163).
The frozen NO-PURITY model (drop b_purity) is fit on the FULL cohort n=507 as a
pre-labelled robustness sensitivity (B2.6/B2.10), NOT a substitute for primary.
Obtain the 4 COVARIATE-ADJUSTED subtype marginal means mu = (mu_POLE, mu_MMRd,
mu_NSMP, mu_p53abn) (EMMs at reference covariate values).

CONTRASTS are applied to mu via the FULL-MODEL contrast matrix with FIXED weights --
NOT nested-subset refits. Group sizes do NOT change the biologically-set weights: each
side of a contrast pools its subtypes by EQUAL SUBTYPE WEIGHT (unweighted mean of the
subtype EMMs), never n-weighted. Integer weight vectors on [POLE,MMRd,NSMP,p53abn]:
    C1 = [-1,-1,-1,+3]  p53abn vs pooled non-p53abn
    C2 = [+1,+1,-2, 0]  (POLE+MMRd) vs NSMP, within non-p53abn (p53abn weight 0)
    C3 = [+1,-1, 0, 0]  POLE vs MMRd
Mutually orthogonal; span the full 3 df. For the standardized effect, weights are
normalized to a difference-of-means form (positive side and negative side each an
EQUAL-weight mean of their subtype EMMs):
    d_C1 = (mu_p53abn - mean(mu_POLE,mu_MMRd,mu_NSMP)) / sigma_resid
    d_C2 = (mean(mu_POLE,mu_MMRd) - mu_NSMP)          / sigma_resid
    d_C3 = (mu_POLE - mu_MMRd)                        / sigma_resid
C2/C3 are computed IN THE FULL MODEL (p53abn present, weight 0), NOT by refitting on a
non-p53abn subset -- so the subset composition cannot re-weight the contrast. Report
each contrast estimate + PATIENT-level bootstrap CI + d. The four individual
subtype-vs-rest comparisons are retained DESCRIPTIVE-ONLY.

--------------------------------------------------------------------------------
B2.5 OMNIBUS = FORMAL HIERARCHICAL GATE (Vladimir point 1)
--------------------------------------------------------------------------------
Per primary module, a 3-df OMNIBUS subtype test (adjusted Wald/LRT on the subtype
factor; Kruskal-Wallis companion) is tested in family F1.
  - C1-C3 are ALWAYS evaluated and reported with effect size + CI (transparency).
  - A CONFIRMATORY subtype-category claim for a module is admissible ONLY if that
    module's omnibus passes the F1 threshold (BH q<=0.05).
  - If the omnibus is NEGATIVE, C1-C3 for that module remain DESCRIPTIVE/EXPLORATORY,
    regardless of any nominal contrast p-value.
The omnibus is a formal gate (gatekeeping / fixed-sequence), NOT a ban on evaluation.

--------------------------------------------------------------------------------
B2.6 COVARIATES (Vladimir point 4)
--------------------------------------------------------------------------------
PROLIFERATION: M4 score, partialed. A call that vanishes on M4 adjustment is
  proliferation-driven, not credited (R-10). M4 is COVARIATE-ONLY, in no testing family.
PURITY (v3 amendment -- the CPE 507/507 premise was disproven at value level; AMENDMENT_v3.md):
  PRIMARY = Aran et al. 2015 consensus CPE (PMID 26634437), COMPLETE-CASE n=506. The one
  patient with a non-finite CPE (TCGA-BS-A0TG, NSMP; CPE/ESTIMATE/ABSOLUTE all NaN) is
  EXCLUDED ONLY from models that include CPE, reason recorded exactly as "primary purity
  covariate non-finite"; in the CPE-adjusted model NSMP=146, the other three unchanged
  (POLE49/MMRd148/p53abn163). NO imputation, NO missingness indicator, NO estimator
  substitution. TCGA-BS-A0TG STAYS in the frozen cohort registry, in scoring, in the
  technical ledgers, and in the no-purity sensitivity; it is NOT replaced by IHC.
  SENSITIVITY (full-cohort robustness) = the frozen NO-PURITY model on all n=507, compared
  to the n=506 CPE-adjusted primary for concordance of effect DIRECTION, MAGNITUDE, and
  CATEGORY -- a robustness analysis for purity adjustment, NOT a test of the excluded
  patient's effect.
  SENSITIVITY (alternative estimator) = TCGA PanCanAtlas ABSOLUTE (purity/ploidy master
  calls, TCGA_mastercalls.abs_tables_JSedit.fixed.txt, sha256 f430a975..., column
  "purity", GDC UUID 4f277128-f793-4354-a13d-30cc7fe9f6b5), COMPLETE-CASE n=502 (POLE48/
  MMRd148/NSMP144/p53abn162). Chosen by PROVENANCE/DEFINITION (matches the spec name "TCGA
  PanCanAtlas ABSOLUTE"), NOT by coverage/result; the Aran-2015 ABSOLUTE column (n=366,
  older partial tabulation) is rejected on provenance. See ABSOLUTE_PROVENANCE.md.
  IHC (0.85 for TCGA-BS-A0TG) may enter ONLY a separate, pre-specified estimator-
  HARMONIZATION sensitivity with a validated cross-estimator mapping -- NEVER a single raw
  IHC substitution into the CPE column.
  REJECTED: mean-imputation + missingness indicator as primary; a "best-available-per-
  patient" composite purity variable; a heterogeneous per-patient estimator fallback.
COMPOSITION = exact ESTIMATE stromal + immune signatures (Yoshihara 2013, PMID 24113773;
  282 genes, frozen). Gene overlap with the modules is REMOVED PRE-SCORING and logged
  (B1 exclusion ledger: M1 -3, M3 -26). EXPLICIT SCOPE: ESTIMATE controls BULK STROMAL +
  IMMUNE admixture only; it does NOT measure the full epithelial / cell-type composition
  (epithelial fraction and specific cell-type/immune-subset proportions remain
  UNMEASURED -- a disclosed limitation).
COLLINEARITY (frozen VIF rule, evaluated at analysis on the covariate design columns
  only, WITHOUT examining subtype effects): VIF <= 5 acceptable; VIF > 5 flagged;
  VIF > 10 -> model STOP / revision before any coefficient interpretation. Remedy on
  collinearity: NEVER drop the subtype factor or a primary module score; remove ONLY a
  pre-labelled REDUNDANT covariate by the recorded priority order (lowest-priority
  first: composition-immune, then composition-stromal, then M4; purity retained), and
  only if VIF>10 persists after flagging. If still >10, revise the spec and report --
  do not interpret coefficients.

--------------------------------------------------------------------------------
B2.7 MULTIPLICITY FAMILIES (Vladimir point 7 / F4 restriction)
--------------------------------------------------------------------------------
Disjoint BH-FDR families (q<=0.05 confirmatory; q<=0.10 sensitivity-only):
  F1  module OMNIBUS subtype tests (M1_full; each M3 CORE_TF regulon)
  F2  module x planned-contrast {C1,C2,C3} tests (only for F1-gated modules)
  F3  TF-MODULE DISCORDANCE (D_t) tests (permutation-based, controlled within F3)
  F4  CONFIRMATORY EQUIVALENCE tests -- contains ONLY the pre-admissible equivalence
      tests = {module x C1} (per B2.9 H3, only C1 is confirmatory-equivalence-eligible).
      C2/C3 equivalence are SENSITIVITY-ONLY and are NOT in F4.
M4 enters NO developmental-testing family. All sensitivity analyses (compact-M1,
M3-sensitivity regulons, signed-ssGSEA, C2/C3 & POLE equivalence, interaction-framed
discordance, ABSOLUTE-purity, no-purity) are reported OUTSIDE the confirmatory families,
explicitly labelled sensitivity.

--------------------------------------------------------------------------------
B2.8 H1-H4 (Vladimir point 6) + power
--------------------------------------------------------------------------------
A-priori power (data-free; power_contrasts_table.tsv). Difference power at d=0.5:
C1 0.9995, C2 0.9955, C3 0.8551. TOST equivalence power at margin 0.30 (true d=0):
C1 0.869, C2 0.732, C3 0.139.
  H1  biological effect floor |d| >= 0.50 for C1-C3, reported with effect estimate +
      patient-bootstrap CI. q<=0.05 is the STATISTICAL (not biological) criterion.
      A NON-SIGNIFICANT difference test is NEVER called a null.
  H2  BH q<=0.05 confirmatory per family; q<=0.10 sensitivity.
  H3  equivalence / clean-null margin |d| < 0.30 is CONFIRMATORY for C1 ONLY (power
      0.869). C2 and C3 are equivalence SENSITIVITY-ONLY -- reported via CI/interval
      width, with NO "no effect" conclusion unless the CI lies ENTIRELY within the SESOI
      zone. There is NO wider POLE-specific SESOI. Clean-null = a TOST/CI-in-SESOI PASS,
      never a non-significant p; the inconclusive/underpowered bucket is kept DISTINCT
      from equivalence-null (R-19).
  H4  TF-MODULE DISCORDANCE D_t (Vladimir point 2):
      Residualization model (NO molecular subtype): for TF t,
         A_{t,i} = b0 + b1*z(TFexpr)_{t,i} + b2*M4_i + b3*purity_i + b4*composition_i + R_{t,i}
      where A = signed regulon activity (VIPER/aREA) z-standardized across 507, z(TFexpr)
      = z-standardized log2(TPM+1) of TF t. R is in units of SD(A). SUBTYPE IS NOT IN
      THIS MODEL.
      D_{t,c} = standardized contrast c on the residual R (same full-model contrast form
      as B2.4), unit = SD of the residual; |D_t|>=0.5 = ">= half-SD subtype shift in
      regulon activity unexplained by the TF's own abundance + frozen covariates."
      DECLARE TF-MODULE DISCORDANCE for (t,c) iff ALL FOUR hold:
         (i)   |D_{t,c}| >= 0.50;
         (ii)  patient-bootstrap CI on D_{t,c} EXCLUDES 0;
         (iii) permutation p passes the F3 multiplicity correction;
         (iv)  DIRECTION: the TF-EXPRESSION contrast e_{X,t,c} (same contrast c on the
               covariate-adjusted subtype means of z(TFexpr)) is itself non-trivial
               (|e_X|>=0.50 with its bootstrap CI excluding 0) AND sign(D_{t,c}) is
               OPPOSITE to sign(e_{X,t,c}).
      If (i)-(iii) hold but (iv) fails, the result is reported as "regulon activity
      UNEXPLAINED BY TF ABUNDANCE," NOT as TF-module discordance.
      Interaction framing (subtype x z(TF) in A ~ subtype*z(TF)+covariates) = sensitivity.
      EXPRESSION-level only -> consistent with, NOT demonstrating, retargeting/cofactor-
      switching (B3 firewall).

--------------------------------------------------------------------------------
B2.9 OUTLIER / MISSINGNESS
--------------------------------------------------------------------------------
No post-hoc SAMPLE exclusion from the FROZEN COHORT (507; registry/scoring unchanged).
Module gene missing after mapping -> dropped+logged; a module below the frozen minimum
(>=8) is DISCLOSED (R-18), not patched. PURITY MISSINGNESS (v3): the CPE-adjusted PRIMARY
is COMPLETE-CASE n=506 -- TCGA-BS-A0TG (non-finite CPE) is dropped from CPE-INCLUDING
MODELS ONLY (reason "primary purity covariate non-finite"), NOT from the cohort, scoring,
the technical ledgers, or the no-purity / other-covariate models. NO imputation, NO
missingness indicator, NO estimator substitution, NO tumour-value imputation.

--------------------------------------------------------------------------------
B2.10 SENSITIVITY / REPLICATION + P3 + SECONDARY-B (Vladimir points 8, 9)
--------------------------------------------------------------------------------
EXTERNAL REPLICATION is NOT a prerequisite; its ABSENCE is recorded a priori. Results
are "INTERNALLY BOOTSTRAPPED," never "validated." Internal robustness = patient-level
bootstrap + leave-one-DONOR-out / PCW-window split on the M1 selection (labelled
stability, not replication). K_STUDIES=1 is stamped on every M1 artifact (recon found
no suitable independent second fetal Mullerian reference; DEFER). An external cohort
(CPTAC-UCEC/other) is flagged as a future option, not a gate.
P3 (doctrine ecological prohibition) -- DEFINED: ALL primary and secondary CLAIMS are
SUBTYPE-LEVEL. NO per-patient biomarker claim is made; per-patient module scores are
intermediate quantities only, never an individual-level claim. The cohort is
one-per-patient: NO aliquot is treated as an independent observation and NO duplicate/
aliquot sensitivity inflates n. (There is no ambiguous "patient-level secondary" -- it
is excluded as a claim; only subtype-level inference is made.)
SECONDARY-B (35 peri-tumoral normals): post-primary-freeze DESCRIPTIVE sensitivity only;
not in selection/thresholds/category assignment/decision gate; field-effect-stamped.

--------------------------------------------------------------------------------
B2.11 DETERMINISTIC SEEDS
--------------------------------------------------------------------------------
Master SEED = 20260713, recorded here. Any bootstrap/permutation sub-seed is derived
DETERMINISTICALLY from (master_seed, analysis_step_id) and logged; no wall-clock or
nondeterministic seeding. Every stochastic result is reproducible from the master seed.

--------------------------------------------------------------------------------
B2.12 DECISION GATE
--------------------------------------------------------------------------------
Per (module, contrast): CONFIRMATORY category (enrichment/depletion CAT1/2, or clean-
null CAT5 via H3, or discordance CAT4 for M3) requires: module omnibus F1-pass (gate) ->
contrast |d|>=0.50 (H1) AND BH q<=0.05 (H2, family F2/F3/F4) AND survival of M4/purity/
composition adjustment AND (for a positive) the gene leave-one-out. Clean-null only via
the H3 equivalence rule (confirmatory C1 only). Escalation to the DEFERRED Level-1 ATAC
design only on the pre-registered trigger (subtype-specific AND robust-after-correction
AND TF-module discordance AND not reducible to a single TF/subtype AND internally-
bootstrapped-stable). Predictable/null -> recorded honestly, no escalation.
