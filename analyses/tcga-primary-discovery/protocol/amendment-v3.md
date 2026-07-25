TASK-028 FREEZE B -- AMENDMENT v3 NOTE (primary-purity completeness; OPTION 2 complete-case)
================================================================================
Directed by Vladimir 2026-07-13, after the Phase-1 locked-execution STOP (data-version
mismatch on primary purity). v2 (sealed_v2/) is preserved UNCHANGED as read-only archive
version 2; this v3 supersedes it for the resumed run. No tumour value informs any threshold;
git 83503bad untouched; SEED 20260713.

--------------------------------------------------------------------------------
REASON
--------------------------------------------------------------------------------
The v2 pre-seal covariate audit measured purity coverage by ROW PRESENCE only (it ran in
the blind, pre-tumour-read state and could not read purity VALUES), and recorded Aran CPE
507/507 and TCGA-PanCanAtlas-ABSOLUTE 503/507. The FIRST authorized value-level QC (Phase 1
of the locked run) found the CPE cell for exactly one patient is a literal NaN, so finite CPE
coverage is 506/507 -- the sealed premise "CPE complete -> no imputation" is false at the
value level. Every in-spec resolution was forbidden (imputation / missingness-indicator /
sample-drop / substitution), so the experimenter correctly STOPPED before model fitting and
returned a BLOCKED data-version mismatch. Sophia independently reconfirmed 506/507.

PATIENT: TCGA-BS-A0TG, subtype NSMP. Its CPE, ESTIMATE, and ABSOLUTE cells are all NaN;
only IHC=0.85 exists.

--------------------------------------------------------------------------------
DECISION (Vladimir): OPTION 2 -- COMPLETE-CASE (not a heterogeneous estimator fallback)
--------------------------------------------------------------------------------
Why complete-case over a CPE->IHC->ESTIMATE->ABSOLUTE fallback: a fallback would give this
one patient a DIFFERENT purity ESTIMATOR (IHC) than the other 506 (CPE), producing a
heterogeneous, cross-estimator-mixed purity covariate on a single scale -- unnecessary and
not validated. Complete-case keeps ONE homogeneous estimator in the primary model; the
patient is NOT lost (stays in the cohort registry, scoring, ledgers, and the no-purity
model). Resolution:
  - PRIMARY purity-adjusted model = Aran CPE COMPLETE-CASE n=506. TCGA-BS-A0TG excluded ONLY
    from CPE-including models, reason recorded exactly "primary purity covariate non-finite".
    CPE-adjusted subtype counts: POLE=49 / MMRd=148 / NSMP=146 / p53abn=163 (only NSMP -1).
  - FULL-COHORT NO-PURITY model = n=507, a PRE-LABELLED robustness sensitivity comparing
    effect DIRECTION / MAGNITUDE / CATEGORY vs the n=506 primary. It is NOT a substitute for
    primary and NOT a test of the single excluded patient's effect (that difference is a
    purity-adjustment robustness analysis, not a per-patient effect test).
  - NO imputation, NO missingness indicator, NO estimator substitution; TCGA-BS-A0TG is NOT
    replaced by IHC in the CPE column. IHC may enter ONLY a separate pre-specified estimator-
    harmonization sensitivity with a validated cross-estimator mapping.

--------------------------------------------------------------------------------
ABSOLUTE SENSITIVITY -- provenance resolved (before this seal)
--------------------------------------------------------------------------------
The v2 "ABSOLUTE 503/507" was also a ROW-PRESENCE figure; value-level finite = 502 (the extra
non-finite is TCGA-A5-A1OH / p53abn, an uncalled row). Two ABSOLUTE sources exist on disk:
the Aran-2015 ABSOLUTE COLUMN (n=366; an older, partial tabulation with whole cancer types
unpopulated) and the TCGA PanCanAtlas ABSOLUTE master calls (n=502). Chosen by PROVENANCE/
DEFINITION (matches the sealed spec name "TCGA PanCanAtlas ABSOLUTE"; single named file,
single "purity" column, pinnable checksum, GDC UUID) -- NOT by coverage or result. The Aran
column is rejected on provenance. ABSOLUTE sensitivity is FROZEN as complete-case n=502
(POLE48/MMRd148/NSMP144/p53abn162). Full detail: ABSOLUTE_PROVENANCE.md.
  pin: TCGA_mastercalls.abs_tables_JSedit.fixed.txt, sha256 f430a975..., column "purity",
       GDC UUID 4f277128-f793-4354-a13d-30cc7fe9f6b5.

--------------------------------------------------------------------------------
OLD -> NEW FILE MAP
--------------------------------------------------------------------------------
AMENDED (SHA differs from v2):
  B2_MODEL_SPEC.md   ONLY the purity/missingness/sample-flow lines: B2.4 (n=506 primary /
                     n=507 no-purity), B2.6 PURITY block, B2.9 purity-missingness.
  B1_PROVENANCE.md   ONLY the covariate-provenance line (CPE 506 / ABSOLUTE 502 pins).
NEW (not in v2):
  ABSOLUTE_PROVENANCE.md (the provenance resolution, byte-copied from execution/absolute_provenance/)
  AMENDMENT_v3.md (this note)
  SEAL_MANIFEST.sha256 (v3 seal manifest, new timestamp)
UNCHANGED -- BYTE-IDENTICAL TO v2 (proven in SEAL_MANIFEST.sha256 identity table):
  ALL B1/* (M1*, compact-M1*, M3 primary/analysis-ready members 603/577, m3_primary_edge_
  ledger 761/S1/S2, orphan_target_ledger, M3_sensitivity, M4, gene-level + exclusion ledgers,
  build script, B1_SHA256.manifest); B3_CLAIM_MATRIX.md; m3_edge_reconciliation.md;
  AMENDMENT_v2.md; amend_scripts/*.

--------------------------------------------------------------------------------
SCOPE GUARANTEE (UNCHANGED by v3)
--------------------------------------------------------------------------------
Module membership + scoring; subtype definitions; contrasts C1-C3 + weights; thresholds,
families F1-F4, SESOI/H1-H3, D_t; the all-zero rule; seeds; the claim firewall (B3). The v3
change is PURELY the purity/missingness/sample-flow specification + the ABSOLUTE-sensitivity
provenance pin. v1 and v2 remain archived, read-only.
