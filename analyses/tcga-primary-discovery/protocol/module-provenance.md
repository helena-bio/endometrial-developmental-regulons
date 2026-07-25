TASK-028 -- B1 MODULE PROVENANCE (FROZEN CANDIDATE for seal)
================================================================================
Read-only on seal; SHA-256 of every B1 file in B1/B1_SHA256.manifest and the top-level
SEAL_MANIFEST.sha256. All reference computation is tumour-blind. Deterministic (SEED
20260713; re-run -> identical SHAs, verified twice).

--------------------------------------------------------------------------------
FROZEN MODULE SET (sealed files under B1/)
--------------------------------------------------------------------------------
  M1  "E-MTAB-15475-derived fetal Mullerian epithelial module" (NOT a universal fetal-
      uterus program; K_STUDIES=1):
        M1_source.txt              239   sha256 02be0133...b98a1
        M1_analysis_ready.txt      236   sha256 fb32ca4f...e64ce   (239 minus composition
                                                                    overlap {CD74,SAMHD1,SPON1})
        compact_M1_source.txt      188   sha256 e96de995...f533b   SENSITIVITY-only
        compact_M1_analysis_ready  185   sha256 f517fb95...b33d96  (188 minus the same 3;
                                                                    all 3 were present; no replacement)
  M3  developmental-TF regulon, CollecTRI-SIGNED (PRIMARY):
        M3_primary_source.txt      603   sha256 79f681c9...8171f4  (20 CORE_TF + 583 targets)
        M3_primary_analysis_ready  577   sha256 49ff8951...dd6da5  (603 minus 26 composition)
        m3_primary_edge_ledger.tsv 843   sha256 e3044c91...80aaa1  (per-edge sign+weight)
        M3_sensitivity.txt         809   sha256 b67cbf62...575d69c  PRE-FROZEN SENSITIVITY
                                                                    (DoRothEA A/B + CollecTRI ChIP/motif)
  M4  proliferation control (COVARIATE-ONLY):
        M4_covariate.txt            30   sha256 055c9fef...0dd7c98
  Ledgers: B1_gene_level_ledger.tsv (1869 rows) sha256 64cd156b...6f3e7f ;
           B1_exclusion_ledger.tsv  (32 rows)   sha256 1610c3b7...f3f4e0
  Generator: build_sealed_B1.py sha256 8a7460c0...1fc28 (python3 stdlib; run twice, identical).

--------------------------------------------------------------------------------
M3 SIGNED-EDGE RULE -- RESOLVED (v2 amendment, DEC-M3-DUAL-EDGE, Vladimir 2026-07-13)
--------------------------------------------------------------------------------
Edge ledger m3_primary_edge_ledger.tsv: VIPER/aREA mode-of-regulation = CollecTRI consensus
sign; per-edge weight = CollecTRI consensus confidence (1.0 consensus / 0.5 non-consensus).
SOURCE-TRUE ACCOUNTING (from the pinned CollecTRI; NOT the earlier imprecise figure):
  802 admitted edges = 723 single-flag + 79 both-flagged.
  By CONSENSUS resolution: 761 consensus-resolved (PRIMARY) + 41 no-consensus (EXCLUDED).
  79 both-flagged = 59 resolved + 20 unresolved.  723 single-flag = 702 resolved + 21 unresolved.
  Excluded-41 = 20 both-flagged-unresolved + 21 single-flag-no-consensus.
PRIMARY POLICY (consensus-only): the 761 consensus-resolved edges enter primary ONCE with
  their consensus sign; the 41 no-consensus edges are EXCLUDED (not split, no raw-flag sign,
  no inferred sign). S1 (flag-split, all 79 both-flagged) and S2 (exclude all 79 both-flagged)
  are pre-frozen sensitivity sets (see B2.3). ORPHANS: 21 unique non-CORE targets are primary-
  edgeless (22 excluded incoming edges; NFKB1 has 2), receive NO contribution, listed in
  B1/orphan_target_ledger.tsv with reason "no consensus-signed primary edge"; PAX2/TBX18 are
  CORE_TFs that lose all incoming edges as targets but remain regulators (separate note).
  Target counts: source 583 / primary-with-edge 562 / orphan 21. Member files 603/577 UNCHANGED.
NOTE ON v1: version-1 (archived, read-only) used "dual" as an IMPRECISE label for the 41 no-
  consensus edges and split each into +/- (843-row ledger). That was superseded here: the true
  both-flagged set is 79 (not 41), and no-consensus edges are EXCLUDED from primary, not split.
  This amendment changes ONLY edge eligibility + the ledger; the admitted edge SET and all
  module membership are byte-identical to v1 (proven in AMENDMENT_v2.md).

--------------------------------------------------------------------------------
SOURCE PROVENANCE (pinned; tumour-blind; selection/freezing only)
--------------------------------------------------------------------------------
FETAL (M1): E-MTAB-15475 (Nature 2025, PMID 41407855, DOI 10.1038/s41586-025-09875-2);
  corrected canonical pins (D10 erratum): study a1c36710... (104223 B, ReleaseDate
  2025-08-07); SDRF e7f7ebc8...769007 (82933 B); IDF fceaec52...0401d0 (7952 B). Annotated
  object post10pcw_females.h5ad (open Reproductive Cell Atlas; 227932 cells / 28 donors)
  -> 15158 Mullerian-epithelial cells / 25 Mullerian-epi donors. D10 forensic report kept
  as PROVENANCE ERRATUM: experimenter_d10/D10_RESOLUTION.md (mis-attribution, not drift;
  byte-stable + upstream-consistent; no rebuild).
M3 regulons: CollecTRI (Muller-Dott NAR 2023, PMC10639077) SIGNED = PRIMARY curated
  source; DoRothEA A/B + GTRD/ChIP/motif = PRE-FROZEN SENSITIVITY (not merged). Pinned in
  MODULE_FINAL_MANIFEST.json.
Filters: proliferation = MSigDB Hallmark E2F+G2M / Tirosh cc.genes; housekeeping = HRT
  Atlas v1.0; non-specific = in-atlas pan-compartment set + blood/immune guard.
Covariate provenance (for B2; v3 amendment): purity PRIMARY = Aran 2015 CPE (PMID 26634437)
  COMPLETE-CASE n=506 (1 non-finite: TCGA-BS-A0TG/NSMP; excluded from CPE-including models
  only). Purity SENSITIVITIES = full-cohort NO-PURITY model n=507; TCGA PanCanAtlas ABSOLUTE
  (TCGA_mastercalls.abs_tables_JSedit.fixed.txt, sha256 f430a975..., column "purity", GDC
  UUID 4f277128-f793-4354-a13d-30cc7fe9f6b5) COMPLETE-CASE n=502; IHC only in a pre-specified
  estimator-harmonization sensitivity. (v2's "CPE 507/507 / ABSOLUTE 503/507" were ROW-
  PRESENCE counts; value-level finite is 506 / 502. See AMENDMENT_v3.md + ABSOLUTE_PROVENANCE.md.)
  composition = ESTIMATE (Yoshihara 2013, PMID 24113773; gmt bfd34f3f..., 282 genes),
  overlap removed pre-scoring (M1 -3, M3 -26; logged).

--------------------------------------------------------------------------------
K_STUDIES=1 (recon before seal -- Vladimir point 8)
--------------------------------------------------------------------------------
Independent source-level reconnaissance for a SECOND open, dissociated, fetal-Mullerian-
epithelial scRNA reference found NONE suitable (RECON_SECOND_FETAL_REF.md): Taelman
GSE181558 = ovary/testis/mesonephros, no Mullerian uterus; He 2026 = single-nucleus +
Stereo-seq + controlled-access; Descartes/HDCA = no uterus; others mouse/spatial/adult.
VERDICT: K_STUDIES=1 CONFIRMED -> DEFER external replication; proceed INTERNALLY
BOOTSTRAPPED (donor bootstrap / leave-one-donor-out / PCW split as stability, not
replication). No weak surrogate constructed (a surrogate would falsely present K=1 as
K=2). "K_STUDIES=1" is stamped on every M1 artifact.

--------------------------------------------------------------------------------
M2 -- DEFERRED (not part of B1)
--------------------------------------------------------------------------------
M2 (adult differentiated endometrial) is NOT frozen (DEC-M2): as built it was 96.4%
menstrual/secretory-state driven (measures a cyclic state, not a stable adult identity);
the composition/menstrual removal collapsed it below the R-16 minimum. NOT rescued by
post-hoc regression or the residual. Reinstatement ONLY via a future independent stage-
balanced adult-endometrium reference + a cycle-invariant module (pre-specified stage-
control, sufficient donor support, pre-fixed minimum size, independent not-hormone/
menstrual check) built BEFORE any tumour read.

--------------------------------------------------------------------------------
INTEGRITY
--------------------------------------------------------------------------------
Reference-atlas + metadata computation only; NO tumour value read; 35 normals off-limits;
git 83503bad untouched; ASCII; determinism re-verified. Anti-fabrication note: 17/27 M1
literature anchors were REJECTED on the atlas data (the reference overrode the literature).
