TASK-028 FREEZE B -- AMENDMENT v2 NOTE (DEC-M3-DUAL-EDGE)
================================================================================
Directed by Vladimir 2026-07-13. Narrow, versioned amendment closing the v1 critic's
dual-edge MINOR. v1 (sealed/) is preserved UNCHANGED as the read-only archive (version 1);
this v2 (sealed_v2/) supersedes it for the eventual run. No tumour value read; git 83503bad
untouched; determinism re-verified; SEED 20260713.

--------------------------------------------------------------------------------
REASON
--------------------------------------------------------------------------------
The v1 seal split "41 dual" edges into +/- signed contributions and used "dual" as an
IMPRECISE label. On amendment, the pinned CollecTRI source shows the true structure:
the both-flagged set is 79 (not 41), and the 41 is the NO-CONSENSUS set (20 both-flagged +
21 single-flag). Vladimir's stated 38/41 split was NOT reproducible from the source; the
experimenter reported the source-true numbers and did not force them. Vladimir then fixed a
STRICTLY CONSENSUS-SIGN-BASED primary policy + two pre-frozen sensitivity sets.

--------------------------------------------------------------------------------
SOURCE-TRUE ACCOUNTING (frozen)
--------------------------------------------------------------------------------
  802 admitted edges = 723 single-flag + 79 both-flagged
  consensus resolution: 761 resolved (PRIMARY) + 41 no-consensus (EXCLUDED)
  79 both-flagged = 59 resolved + 20 unresolved
  723 single-flag  = 702 resolved + 21 unresolved
  excluded-41      = 20 both-flagged-unresolved + 21 single-flag-no-consensus
  PRIMARY = 761 (consensus sign, once each); S1 = all 79 both-flagged split; S2 = exclude all 79
  target counts: source 583 / primary-with-edge 562 / orphan 21 (22 excluded incoming edges)

--------------------------------------------------------------------------------
OLD -> NEW FILE MAP
--------------------------------------------------------------------------------
AMENDED (SHA differs from v1):
  B1/m3_primary_edge_ledger.tsv    v1: 843 rows (41 no-consensus split +/-)
                                   v2: policy-tagged; 761 primary-included (once, consensus sign),
                                       41 excluded, per-edge S1/S2 fields; sha 96d3dd9a...607347
  B1/B1_SHA256.manifest            regenerated (edge-ledger SHA updated; orphan ledger added)
                                       sha 5915486b... -> after orphan add
  B2_MODEL_SPEC.md                 ONLY section B2.3 (M3 scoring) amended
  B1_PROVENANCE.md                 ONLY the M3 signed-edge section amended
NEW (not in v1):
  B1/orphan_target_ledger.tsv      sha c4391216...e9352 (21 orphan targets + 562 primary + CORE note)
  m3_edge_reconciliation.md        reconciliation + orphan section
  amend_scripts/                   deterministic generators
  AMENDMENT_v2.md                  this note
  SEAL_MANIFEST.sha256             v2 seal manifest (new timestamp)
UNCHANGED -- BYTE-IDENTICAL TO v1 (proven in SEAL_MANIFEST.sha256 identity table):
  B1/M1_source.txt, M1_analysis_ready.txt, compact_M1_source.txt, compact_M1_analysis_ready.txt,
  B1/M3_primary_source.txt (603), M3_primary_analysis_ready.txt (577), M3_sensitivity.txt,
  M4_covariate.txt, B1_gene_level_ledger.tsv, B1_exclusion_ledger.tsv, build_sealed_B1.py,
  B3_CLAIM_MATRIX.md.

--------------------------------------------------------------------------------
SCOPE GUARANTEE (unchanged by this amendment)
--------------------------------------------------------------------------------
M1/M4 membership; M3 TARGET membership (603/577 provenance universe); contrasts (C1-C3);
covariates; thresholds; multiplicity families; SESOI/H1-H3; D_t; seeds; claim firewall (B3).
The change is PURELY M3 primary-edge ELIGIBILITY + the sensitivity-set definitions + the
orphan ledger. The admitted edge SET is byte-identical to v1; membership is unchanged.
v1 remains the archived, read-only version-1 record.
