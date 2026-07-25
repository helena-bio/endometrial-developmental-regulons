TASK-029 (B): SIMULATED BH-of-6 POWER for the F1 POSITIVE-REGULON FAMILY
========================================================================
DATA-FREE Monte-Carlo. No expression, no scoring. Arm sizes + true effect +
test structure only. Seed pinned = 20260714. N_rep = 200000 (>= 20000 required).
Script: simulate_power_f1.py . Full JSON: simulated_power_f1_results.json .
Git HEAD unchanged: 83503bad47b60193598b2b9ebe819c22c83e8ac1.

------------------------------------------------------------------------------
MODEL (pre-specified, no tuning)
------------------------------------------------------------------------------
F1 = SIX tests, all with true standardized effect d = 0.5 (the |d|>=0.5 floor,
used as the power alternative):
  C1 (p53abn vs pooled {POLE,MMRd,NSMP}): LHX1, PAX8
  C2 (pooled {POLE,MMRd} vs NSMP):        GATA2, HOXA9, SOX9, WT1
Per test i (two independent groups, Cohen's d):
  SE(d) = sqrt( (n1+n2)/(n1*n2) + d^2/(2*(n1+n2)) )
  d_hat_i ~ Normal(0.5, SE_i);  z_i = d_hat_i/SE_i;  two-sided p_i = 2*(1-Phi(|z_i|)).
Rejection = BH-significant AND correct direction (d_hat_i > 0).
Multiplicity: Benjamini-Hochberg q<=0.05 across all six (primary); Holm and
Bonferroni q<=0.05 reported as conservative sensitivity bounds; nominal alpha=0.05
per-test as the upper reference.
Tests treated as INDEPENDENT (default). CAVEAT: positive correlation among the four
C2 regulons would generally RAISE effective BH power (a stated caveat, NOT a change).

VALIDATION: the vectorized BH step-up matches statsmodels multipletests(fdr_bh)
EXACTLY (per-target power identical to 4 dp); RNG draws are reproducible under the
pinned seed (verified).

RULE: simulated BH power < 0.80 -> `underpowered` (even if all subtype arms >= 20).

------------------------------------------------------------------------------
PRIMARY SCENARIO: REAL RESOLVED POOLED SIZES (from part A)
  C1 p53abn=36 vs rest=194 ;  C2 POLE+MMRd=85 vs NSMP=109
------------------------------------------------------------------------------
  target        nominal   BH-sim    Holm    Bonf    status
  LHX1  (C1)     0.780    0.766    0.706   0.538    UNDERPOWERED
  PAX8  (C1)     0.780    0.766    0.706   0.538    UNDERPOWERED
  GATA2 (C2)     0.925    0.916    0.876   0.777    OK (>=0.80)
  HOXA9 (C2)     0.926    0.916    0.876   0.777    OK (>=0.80)
  SOX9  (C2)     0.926    0.917    0.876   0.778    OK (>=0.80)
  WT1   (C2)     0.926    0.917    0.877   0.778    OK (>=0.80)
  family all-6 BH power = 0.446 ; mean #rejected(BH) = 5.20 / 6

VERDICT (primary/real): 4 of 6 targets are adequately powered under BH-of-6
(all four C2 regulons, ~0.92). The TWO C1 targets (LHX1, PAX8) are UNDERPOWERED
(BH 0.766 < 0.80), driven by the small p53abn arm (n=36 real vs ~50 anticipated).

------------------------------------------------------------------------------
SENSITIVITY 1: SPEC EXPECTED SIZES  (C1 50/183 ; C2 78/95)
------------------------------------------------------------------------------
  target        nominal   BH-sim    Holm    Bonf    status
  LHX1  (C1)     0.874    0.863    0.815   0.677    OK (>=0.80)
  PAX8  (C1)     0.873    0.862    0.813   0.677    OK (>=0.80)
  GATA2 (C2)     0.896    0.886    0.842   0.720    OK (>=0.80)
  HOXA9 (C2)     0.896    0.886    0.841   0.718    OK (>=0.80)
  SOX9  (C2)     0.897    0.887    0.842   0.721    OK (>=0.80)
  WT1   (C2)     0.897    0.887    0.843   0.721    OK (>=0.80)
  family all-6 BH power = 0.492
Under the ANTICIPATED sizes all six clear 0.80. The only reason C1 fails on the
REAL data is the p53abn arm being 36, not ~50. This is the load-bearing gap.

------------------------------------------------------------------------------
SENSITIVITY 2 & 3: single-cohort arms (NOT the pooled design; for context)
------------------------------------------------------------------------------
Confirmatory-only (C1 16/119, C2 53/66): ALL SIX underpowered (BH 0.40-0.70).
Discovery-only    (C1 20/75,  C2 32/43): ALL SIX underpowered (BH 0.40-0.45).
=> Neither cohort alone can carry F1; the design MUST pool Discovery+Confirmatory,
   and even pooled the C1 pair is short of 0.80 on the real counts.

------------------------------------------------------------------------------
FINAL SIMULATED-BH POWER TABLE (PRIMARY = REAL pooled)
------------------------------------------------------------------------------
  target   nominal  BH-sim  Holm   Bonf   status
  LHX1     0.780    0.766   0.706  0.538  underpowered
  PAX8     0.780    0.766   0.706  0.538  underpowered
  GATA2    0.925    0.916   0.876  0.777  ok
  HOXA9    0.926    0.916   0.876  0.777  ok
  SOX9     0.926    0.917   0.876  0.778  ok
  WT1      0.926    0.917   0.877  0.778  ok

DESIGN IMPLICATION (recon, not a freeze): under the REAL resolved CPTAC counts and
the frozen d=0.5 floor, the C2 regulon replication (GATA2/HOXA9/SOX9/WT1) is
adequately powered under BH-of-6 (~0.92 each); the C1 pair (LHX1/PAX8) is
underpowered (~0.77) and should be flagged accordingly before the interval is
replaced/frozen. Positive C2-regulon correlation would only help; the C1 shortfall
would remain. NO retuning of d, alpha, or margins was done to move these numbers.
