TASK-029 (C): M1 (F2 NEGATIVE) EQUIVALENCE POWER / PRECISION
============================================================
DATA-FREE. Exact Schuirmann two-one-sided-t (TOST) equivalence power at true d=0,
plus 95% CI half-width on d. Script: m1_equivalence_power.py .
JSON: m1_equivalence_power_results.json . Seed 20260714 (calc is deterministic).
Git HEAD unchanged: 83503bad47b60193598b2b9ebe819c22c83e8ac1.

M1 (fetal Mullerian-epithelial module) is the F2 NEGATIVE target, tested as an
EQUIVALENCE (TOST / CI-in-SESOI) hypothesis, SEPARATE from the F1 BH-of-6 family.
Proposed margin = TASK-028 SESOI 0.30; sensitivity margins 0.20 and 0.40.

------------------------------------------------------------------------------
METHOD CONSISTENCY (validated against the frozen TASK-028 numbers)
------------------------------------------------------------------------------
The exact-t TOST used here REPRODUCES TASK-028's reported equivalence table:
  TASK-028 C1 (163/344) @0.30: this method 0.868  vs  028 reported 0.869
  TASK-028 C2 (197/147) @0.30: this method 0.729  vs  028 reported 0.732
  CI half-width: matches exactly (C1 0.186, C2 0.214).
=> same method, no optimistic normal-approx inflation. (A naive normal-approx TOST
   would have over-stated these by ~0.06-0.13; not used.)

------------------------------------------------------------------------------
RESULTS -- REAL RESOLVED POOLED CPTAC SIZES (primary)
------------------------------------------------------------------------------
  contrast                 SE0     pw@0.20  pw@0.30  pw@0.40   95%CI half-width
  C1 real (36 vs 194)     0.1815    0.000    0.001    0.419        0.356
  C2 real (85 vs 109)     0.1447    0.000    0.325    0.732        0.284

  C1 SPEC-expected (50/183) 0.1596  0.000    0.181    0.607        0.313
  C2 SPEC-expected (78/95)  0.1528  0.000    0.243    0.664        0.299

Single-cohort context (NOT the pooled design):
  C1 confirmatory (16/119): @0.30 = 0.000 ; C1 discovery (20/75): @0.30 = 0.000
  C2 confirmatory (53/66):  @0.30 = 0.000 ; C2 discovery (32/43): @0.30 = 0.038

------------------------------------------------------------------------------
INTERPRETATION
------------------------------------------------------------------------------
- At margin 0.30 on the REAL pooled counts, M1 equivalence power is ~0.001 on C1
  and 0.325 on C2. Neither meets the >=0.80 confirmatory bar.
  (C1's ~0 is exact, not an artifact: a = m/SE0 = 0.30/0.1815 = 1.653 barely exceeds
   t_crit = 1.652, so the TOST acceptance window is essentially empty at n=230.)
- Margin needed for 0.80 power on the REAL counts: 0.533 (C1) and 0.425 (C2) --
  BOTH wider than the SESOI. Widening the margin to manufacture a PASS is forbidden
  (mirrors TASK-028 R: "NO wider margin to manufacture a PASS").
- CI half-widths on d are 0.356 (C1) and 0.284 (C2) -- wider than TASK-028's 0.186 /
  0.214 because CPTAC is smaller than the TCGA cohort.

FEASIBILITY VERDICT
  C1 (p53abn vs rest): M1 confirmatory equivalence NOT feasible at CPTAC sizes
      (power ~0 @0.30). SENSITIVITY / inconclusive ONLY.
  C2 (POLE+MMRd vs NSMP): M1 confirmatory equivalence NOT feasible @0.30
      (power 0.325); reaches ~0.73 only at margin 0.40. SENSITIVITY ONLY.
  => On CPTAC, M1 can be reported as an equivalence SENSITIVITY on both contrasts,
     but CANNOT carry a CONFIRMATORY replicated-negative at the SESOI 0.30.

REINFORCE (anti-fabrication): a non-significant DIFFERENCE test is NOT a replicated
negative. A replicated M1 negative requires a TOST/CI-in-SESOI PASS (CI entirely
inside +/-margin). The underpowered/inconclusive bucket is kept DISTINCT from an
equivalence-null. On CPTAC, an M1 "no difference" p>0.05 would be INCONCLUSIVE,
not a clean null.

FREEZE-DECISION ITEM (flagged, NOT authored here): the "all-subtype-contrasts vs
pooled" choice for M1 (i.e. whether M1 equivalence is tested per-contrast C1/C2 or
on a single pooled 4-class omnibus) is a design freeze decision. Recon only.
