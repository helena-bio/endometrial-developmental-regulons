#!/usr/bin/env python3
"""
TASK-029 a-priori DIFFERENCE-DETECTION power + CI-half-width on Cohen's d.

DATA-FREE by construction. This script reads NO cohort data, NO expression
values, NO subtype labels, NO matrix. Its ONLY inputs are integer arm sizes
(n1, n2) and a target standardized effect size |d|. It computes, from standard
closed-form two-sample formulas:

  (1) DIFFERENCE-detection power: probability that a two-sample test rejects
      H0: d=0 at a two-sided level alpha, when the true effect is |d|=0.5.
      Uses the noncentral-t distribution (exact small-sample form), which is
      the correct reference for a two-sample t/standardized-mean-difference
      test; the large-sample normal approximation is printed alongside as a
      cross-check.

  (2) EXPECTED 95% CI HALF-WIDTH on the estimated d (Cohen's d), from the
      large-sample variance of d (Hedges & Olkin 1985):
         Var(d) = (n1+n2)/(n1*n2) + d^2 / (2*(n1+n2))
         SE(d)  = sqrt(Var(d));  half-width = z_{0.975} * SE(d)
      This half-width is the a-priori PRECISION statement: how tightly the
      replication would pin d, independent of the observed value.

alpha handling: the frozen replication uses BH q<=0.05 across the six targets.
For an a-priori per-target power floor we report power at BOTH:
  - alpha = 0.05 (nominal per-test), and
  - alpha = 0.05/6 ~= 0.00833 (Bonferroni-conservative bound on the BH-adjusted
    family of six targets; BH is uniformly at least as powerful as Bonferroni,
    so the true frozen-procedure power lies between these two -- Bonferroni is
    the conservative floor).
No threshold is tuned; alpha and d are fixed constants from the frozen design.

Standard references:
  Cohen J. Statistical Power Analysis for the Behavioral Sciences, 2nd ed 1988.
  Hedges LV, Olkin I. Statistical Methods for Meta-Analysis, 1985 (Var(d)).
"""

import numpy as np
from scipy import stats

Z975 = stats.norm.ppf(0.975)


def se_d(n1, n2, d):
    """Large-sample SE of Cohen's d (Hedges-Olkin)."""
    var = (n1 + n2) / (n1 * n2) + (d * d) / (2.0 * (n1 + n2))
    return np.sqrt(var)


def ci_halfwidth(n1, n2, d, z=Z975):
    return z * se_d(n1, n2, d)


def power_noncentral_t(n1, n2, d, alpha):
    """Exact two-sample difference-detection power via noncentral t.

    df = n1 + n2 - 2 ; noncentrality ncp = d * sqrt(n1*n2/(n1+n2)).
    Two-sided test at level alpha.
    """
    df = n1 + n2 - 2
    ncp = d * np.sqrt(n1 * n2 / (n1 + n2))
    tcrit = stats.t.ppf(1.0 - alpha / 2.0, df)
    # P(reject) = P(T > tcrit) + P(T < -tcrit) under noncentral t(df, ncp)
    upper = stats.nct.sf(tcrit, df, ncp)
    lower = stats.nct.cdf(-tcrit, df, ncp)
    return upper + lower


def power_normal(n1, n2, d, alpha):
    """Large-sample normal-approx difference-detection power (cross-check)."""
    ncp = d * np.sqrt(n1 * n2 / (n1 + n2))
    zcrit = stats.norm.ppf(1.0 - alpha / 2.0)
    return stats.norm.sf(zcrit - ncp) + stats.norm.cdf(-zcrit - ncp)


ALPHA_NOMINAL = 0.05
ALPHA_BONF6 = 0.05 / 6.0
D = 0.5

# ---------------------------------------------------------------------------
# Scenarios. Arm sizes are DESIGN NUMBERS (expected counts from the recon),
# NOT read from any cohort. They are the contrast-group sizes the frozen design
# would face. POLE-thin single-cohort scenarios included to show the taxonomy.
# ---------------------------------------------------------------------------
scenarios = [
    # label, n1 (focus arm), n2 (comparison arm), note
    ("POOLED CPTAC  C1  p53abn vs pooled non-p53abn", 50, 183,
     "expected p53abn ~45-55; rest ~180 (POLE+MMRd+NSMP)"),
    ("POOLED CPTAC  C2  POLE+MMRd vs NSMP",           78,  95,
     "expected POLE+MMRd ~70-85; NSMP ~90-100"),
    ("POOLED CPTAC  M1-neg  p53abn vs rest (C1 axis)", 50, 183,
     "M1 must be powered to DETECT d=0.5 to certify a replicated negative"),
    ("POOLED CPTAC  M1-neg  POLE+MMRd vs NSMP (C2 axis)", 78, 95,
     "M1 second axis; same power question"),
    ("DISCOVERY only  C1  p53abn(20) vs rest(75)",    20,  75,
     "Dou2020 published CNV-high=20 vs 75"),
    ("DISCOVERY only  C2  POLE+MMRd(32) vs NSMP(43)", 32,  43,
     "Dou2020 POLE7+MSIh25=32 vs CNV-low43"),
    ("CONFIRMATORY only C1  p53abn(~30) vs rest(~108)", 30, 108,
     "expected only; not published per-subtype at recon"),
    ("CONFIRMATORY only C2  POLE+MMRd(~47) vs NSMP(~57)", 47, 57,
     "expected only"),
    # sensitivity across the plausible p53abn / pooled ranges for pooled C1
    ("POOLED CPTAC  C1 low end  p53abn(45) vs rest(188)", 45, 188, "range low"),
    ("POOLED CPTAC  C1 high end p53abn(55) vs rest(178)", 55, 178, "range high"),
    ("POOLED CPTAC  C2 low end  POLE+MMRd(70) vs NSMP(100)", 70, 100, "range low"),
    ("POOLED CPTAC  C2 high end POLE+MMRd(85) vs NSMP(90)",  85,  90, "range high"),
]

print("=" * 100)
print("TASK-029 A-PRIORI POWER / PRECISION  (DATA-FREE: inputs = n1, n2, |d|=%.2f only)" % D)
print("power = two-sample difference-detection power (noncentral-t exact; normal-approx cross-check)")
print("CI half-width = z0.975 * SE(d), SE(d)=sqrt((n1+n2)/(n1 n2) + d^2/(2(n1+n2)))  [Hedges-Olkin]")
print("alpha_nominal=%.4f ; alpha_Bonf6=%.5f (conservative BH-of-6 floor)" % (ALPHA_NOMINAL, ALPHA_BONF6))
print("=" * 100)
hdr = "%-52s %5s %5s | %7s %7s %7s | %8s" % (
    "scenario", "n1", "n2", "pw@.05", "pw@B6", "pw~N.05", "CIhalf_d")
print(hdr)
print("-" * 100)
for label, n1, n2, note in scenarios:
    pw05 = power_noncentral_t(n1, n2, D, ALPHA_NOMINAL)
    pwb6 = power_noncentral_t(n1, n2, D, ALPHA_BONF6)
    pwN = power_normal(n1, n2, D, ALPHA_NOMINAL)
    ch = ci_halfwidth(n1, n2, D)
    print("%-52s %5d %5d | %7.3f %7.3f %7.3f | %8.3f" %
          (label, n1, n2, pw05, pwb6, pwN, ch))
print("-" * 100)
print("NOTES:")
for label, n1, n2, note in scenarios:
    print("  %-52s : %s" % (label, note))
print()
print("INTERPRETATION KEY (feasibility taxonomy, applied in the v2 memo):")
print("  feasible          : both required arms >=20 AND difference-power >=0.80 at d=0.5")
print("                      (propose alpha per the frozen BH; report both .05 and Bonf6 bounds)")
print("  underpowered      : any required arm <20 OR power <0.80 at d=0.5")
print("  CI half-width      : the a-priori precision; a narrower half-width = tighter pinning of d")
print("DATA-FREE CONFIRMATION: this script opened NO cohort file; all n are design expectations.")
