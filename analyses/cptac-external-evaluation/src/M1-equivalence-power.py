#!/usr/bin/env python3
"""
TASK-029 (C): M1 (F2 negative) EQUIVALENCE power/precision. DATA-FREE.
TOST / CI-in-SESOI equivalence on a standardized effect (Cohen's d), true d = 0.

EXACT Schuirmann two-one-sided-t (TOST) power, df = n1+n2-2, at true d=0:
  se = SE0 of d = sqrt((n1+n2)/(n1*n2)); a = m/se; tcrit = t_{1-alpha, df}.
  S = dhat/se ~ noncentral-t(df, ncp = d_true/se).
  TOST rejects iff (tcrit - a) < S < (a - tcrit); requires a > tcrit.
  power = F_nct(a - tcrit) - F_nct(tcrit - a).
This exact-t form REPRODUCES the frozen TASK-028 values (C1 163/344 @0.30 = 0.868
vs 028's 0.869; C2 197/147 @0.30 = 0.729 vs 028's 0.732), i.e. same method, no
optimistic normal-approx inflation.
95% CI half-width on d = 1.96 * SE0 (matches TASK-028 exactly: C1 0.186, C2 0.214).

Sizes are the REAL resolved pooled contrast group sizes from part A (and single
cohorts for context). Margin proposal 0.30 (TASK-028 SESOI); 0.20 and 0.40 sensitivity.
Seed recorded for manifest completeness (calc is closed-form/deterministic).
"""
import numpy as np
from scipy.stats import norm, nct, t as tdist
import json

SEED = 20260714
ALPHA = 0.05
Z_CI = norm.ppf(0.975)         # 1.9600, for the reported 95% CI half-width

def se0(n1, n2):
    n1, n2 = float(n1), float(n2)
    return np.sqrt((n1 + n2) / (n1 * n2))   # SE of d at d=0

def tost_power_d0(n1, n2, m, d_true=0.0):
    df = n1 + n2 - 2
    se = se0(n1, n2)
    a = m / se
    tcrit = tdist.ppf(1 - ALPHA, df)
    if (a - tcrit) <= (tcrit - a):
        return 0.0
    ncp = d_true / se
    return float(nct.cdf(a - tcrit, df, ncp) - nct.cdf(tcrit - a, df, ncp))

def ci_halfwidth(n1, n2):
    return float(Z_CI * se0(n1, n2))

CONTRASTS = {
    # REAL resolved pooled (primary)
    "C1_pooled_real (36 vs 194)": (36, 194),
    "C2_pooled_real (85 vs 109)": (85, 109),
    # single-cohort context
    "C1_confirmatory (16 vs 119)": (16, 119),
    "C2_confirmatory (53 vs 66)": (53, 66),
    "C1_discovery (20 vs 75)": (20, 75),
    "C2_discovery (32 vs 43)": (32, 43),
    # spec expected (sensitivity)
    "C1_expected_spec (50 vs 183)": (50, 183),
    "C2_expected_spec (78 vs 95)": (78, 95),
}
MARGINS = [0.20, 0.30, 0.40]

if __name__ == "__main__":
    out = {"seed": SEED, "alpha": ALPHA, "margins": MARGINS, "results": {}}
    print("%-30s %8s | %s | %10s" %
          ("contrast", "SE0", "  ".join("pw@%.2f" % m for m in MARGINS), "CIhalfwid"))
    for name, (n1, n2) in CONTRASTS.items():
        se = se0(n1, n2)
        powers = {m: tost_power_d0(n1, n2, m) for m in MARGINS}
        hw = ci_halfwidth(n1, n2)
        out["results"][name] = {
            "n1": n1, "n2": n2, "SE0": se,
            "tost_power_true_d0": powers,
            "ci95_halfwidth_on_d": hw,
            "confirmatory_at_0.30": powers[0.30] >= 0.80,
        }
        print("%-30s %8.4f | %s | %10.3f" %
              (name, se,
               "  ".join("%6.3f" % powers[m] for m in MARGINS), hw))
    json.dump(out, open("m1_equivalence_power_results.json", "w"), indent=2)
