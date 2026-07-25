#!/usr/bin/env python3
"""
TASK-029 (B): SIMULATED BH-of-6 power for the F1 positive-regulon family.
DATA-FREE. No expression, no scoring. Arm sizes + true effect + test structure only.

F1 = SIX tests, effect floor |d|>=0.5, true d assumed = 0.5 (the floor) for power:
  C1 (p53abn vs pooled non-p53abn): LHX1, PAX8
  C2 (POLE+MMRd vs NSMP):           GATA2, HOXA9, SOX9, WT1
Primary multiplicity = BH FDR q<=0.05 ACROSS all six.

Per-test standardized-effect estimate model (two independent groups, Cohen's d):
  SE(d) = sqrt( (n1+n2)/(n1*n2) + d^2 / (2*(n1+n2)) )
  d_hat_i ~ Normal(true_d, SE_i);  z_i = d_hat_i / SE_i;  two-sided p_i = 2*(1-Phi(|z_i|)).
Rejection requires BH-significance AND correct direction (d_hat_i > 0).

Outputs: per-target and family power under (i) nominal alpha=0.05,
(ii) simulated BH-of-6, (iii) Holm, (iv) Bonferroni. Independence assumed
across the six tests (stated caveat: positive C2-regulon correlation would
generally RAISE effective BH power).

Rule: any target with simulated BH power < 0.80 is `underpowered`.
"""
import numpy as np
from scipy.stats import norm
import json, sys, hashlib

SEED = 20260714           # pinned
N_REP = 200_000           # >= 20000 required; use 200k for stable 3rd decimal
TRUE_D = 0.5              # effect floor as the power alternative
Q = 0.05                 # BH / Holm / Bonferroni family level
ALPHA = 0.05             # nominal per-test

def se_d(n1, n2, d):
    n1 = float(n1); n2 = float(n2)
    return np.sqrt((n1 + n2) / (n1 * n2) + d * d / (2.0 * (n1 + n2)))

def run(contrast_sizes, label):
    """contrast_sizes: dict target -> (n1, n2). Returns result dict."""
    rng = np.random.default_rng(SEED)
    targets = list(contrast_sizes.keys())
    k = len(targets)
    ses = np.array([se_d(*contrast_sizes[t], TRUE_D) for t in targets])  # (k,)

    # draw d_hat for all reps/targets
    d_hat = rng.normal(loc=TRUE_D, scale=ses, size=(N_REP, k))          # (rep, k)
    z = d_hat / ses
    p = 2.0 * norm.sf(np.abs(z))                                        # two-sided
    correct_dir = d_hat > 0.0

    # nominal per-test
    rej_nominal = (p < ALPHA) & correct_dir

    # Bonferroni
    rej_bonf = (p < (Q / k)) & correct_dir

    # Holm (step-down) per rep
    rej_holm = np.zeros_like(p, dtype=bool)
    order = np.argsort(p, axis=1)
    sorted_p = np.take_along_axis(p, order, axis=1)
    thr = Q / (k - np.arange(k))                                        # k, k-1, ..., 1
    passes = sorted_p <= thr
    # Holm stops at first failure: cumulative AND
    holm_sorted = np.logical_and.accumulate(passes, axis=1)
    # scatter back
    rej_holm_sorted = holm_sorted
    inv = np.argsort(order, axis=1)
    rej_holm = np.take_along_axis(rej_holm_sorted, inv, axis=1) & correct_dir

    # Benjamini-Hochberg (step-up) per rep
    bh_thr = (np.arange(1, k + 1) / k) * Q                              # i/k * q
    below = sorted_p <= bh_thr
    # largest i with sorted_p_i <= (i/k)q ; reject all up to that i
    rej_bh_sorted = np.zeros_like(sorted_p, dtype=bool)
    any_below = below.any(axis=1)
    # index of largest True
    idx_max = np.where(below, np.arange(k), -1).max(axis=1)             # -1 if none
    for j in range(k):
        rej_bh_sorted[:, j] = (idx_max >= j)
    rej_bh_sorted &= any_below[:, None]
    rej_bh = np.take_along_axis(rej_bh_sorted, inv, axis=1) & correct_dir

    def power_by_target(rejmat):
        return {targets[j]: float(rejmat[:, j].mean()) for j in range(k)}

    res = {
        "label": label,
        "seed": SEED,
        "n_rep": N_REP,
        "true_d": TRUE_D,
        "family_level_q": Q,
        "contrast_sizes": {t: list(map(int, contrast_sizes[t])) for t in targets},
        "se_d": {targets[j]: float(ses[j]) for j in range(k)},
        "nominal": power_by_target(rej_nominal),
        "bh_sim": power_by_target(rej_bh),
        "holm": power_by_target(rej_holm),
        "bonferroni": power_by_target(rej_bonf),
        # family-level: all six rejected
        "family_all6_nominal": float((rej_nominal.sum(axis=1) == k).mean()),
        "family_all6_bh": float((rej_bh.sum(axis=1) == k).mean()),
        "family_all6_holm": float((rej_holm.sum(axis=1) == k).mean()),
        "family_all6_bonf": float((rej_bonf.sum(axis=1) == k).mean()),
        "family_mean_rejected_bh": float(rej_bh.sum(axis=1).mean()),
    }
    res["status_bh"] = {t: ("OK>=0.80" if res["bh_sim"][t] >= 0.80 else "UNDERPOWERED")
                        for t in targets}
    return res


def build_contrasts(c1_p53, c1_rest, c2_pm, c2_nsmp):
    return {
        "LHX1(C1)":  (c1_p53, c1_rest),
        "PAX8(C1)":  (c1_p53, c1_rest),
        "GATA2(C2)": (c2_pm, c2_nsmp),
        "HOXA9(C2)": (c2_pm, c2_nsmp),
        "SOX9(C2)":  (c2_pm, c2_nsmp),
        "WT1(C2)":   (c2_pm, c2_nsmp),
    }


if __name__ == "__main__":
    scenarios = {
        # PRIMARY: real resolved pooled sizes from part A
        "REAL_pooled":       build_contrasts(36, 194, 85, 109),
        # sensitivity: spec's EXPECTED pooled sizes
        "EXPECTED_spec":     build_contrasts(50, 183, 78, 95),
        # sensitivity: Confirmatory-only (if run as a standalone confirmatory arm)
        "Confirmatory_only": build_contrasts(16, 119, 53, 66),
        # sensitivity: Discovery-only
        "Discovery_only":    build_contrasts(20, 75, 32, 43),
    }
    out = {}
    for name, cs in scenarios.items():
        out[name] = run(cs, name)

    with open("simulated_power_f1_results.json", "w") as fh:
        json.dump(out, fh, indent=2)

    # human-readable table for the PRIMARY scenario
    def table(res):
        rows = []
        for t in res["nominal"]:
            rows.append((t,
                         res["nominal"][t],
                         res["bh_sim"][t],
                         res["holm"][t],
                         res["bonferroni"][t],
                         res["status_bh"][t]))
        return rows

    print("SEED", SEED, "N_REP", N_REP, "true_d", TRUE_D)
    for name in scenarios:
        r = out[name]
        print("\n=== %s ===" % name)
        print("  sizes:", r["contrast_sizes"]["LHX1(C1)"], "(C1)  ",
              r["contrast_sizes"]["GATA2(C2)"], "(C2)")
        print("  %-11s %8s %8s %8s %8s  %s" %
              ("target", "nominal", "BH-sim", "Holm", "Bonf", "status"))
        for t, nom, bh, hol, bon, st in table(r):
            print("  %-11s %8.3f %8.3f %8.3f %8.3f  %s" % (t, nom, bh, hol, bon, st))
        print("  family all-6 BH power: %.3f | mean #rejected(BH): %.3f" %
              (r["family_all6_bh"], r["family_mean_rejected_bh"]))
