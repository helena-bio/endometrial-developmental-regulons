#!/usr/bin/env python3
"""
TASK-029 (d) PURITY: complete-case BH-of-6 power recompute.
DATA-FREE. No expression, no scoring, no purity VALUE enters the power model --
only the complete-case ARM SIZES (integer counts) after purity-missingness.
Reuses the IDENTICAL power model as simulate_power_f1.py (same SE(d), same BH-of-6,
same seed, same N_REP, same true_d floor). ONLY the arm sizes change.

Purpose: does the complete-case restriction (dropping patients with no finite
published purity value) drop C2 below the frozen evaluability criterion
(simulated BH-of-6 power >= 0.80)?

Complete-case counts (from covariate-QC read of published CPTAC-UCEC purity):
  Discovery PURITY_CANCER (proteogenomic decomposition, cBioPortal ucec_cptac_2020)
  per-subtype FINITE-value missingness (row-present but non-finite):
    POLE 0/7, MMRd 1/25, NSMP 0/43, p53abn 3/20  -> Discovery complete-case
    POLE 7, MMRd 24, NSMP 43, p53abn 17.
  Confirmatory (Dou 2023) per-case purity is behind the PMC gate at recon
  (missingness UNKNOWN) -> two bounding scenarios:
    A: Confirmatory FULL (only Discovery restricted)  -> C1 33/193, C2 84/109
    B: Confirmatory restricted at the Discovery per-subtype rate -> C1 31/191, C2 82/109
  Baseline (no purity restriction, from independence check): C1 36/194, C2 85/109.
"""
import numpy as np
from scipy.stats import norm
import json

SEED = 20260714           # pinned, IDENTICAL to simulate_power_f1.py
N_REP = 200_000
TRUE_D = 0.5
Q = 0.05
ALPHA = 0.05


def se_d(n1, n2, d):
    n1 = float(n1); n2 = float(n2)
    return np.sqrt((n1 + n2) / (n1 * n2) + d * d / (2.0 * (n1 + n2)))


def run(contrast_sizes, label):
    rng = np.random.default_rng(SEED)
    targets = list(contrast_sizes.keys())
    k = len(targets)
    ses = np.array([se_d(*contrast_sizes[t], TRUE_D) for t in targets])
    d_hat = rng.normal(loc=TRUE_D, scale=ses, size=(N_REP, k))
    z = d_hat / ses
    p = 2.0 * norm.sf(np.abs(z))
    correct_dir = d_hat > 0.0

    rej_nominal = (p < ALPHA) & correct_dir
    rej_bonf = (p < (Q / k)) & correct_dir

    order = np.argsort(p, axis=1)
    sorted_p = np.take_along_axis(p, order, axis=1)
    inv = np.argsort(order, axis=1)

    # Holm
    thr = Q / (k - np.arange(k))
    passes = sorted_p <= thr
    holm_sorted = np.logical_and.accumulate(passes, axis=1)
    rej_holm = np.take_along_axis(holm_sorted, inv, axis=1) & correct_dir

    # BH
    bh_thr = (np.arange(1, k + 1) / k) * Q
    below = sorted_p <= bh_thr
    any_below = below.any(axis=1)
    idx_max = np.where(below, np.arange(k), -1).max(axis=1)
    rej_bh_sorted = np.zeros_like(sorted_p, dtype=bool)
    for j in range(k):
        rej_bh_sorted[:, j] = (idx_max >= j)
    rej_bh_sorted &= any_below[:, None]
    rej_bh = np.take_along_axis(rej_bh_sorted, inv, axis=1) & correct_dir

    def pbt(m):
        return {targets[j]: float(m[:, j].mean()) for j in range(k)}

    res = {
        "label": label, "seed": SEED, "n_rep": N_REP, "true_d": TRUE_D,
        "family_level_q": Q,
        "contrast_sizes": {t: list(map(int, contrast_sizes[t])) for t in targets},
        "se_d": {targets[j]: float(ses[j]) for j in range(k)},
        "nominal": pbt(rej_nominal), "bh_sim": pbt(rej_bh),
        "holm": pbt(rej_holm), "bonferroni": pbt(rej_bonf),
        "family_all6_bh": float((rej_bh.sum(axis=1) == k).mean()),
        "family_mean_rejected_bh": float(rej_bh.sum(axis=1).mean()),
    }
    res["status_bh"] = {t: ("OK>=0.80" if res["bh_sim"][t] >= 0.80 else "UNDERPOWERED")
                        for t in targets}
    return res


def build(c1_p53, c1_rest, c2_pm, c2_nsmp):
    return {
        "LHX1(C1)": (c1_p53, c1_rest), "PAX8(C1)": (c1_p53, c1_rest),
        "GATA2(C2)": (c2_pm, c2_nsmp), "HOXA9(C2)": (c2_pm, c2_nsmp),
        "SOX9(C2)": (c2_pm, c2_nsmp), "WT1(C2)": (c2_pm, c2_nsmp),
    }


if __name__ == "__main__":
    scenarios = {
        # baseline, no purity restriction (matches simulate_power_f1.py REAL_pooled)
        "BASELINE_no_purity":       build(36, 194, 85, 109),
        # complete-case A: Discovery restricted, Confirmatory full
        "CC_A_ConfFull":            build(33, 193, 84, 109),
        # complete-case B: both restricted at Discovery per-subtype rate
        "CC_B_ConfSameRate":        build(31, 191, 82, 109),
        # Discovery-only complete-case (single-stratum floor)
        "DISCOVERY_CC_only":        build(17, 74, 31, 43),
    }
    out = {}
    for name, cs in scenarios.items():
        out[name] = run(cs, name)
    with open("simulate_power_f1_completecase_results.json", "w") as fh:
        json.dump(out, fh, indent=2)

    print("SEED", SEED, "N_REP", N_REP, "true_d", TRUE_D,
          "| evaluability criterion: BH-of-6 power >= 0.80")
    for name in scenarios:
        r = out[name]
        c1 = r["contrast_sizes"]["LHX1(C1)"]
        c2 = r["contrast_sizes"]["GATA2(C2)"]
        print("\n=== %s ===  C1 %s  C2 %s" % (name, c1, c2))
        print("  %-11s %8s %8s  %s" % ("target", "nominal", "BH-sim", "status"))
        for t in r["bh_sim"]:
            print("  %-11s %8.3f %8.3f  %s" %
                  (t, r["nominal"][t], r["bh_sim"][t], r["status_bh"][t]))
        # C2 min BH power (the evaluability-decisive number)
        c2t = [t for t in r["bh_sim"] if "(C2)" in t]
        c1t = [t for t in r["bh_sim"] if "(C1)" in t]
        c2min = min(r["bh_sim"][t] for t in c2t)
        c1min = min(r["bh_sim"][t] for t in c1t)
        print("  C2 min BH power = %.3f (%s 0.80) | C1 min BH power = %.3f (%s 0.80)" %
              (c2min, ">=" if c2min >= 0.80 else "<", c1min, ">=" if c1min >= 0.80 else "<"))
