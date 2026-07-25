"""
TASK-028 Freeze B v3 -- PHASE 2a: SCORE all modules (deterministic, no randomness).
Builds the per-sample score matrices required by Phase 2b models:

PRIMARY:
  - M1_analysis_ready ssGSEA (236 minus all-zero) -> M1 score (samples,)
  - 20 M3 CORE_TF SIGNED VIPER/aREA regulons over the 761 consensus-signed edges
    (weight = consensus confidence 1.0; excluded 41 no-consensus; 21 orphans contribute
    nothing) -> M3_<TF> score (samples,)
SENSITIVITY score-sets (scored here, modelled in 2b):
  - compact-M1 ssGSEA (185 minus all-zero)
  - signed-ssGSEA for each of the 20 regulons (primary edge set)
  - S1 regulons: flag-split all 79 both-flagged edges (activation + repression contributions)
  - S2 regulons: exclude all 79 both-flagged edges

Also stores the per-TF z(TFexpr) vector (log2 TPM+1 of the TF gene) for the D_t H4 model.
All scores written to intermediate/scores_v3.npz + a manifest of what was scored.
"""
import os
import sys
import json
import numpy as np
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(__file__))
from common_v3 import (B1, INTER, RESULTS_V3, read_module_genes, sha256_file)
from scoring import ssgsea_matrix, area_regulon, signed_ssgsea_regulon


def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


print("=== PHASE 2a (score) start", utcnow(), "===", flush=True)

LOG = np.load(os.path.join(INTER, "log2tpm_v3.npy"))          # 60660 x 507
gene_names = json.load(open(os.path.join(INTER, "gene_names_v3.json")))
patient_order = json.load(open(os.path.join(INTER, "patient_order_v3.json")))
TPM_allzero = None  # not needed; all-zero handled by dropping genes below
n_genes, n_samples = LOG.shape
assert n_samples == 507

# symbol -> rows; best row = highest mean log2 among duplicate symbols
sym_to_rows = {}
for i, nm in enumerate(gene_names):
    sym_to_rows.setdefault(nm, []).append(i)


def best_row(sym):
    rows = sym_to_rows.get(sym)
    if not rows:
        return None
    if len(rows) == 1:
        return rows[0]
    means = [LOG[r].mean() for r in rows]
    return rows[int(np.argmax(means))]


# all-zero mask (a gene is all-zero if ALL its rows are 0 across 507; log2(0+1)=0)
allzero_sym = set()
for sym, rows in sym_to_rows.items():
    if all((LOG[r] == 0).all() for r in rows):
        allzero_sym.add(sym)

manifest = {"scored_at": utcnow(), "n_samples": n_samples, "modules": {}}

# ---------------- M1 ssGSEA (primary) ----------------
def score_ssgsea_module(module_path, name):
    genes = read_module_genes(module_path)
    kept = [g for g in genes if g not in allzero_sym]
    dropped = [g for g in genes if g in allzero_sym]
    rows, used = [], set()
    for g in kept:
        r = best_row(g)
        if r is not None and r not in used:
            rows.append(r); used.add(r)
    sc = ssgsea_matrix(LOG, {name: rows})[name]
    manifest["modules"][name] = {"method": "ssGSEA", "n_members": len(genes),
                                 "n_after_allzero": len(kept), "n_scored_rows": len(rows),
                                 "allzero_dropped": dropped}
    return sc

m1 = score_ssgsea_module(os.path.join(B1, "M1_analysis_ready.txt"), "M1")
compact_m1 = score_ssgsea_module(os.path.join(B1, "compact_M1_analysis_ready.txt"), "compactM1")

# ---------------- M3 signed VIPER/aREA regulons (primary + S1 + S2) ----------------
# Parse edge ledger; build per-TF regulon lists for primary, S1, S2.
edge_path = os.path.join(B1, "m3_primary_edge_ledger.tsv")
edges = []
with open(edge_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        edges.append(line.rstrip("\n").split("\t"))
hdr = edges[0]; edges = edges[1:]
ix = {c: i for i, c in enumerate(hdr)}

from collections import defaultdict
prim = defaultdict(lambda: {"targets": [], "weights": [], "mor": []})
s1 = defaultdict(lambda: {"targets": [], "weights": [], "mor": []})
s2 = defaultdict(lambda: {"targets": [], "weights": [], "mor": []})

for r in edges:
    tf = r[ix["TF"]]
    tgt = r[ix["target"]]
    w = r[ix["weight"]]
    incl = r[ix["primary_inclusion"]]
    psign = r[ix["primary_sign"]]
    both = r[ix["both_flagged"]]
    s1rep = r[ix["s1_representation"]]
    s2inc = r[ix["s2_inclusion"]]
    # PRIMARY: consensus-resolved edges once, with consensus sign + weight 1.0
    if incl == "Y":
        prim[tf]["targets"].append(tgt)
        prim[tf]["weights"].append(float(w))
        prim[tf]["mor"].append(int(psign))
    # S1: every both-flagged edge (79) -> TWO contributions (+1 and -1) each with weight;
    #     non-both-flagged consensus edges enter once with consensus sign.
    #     s1_representation is 'split-dual' for the 79 both-flagged, 'single' otherwise.
    if s1rep == "split-dual":
        wt = float(w) if w not in ("", None) else 1.0
        s1[tf]["targets"].extend([tgt, tgt])
        s1[tf]["weights"].extend([wt, wt])
        s1[tf]["mor"].extend([1, -1])
    elif s1rep == "single":
        # single = a non-both-flagged consensus edge; enters once with its consensus sign
        if incl == "Y":
            s1[tf]["targets"].append(tgt)
            s1[tf]["weights"].append(float(w))
            s1[tf]["mor"].append(int(psign))
    # S2: all both-flagged excluded; single-flag consensus edges retained by consensus sign
    if s2inc == "Y" and incl == "Y":
        s2[tf]["targets"].append(tgt)
        s2[tf]["weights"].append(float(w))
        s2[tf]["mor"].append(int(psign))

# best-row gene index for aREA
gene_index = {}
for sym in set([t for d in prim.values() for t in d["targets"]] +
               [t for d in s1.values() for t in d["targets"]] +
               [t for d in s2.values() for t in d["targets"]]):
    r = best_row(sym)
    if r is not None:
        gene_index[sym] = r

CORE_TFS = sorted(prim.keys())
manifest["CORE_TFS"] = CORE_TFS
manifest["n_core_tfs"] = len(CORE_TFS)

m3_primary = {}
m3_signed_ssgsea = {}
m3_s1 = {}
m3_s2 = {}
tf_expr_z = {}
tf_expr_log = {}

for tf in CORE_TFS:
    d = prim[tf]
    nes, nused = area_regulon(LOG, gene_index, d["targets"], d["weights"], d["mor"])
    m3_primary[tf] = nes
    # signed-ssGSEA sensitivity on the SAME primary edge set (sign = mor)
    sig, npos, nneg = signed_ssgsea_regulon(LOG, gene_index, d["targets"], d["mor"])
    m3_signed_ssgsea[tf] = sig
    # S1
    d1 = s1[tf]
    nes1, n1 = area_regulon(LOG, gene_index, d1["targets"], d1["weights"], d1["mor"])
    m3_s1[tf] = nes1
    # S2
    d2 = s2[tf]
    nes2, n2 = area_regulon(LOG, gene_index, d2["targets"], d2["weights"], d2["mor"])
    m3_s2[tf] = nes2
    # TF expression (log2 TPM+1) and z across 507 for D_t H4
    r = best_row(tf)
    if r is not None:
        te = LOG[r].astype(float)
        tf_expr_log[tf] = te
        mu = te.mean(); sd = te.std(ddof=1)
        tf_expr_z[tf] = (te - mu) / sd if sd > 0 else np.zeros_like(te)
    else:
        tf_expr_log[tf] = np.full(n_samples, np.nan)
        tf_expr_z[tf] = np.full(n_samples, np.nan)
    manifest["modules"][f"M3_{tf}"] = {
        "method": "signed-VIPER/aREA",
        "primary_targets_used": int(nused),
        "primary_edges_listed": len(d["targets"]),
        "signed_ssgsea_pos": int(npos), "signed_ssgsea_neg": int(nneg),
        "s1_edges_listed": len(d1["targets"]), "s1_targets_used": int(n1),
        "s2_edges_listed": len(d2["targets"]), "s2_targets_used": int(n2),
        "tf_expr_row": (int(best_row(tf)) if best_row(tf) is not None else -1),
    }

# ---------------- pack + save ----------------
save = {"patient_order": np.array(patient_order, dtype=object),
        "M1": m1, "compactM1": compact_m1,
        "CORE_TFS": np.array(CORE_TFS, dtype=object)}
for tf in CORE_TFS:
    save[f"M3primary__{tf}"] = m3_primary[tf]
    save[f"M3signedssgsea__{tf}"] = m3_signed_ssgsea[tf]
    save[f"M3S1__{tf}"] = m3_s1[tf]
    save[f"M3S2__{tf}"] = m3_s2[tf]
    save[f"TFexprZ__{tf}"] = tf_expr_z[tf]
    save[f"TFexprLog__{tf}"] = tf_expr_log[tf]

np.savez(os.path.join(INTER, "scores_v3.npz"), **save)

# non-finite audit on primary scores
nonfinite = {}
for k, v in [("M1", m1), ("compactM1", compact_m1)] + \
            [(f"M3primary__{tf}", m3_primary[tf]) for tf in CORE_TFS]:
    nf = int(np.sum(~np.isfinite(v)))
    if nf:
        nonfinite[k] = nf
manifest["nonfinite_primary_scores"] = nonfinite
manifest["scores_npz"] = os.path.join(INTER, "scores_v3.npz")
with open(os.path.join(RESULTS_V3, "phase2a_score_manifest.json"), "w") as f:
    json.dump(manifest, f, indent=2)

print("=== PHASE 2a done", utcnow(), "n_core_tfs:", len(CORE_TFS),
      "nonfinite_primary:", nonfinite, "===", flush=True)
