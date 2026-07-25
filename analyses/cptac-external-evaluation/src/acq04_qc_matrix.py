#!/usr/bin/env python3
"""
TASK-029 ACQUISITION step 04 -- gene-identity QC + build log2(TPM+1) per stratum.

Reproduces the TASK-028 phase1 QC logic VERBATIM for gene handling:
  - parse each STAR-Counts TSV (skiprows=1; ENSG rows only; tpm_unstranded),
  - verify the gene universe is byte-identical across all 230 files,
  - verify GENCODE v36 / 60660-gene structure (freeze item h),
  - build TWO stratum matrices (Discovery 95, Confirmatory 135), genes x samples,
  - log2(TPM+1),
  - frozen all-zero-across-cohort gene rule PER STRATUM (drop genes all-zero across
    that stratum's analytic cohort; log; NO replacement) -- freeze item i,
  - record every all-zero MODULE gene (M1/M3/M4/ESTIMATE) per stratum.

Scoring uses the sealed modules by GENCODE-v36 gene_name identity (best-row =
highest mean log2 among duplicate symbols), identical to TASK-028.

NO model, NO contrast, NO biological read here. Deterministic. Git HEAD unchanged.
"""
import json
import os
import sys
import hashlib
import numpy as np
import pandas as pd
from datetime import datetime, timezone

BASE = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
WORK = os.environ.get(
    "CPTAC_WORK_DIR",
    os.path.join(BASE, "work"),
)
INTER = os.environ.get(
    "CPTAC_INTERMEDIATE_DIR",
    os.path.join(WORK, "intermediate"),
)
RESULTS = os.environ.get(
    "CPTAC_RESULTS_DIR",
    os.path.join(BASE, "results"),
)
RAWDIR = os.environ.get(
    "CPTAC_STAR_COUNTS_DIR",
    os.path.join(WORK, "rna-star-counts"),
)
DISCOVERY_ROOT = os.environ.get(
    "TCGA_DISCOVERY_ROOT",
    os.path.join(
        os.path.dirname(BASE),
        "tcga-primary-discovery",
    ),
)
B1 = os.path.join(DISCOVERY_ROOT, "definitions")
ESTIMATE_GMT = os.environ.get(
    "ESTIMATE_GMT",
    os.path.join(
        DISCOVERY_ROOT,
        "external-inputs",
        "ESTIMATE_SI_geneset.gmt",
    ),
)

os.makedirs(RESULTS, exist_ok=True)


def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_module_genes(path):
    g = []
    for line in open(path):
        s = line.strip()
        if s and not s.startswith("#"):
            g.append(s)
    return g


lg = json.load(open(os.path.join(INTER, "acq02_join_ledger.json")))
pick = lg["pick"]
verify = json.load(open(os.path.join(INTER, "acq03_download_verify.json")))["per_case"]

QC = {"phase": "acq04_qc", "git_HEAD": "83503bad47b60193598b2b9ebe819c22c83e8ac1",
      "stops": [], "checks": {}, "started": utcnow()}


def stop(cond, detail=""):
    QC["stops"].append({"condition": cond, "detail": detail})
    print("HARD STOP:", cond, "--", detail, flush=True)


def load_star(path):
    df = pd.read_csv(path, sep="\t", comment=None, skiprows=1,
                     dtype={"gene_id": str, "gene_name": str, "gene_type": str})
    df = df[df["gene_id"].str.startswith("ENSG")].copy()
    return df


strata = {"Discovery": [], "Confirmatory": []}
for c in sorted(pick):
    strata[pick[c]["stratum"]].append(c)
QC["checks"]["n_discovery"] = len(strata["Discovery"])
QC["checks"]["n_confirmatory"] = len(strata["Confirmatory"])

gene_ids_ref = gene_names_ref = gene_types_ref = None
strat_matrices = {}
strat_order = {}
for strat in ["Discovery", "Confirmatory"]:
    cases = strata[strat]
    cols = {}
    for c in cases:
        df = load_star(verify[c]["path"])
        gids = df["gene_id"].tolist()
        if gene_ids_ref is None:
            gene_ids_ref = gids
            gene_names_ref = df["gene_name"].tolist()
            gene_types_ref = df["gene_type"].tolist()
        elif gids != gene_ids_ref:
            stop("gene_universe_differs", f"case {c}")
        cols[c] = df["tpm_unstranded"].to_numpy(dtype=np.float64)
    TPM = np.column_stack([cols[c] for c in cases])  # genes x samples
    strat_matrices[strat] = TPM
    strat_order[strat] = cases
    QC["checks"][f"{strat}_matrix_shape"] = list(TPM.shape)
    QC["checks"][f"{strat}_nonfinite_raw_tpm"] = int(np.sum(~np.isfinite(TPM)))

QC["checks"]["n_genes_gencode_v36"] = len(gene_ids_ref) if gene_ids_ref else 0
QC["checks"]["gene_universe_is_60660"] = (len(gene_ids_ref) == 60660)
QC["checks"]["all_ensg_ids"] = all(g.startswith("ENSG") for g in gene_ids_ref)
if len(gene_ids_ref) != 60660:
    stop("gene_universe_not_60660", f"got {len(gene_ids_ref)}")

# gene_id ordering identity hash (should equal the TASK-028 universe if same release+workflow)
gid_hash = hashlib.sha256("\n".join(gene_ids_ref).encode()).hexdigest()
gname_hash = hashlib.sha256("\n".join(gene_names_ref).encode()).hexdigest()
QC["checks"]["gene_id_order_sha256"] = gid_hash
QC["checks"]["gene_name_order_sha256"] = gname_hash

# compare to TASK-028 frozen gene universe (if intermediate present)
t28_gids_path = os.environ.get(
    "TCGA_GENE_IDS_VERSIONED",
    os.path.join(
        DISCOVERY_ROOT,
        "external-inputs",
        "gene_ids_versioned_v3.json",
    ),
)
if os.path.exists(t28_gids_path):
    t28_gids = json.load(open(t28_gids_path))
    QC["checks"]["task028_gene_universe_identical"] = (t28_gids == gene_ids_ref)
    QC["checks"]["task028_n_genes"] = len(t28_gids)

QC["checks"]["nonfinite_raw_tpm_total"] = (
    QC["checks"]["Discovery_nonfinite_raw_tpm"] + QC["checks"]["Confirmatory_nonfinite_raw_tpm"])
if QC["checks"]["nonfinite_raw_tpm_total"] > 0:
    stop("nonfinite_raw_tpm", str(QC["checks"]["nonfinite_raw_tpm_total"]))

# ---- log2(TPM+1) + per-stratum all-zero rule ----
sym_to_rows = {}
for idx, nm in enumerate(gene_names_ref):
    sym_to_rows.setdefault(nm, []).append(idx)

module_files = {
    "M1_analysis_ready": os.path.join(B1, "M1_analysis_ready.txt"),
    "compact_M1_analysis_ready": os.path.join(B1, "compact_M1_analysis_ready.txt"),
    "M3_primary_analysis_ready": os.path.join(B1, "M3_primary_analysis_ready.txt"),
    "M4_covariate": os.path.join(B1, "M4_covariate.txt"),
}
# ESTIMATE signatures
estimate_sig = {}
for line in open(ESTIMATE_GMT):
    p = line.rstrip("\n").split("\t")
    estimate_sig[p[0]] = p[2:]

allzero_ledger = {"Discovery": {}, "Confirmatory": {}}
LOGmats = {}
for strat in ["Discovery", "Confirmatory"]:
    TPM = strat_matrices[strat]
    LOG = np.log2(TPM + 1.0)
    LOGmats[strat] = LOG
    if int(np.sum(~np.isfinite(LOG))) > 0:
        stop("nonfinite_log2", strat)
    allzero_mask = (TPM == 0).all(axis=1)
    QC["checks"][f"{strat}_allzero_genes_total"] = int(allzero_mask.sum())
    # per-module all-zero removals + unmapped
    mod_removals = {}
    mod_unmapped = {}
    for mod, path in module_files.items():
        genes = read_module_genes(path)
        rem, unm = [], []
        for g in genes:
            rows = sym_to_rows.get(g)
            if rows is None:
                unm.append(g)
            elif all(allzero_mask[r] for r in rows):
                rem.append(g)
        mod_removals[mod] = rem
        mod_unmapped[mod] = unm
    # ESTIMATE all-zero
    est_removals = {}
    for signm, genes in estimate_sig.items():
        rem = [g for g in genes if sym_to_rows.get(g) and all(allzero_mask[r] for r in sym_to_rows[g])]
        est_removals[signm] = rem
    allzero_ledger[strat] = {"module_allzero_removals": mod_removals,
                             "module_unmapped": mod_unmapped,
                             "estimate_allzero_removals": est_removals}
    QC["checks"][f"{strat}_module_allzero_counts"] = {m: len(v) for m, v in mod_removals.items()}
    QC["checks"][f"{strat}_module_unmapped_counts"] = {m: len(v) for m, v in mod_unmapped.items()}

# ---- persist matrices + gene info + ledgers ----
np.save(os.path.join(INTER, "log2tpm_Discovery.npy"), LOGmats["Discovery"])
np.save(os.path.join(INTER, "log2tpm_Confirmatory.npy"), LOGmats["Confirmatory"])
json.dump(gene_ids_ref, open(os.path.join(INTER, "gene_ids_versioned.json"), "w"))
json.dump(gene_names_ref, open(os.path.join(INTER, "gene_names.json"), "w"))
json.dump(gene_types_ref, open(os.path.join(INTER, "gene_types.json"), "w"))
json.dump({s: strat_order[s] for s in strat_order},
          open(os.path.join(INTER, "case_order_by_stratum.json"), "w"))
json.dump(allzero_ledger, open(os.path.join(RESULTS, "gene_flow_ledger.json"), "w"), indent=2)

QC["verdict"] = "STOP" if QC["stops"] else "PASS"
QC["finished"] = utcnow()
json.dump(QC, open(os.path.join(RESULTS, "acq04_qc_verdict.json"), "w"), indent=2)
print("QC verdict:", QC["verdict"], utcnow(), flush=True)
print(json.dumps({k: QC["checks"].get(k) for k in [
    "n_discovery", "n_confirmatory", "Discovery_matrix_shape", "Confirmatory_matrix_shape",
    "n_genes_gencode_v36", "gene_universe_is_60660", "task028_gene_universe_identical",
    "Discovery_allzero_genes_total", "Confirmatory_allzero_genes_total",
    "Discovery_module_allzero_counts", "Confirmatory_module_allzero_counts",
    "Discovery_module_unmapped_counts", "Confirmatory_module_unmapped_counts"]}, indent=2))
