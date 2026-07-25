"""
TASK-028 Freeze B v3 PRIMARY=C -- PHASE 1: RE-VERIFY BLIND TECHNICAL QC (resume point).
Points at sealed_v3. Implements the v3 purity resolution.

PERMITTED (Phase 1): 507 UUID/barcode/patient SET EQUALITY vs frozen cohort; GENCODE v36
60660 universe + identifiers; tpm_unstranded -> log2(TPM+1); all-zero-across-507 rule +
record every all-zero module gene (NO replacement); join subtype + covariate tables;
missingness/dup/non-finite checks; design-matrix RANK + VIF for BOTH the n=506 CPE model
and the n=507 no-purity model; confirm fixed C1-C3 matrix; confirm seeds.

PROHIBITED: subtype score plots; sorting by subtype separation; inspecting subtype
coefficients; module/covariate/threshold changes; biology.

HARD STOPS -> if ANY fires, WRITE the tripped condition and DO NOT proceed to model fit:
VIF>10; rank deficiency; CPE model != 506 OR no-purity != 507; checksum drift;
non-finite scores; any other pre-defined technical stop.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(__file__))
from common_v3 import (MASTER_SEED, GIT_HEAD, substep_seed, sha256_file, read_module_genes,
                       load_cohort, load_star_tpm, ensembl_unversioned,
                       load_cpe_by_patient, load_absolute_by_patient, load_estimate_signatures,
                       B1, ARAN_XLSX, ABSOLUTE_TXT, ESTIMATE_GMT, INTER, RESULTS_V3, LOGS,
                       SUBTYPE_TABLE, COHORT_FILE, SUBTYPES,
                       CPE_NONFINITE_PATIENT, CPE_DROP_REASON, ORIGINALS)
from scoring import ssgsea_matrix

QC = {"phase": 1, "version": "v3", "master_seed": MASTER_SEED, "git_HEAD": GIT_HEAD,
      "stops": [], "checks": {}}


def stop(cond, detail=""):
    QC["stops"].append({"condition": cond, "detail": detail})
    print("HARD STOP CONDITION TRIPPED:", cond, "--", detail, flush=True)


def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


os.makedirs(RESULTS_V3, exist_ok=True)
os.makedirs(INTER, exist_ok=True)
print("=== PHASE 1 (v3) start", utcnow(), "===", flush=True)

# -----------------------------------------------------------------
# 1. Cohort + SET EQUALITY vs frozen subtype table
# -----------------------------------------------------------------
cohort = load_cohort()
QC["checks"]["cohort_n_rows"] = len(cohort)
if len(cohort) != 507:
    stop("cohort_rows_not_507", f"got {len(cohort)}")

n_pat = cohort["patient_barcode"].nunique()
n_uuid = cohort["kept_uuid"].nunique()
n_ali = cohort["kept_aliquot_barcode"].nunique()
QC["checks"]["unique_patients"] = int(n_pat)
QC["checks"]["unique_uuids"] = int(n_uuid)
QC["checks"]["unique_aliquots"] = int(n_ali)
if n_pat != 507 or n_uuid != 507 or n_ali != 507:
    stop("set_not_507_unique", f"pat={n_pat} uuid={n_uuid} aliquot={n_ali}")

sc = cohort["mapped_4way"].value_counts().to_dict()
QC["checks"]["subtype_counts_507"] = {k: int(sc.get(k, 0)) for k in SUBTYPES}
expected_sc = {"POLE": 49, "MMRd": 148, "NSMP": 147, "p53abn": 163}
QC["checks"]["subtype_counts_507_expected"] = expected_sc
if QC["checks"]["subtype_counts_507"] != expected_sc:
    stop("subtype_counts_mismatch", f"{QC['checks']['subtype_counts_507']} vs {expected_sc}")

subt = pd.read_csv(SUBTYPE_TABLE, sep="\t", dtype=str)
QC["checks"]["subtype_table_sha256"] = sha256_file(SUBTYPE_TABLE)
QC["checks"]["cohort_file_sha256"] = sha256_file(COHORT_FILE)
subt_map = dict(zip(subt["patient_barcode"], subt["mapped_4way"]))
mismatch, missing = [], []
for _, r in cohort.iterrows():
    pb = r["patient_barcode"]
    if pb not in subt_map:
        missing.append(pb)
    elif subt_map[pb] != r["mapped_4way"]:
        mismatch.append(pb)
# SET EQUALITY both directions across the 507
cohort_set = set(cohort["patient_barcode"])
subt_507 = set(p for p, s in subt_map.items() if p in cohort_set)
QC["checks"]["set_equality_cohort_vs_subtype_table"] = (
    len(missing) == 0 and len(mismatch) == 0 and cohort_set == subt_507)
QC["checks"]["subtype_label_mismatches"] = len(mismatch)
QC["checks"]["cohort_patients_missing_from_subtype_table"] = len(missing)
if missing:
    stop("cohort_missing_from_subtype_table", f"n={len(missing)}")
if mismatch:
    stop("subtype_label_mismatch", f"n={len(mismatch)}")

# STAR file availability
missing_star = [u for u in cohort["kept_uuid"]
                if not os.path.exists(os.path.join(ORIGINALS, f"{u}__augmented_star_gene_counts.tsv"))]
QC["checks"]["missing_star_files"] = len(missing_star)
if missing_star:
    stop("star_files_missing", f"n={len(missing_star)}")

if QC["stops"]:
    QC["verdict"] = "STOP"
    with open(os.path.join(RESULTS_V3, "phase1_qc_verdict.json"), "w") as f:
        json.dump(QC, f, indent=2)
    print("PHASE 1 STOPPED before matrix load.", flush=True)
    sys.exit(0)

# -----------------------------------------------------------------
# 2. Materialize expression matrix (60660 x 507); verify gene universe identity
# -----------------------------------------------------------------
print("Loading 507 STAR files...", utcnow(), flush=True)
tpm_cols = {}
gene_ids_ref = gene_names_ref = gene_types_ref = None
patient_order = cohort["patient_barcode"].tolist()
uuid_by_pat = dict(zip(cohort["patient_barcode"], cohort["kept_uuid"]))
for i, p in enumerate(patient_order):
    df = load_star_tpm(uuid_by_pat[p])
    gids = df["gene_id"].tolist()
    if gene_ids_ref is None:
        gene_ids_ref = gids
        gene_names_ref = df["gene_name"].tolist()
        gene_types_ref = df["gene_type"].tolist()
    elif gids != gene_ids_ref:
        stop("gene_universe_differs", f"uuid {uuid_by_pat[p]}")
        break
    tpm_cols[p] = df["tpm_unstranded"].to_numpy(dtype=np.float64)
    if (i + 1) % 100 == 0:
        print(f"  loaded {i+1}/507", utcnow(), flush=True)

QC["checks"]["n_genes_gencode_v36"] = len(gene_ids_ref) if gene_ids_ref else 0
QC["checks"]["gene_universe_is_60660"] = (len(gene_ids_ref) == 60660)
QC["checks"]["all_ensg_ids"] = all(g.startswith("ENSG") for g in gene_ids_ref)
if len(gene_ids_ref) != 60660:
    stop("gene_universe_not_60660", f"got {len(gene_ids_ref)}")
if QC["stops"]:
    QC["verdict"] = "STOP"
    with open(os.path.join(RESULTS_V3, "phase1_qc_verdict.json"), "w") as f:
        json.dump(QC, f, indent=2)
    sys.exit(0)

TPM = np.column_stack([tpm_cols[p] for p in patient_order])  # 60660 x 507
QC["checks"]["matrix_shape_genes_x_samples"] = list(TPM.shape)
QC["checks"]["nonfinite_in_raw_tpm"] = int(np.sum(~np.isfinite(TPM)))
if QC["checks"]["nonfinite_in_raw_tpm"] > 0:
    stop("nonfinite_raw_tpm", str(QC["checks"]["nonfinite_in_raw_tpm"]))

ensg_unv = [ensembl_unversioned(g) for g in gene_ids_ref]
QC["checks"]["duplicate_unversioned_ensg"] = int(pd.Series(ensg_unv).duplicated().sum())
QC["checks"]["duplicate_gene_names_in_universe"] = int(pd.Series(gene_names_ref).duplicated().sum())

# -----------------------------------------------------------------
# 3. log2(TPM+1)
# -----------------------------------------------------------------
LOG = np.log2(TPM + 1.0)
QC["checks"]["nonfinite_in_log2"] = int(np.sum(~np.isfinite(LOG)))
if QC["checks"]["nonfinite_in_log2"] > 0:
    stop("nonfinite_log2", str(QC["checks"]["nonfinite_in_log2"]))

# -----------------------------------------------------------------
# 4. all-zero-across-507 rule on MODULE genes (B2.2). Record; NO replacement.
# -----------------------------------------------------------------
sym_to_rows = {}
for idx, nm in enumerate(gene_names_ref):
    sym_to_rows.setdefault(nm, []).append(idx)

allzero_mask = (TPM == 0).all(axis=1)
module_files = {
    "M1_analysis_ready": os.path.join(B1, "M1_analysis_ready.txt"),
    "compact_M1_analysis_ready": os.path.join(B1, "compact_M1_analysis_ready.txt"),
    "M3_primary_analysis_ready": os.path.join(B1, "M3_primary_analysis_ready.txt"),
    "M4_covariate": os.path.join(B1, "M4_covariate.txt"),
}
allzero_removals, unmapped_genes, module_membership = [], [], {}
for mod, path in module_files.items():
    genes = read_module_genes(path)
    module_membership[mod] = genes
    for g in genes:
        rows = sym_to_rows.get(g)
        if rows is None:
            unmapped_genes.append((mod, g))
            continue
        if all(allzero_mask[r] for r in rows):
            allzero_removals.append((mod, g, "all-zero-507"))
QC["checks"]["module_gene_counts"] = {m: len(g) for m, g in module_membership.items()}
QC["checks"]["allzero_removals_n"] = len(allzero_removals)
QC["checks"]["allzero_removals"] = [{"module": m, "gene": g, "reason": r}
                                    for (m, g, r) in allzero_removals]
QC["checks"]["unmapped_module_genes_n"] = len(unmapped_genes)
QC["checks"]["unmapped_module_genes"] = [{"module": m, "gene": g} for (m, g) in unmapped_genes]

# -----------------------------------------------------------------
# 5. Covariates: CPE (n=506), ABSOLUTE (n=502), ESTIMATE composition, M4
# -----------------------------------------------------------------
cpe_by_patient = load_cpe_by_patient()
abs_by_patient = load_absolute_by_patient()

cpe_vals, cpe_missing = [], []
for p in patient_order:
    if p in cpe_by_patient:
        cpe_vals.append(cpe_by_patient[p])
    else:
        cpe_vals.append(np.nan); cpe_missing.append(p)
cpe_vals = np.array(cpe_vals, dtype=float)
QC["checks"]["cpe_finite_coverage_n"] = 507 - len(cpe_missing)
QC["checks"]["cpe_missing_patients"] = [{"patient": p, "subtype": subt_map[p]} for p in cpe_missing]

abs_vals, abs_missing = [], []
for p in patient_order:
    if p in abs_by_patient:
        abs_vals.append(abs_by_patient[p])
    else:
        abs_vals.append(np.nan); abs_missing.append(p)
abs_vals = np.array(abs_vals, dtype=float)
QC["checks"]["absolute_finite_coverage_n"] = 507 - len(abs_missing)
QC["checks"]["absolute_missing_patients"] = [{"patient": p, "subtype": subt_map[p]} for p in abs_missing]

# ESTIMATE composition = stromal + immune ssGSEA (deterministic)
estimate_sig = load_estimate_signatures()
QC["checks"]["estimate_signature_sizes"] = {k: len(v) for k, v in estimate_sig.items()}


def symbol_best_row(sym):
    rows = sym_to_rows.get(sym)
    if not rows:
        return None
    if len(rows) == 1:
        return rows[0]
    means = [LOG[r].mean() for r in rows]
    return rows[int(np.argmax(means))]


def geneset_row_index(genes):
    idx, used = [], set()
    for g in genes:
        r = symbol_best_row(g)
        if r is not None and r not in used:
            idx.append(r); used.add(r)
    return idx


m4_rows = geneset_row_index(module_membership["M4_covariate"])
strom_rows = geneset_row_index(estimate_sig["StromalSignature"])
imm_rows = geneset_row_index(estimate_sig["ImmuneSignature"])
QC["checks"]["m4_genes_scored"] = len(m4_rows)
QC["checks"]["estimate_stromal_scored"] = len(strom_rows)
QC["checks"]["estimate_immune_scored"] = len(imm_rows)

m4_score = ssgsea_matrix(LOG, {"M4": m4_rows})["M4"]
strom_score = ssgsea_matrix(LOG, {"S": strom_rows})["S"]
imm_score = ssgsea_matrix(LOG, {"I": imm_rows})["I"]
composition_score = strom_score + imm_score

for nm, arr in [("m4", m4_score), ("strom", strom_score), ("imm", imm_score),
                ("comp", composition_score)]:
    QC["checks"][f"nonfinite_{nm}"] = int(np.sum(~np.isfinite(arr)))
if QC["checks"]["nonfinite_m4"] or QC["checks"]["nonfinite_comp"]:
    stop("nonfinite_covariate_score", "M4 or composition")

# -----------------------------------------------------------------
# 6. Design-matrix RANK + VIF for BOTH models (v3):
#    (a) n=506 CPE model (drop A0TG); (b) n=507 no-purity model.
#    VIF on covariate design columns only; NO subtype-effect inspection.
# -----------------------------------------------------------------
from statsmodels.stats.outliers_influence import variance_inflation_factor

subt_cat = pd.Categorical(cohort["mapped_4way"].values, categories=SUBTYPES)
subt_dum_full = pd.get_dummies(subt_cat, drop_first=True).astype(float).to_numpy()  # 507 x 3

# mask for CPE complete-case (drop the single non-finite CPE patient)
cpe_finite_mask = np.isfinite(cpe_vals)
QC["checks"]["cpe_complete_case_n"] = int(cpe_finite_mask.sum())
# verify the dropped patient is exactly the expected one
dropped = [patient_order[i] for i in range(507) if not cpe_finite_mask[i]]
QC["checks"]["cpe_dropped_patients"] = dropped
QC["checks"]["cpe_drop_reason"] = CPE_DROP_REASON
if int(cpe_finite_mask.sum()) != 506:
    stop("cpe_model_n_not_506", f"got {int(cpe_finite_mask.sum())}")
if dropped != [CPE_NONFINITE_PATIENT]:
    stop("cpe_dropped_patient_unexpected", f"got {dropped}, expected [{CPE_NONFINITE_PATIENT}]")

# CPE complete-case subtype counts
cc_cpe = {s: 0 for s in SUBTYPES}
for i in range(507):
    if cpe_finite_mask[i]:
        cc_cpe[cohort["mapped_4way"].iloc[i]] += 1
QC["checks"]["cpe_complete_case_subtype_counts"] = cc_cpe
QC["checks"]["cpe_complete_case_subtype_counts_expected"] = {
    "POLE": 49, "MMRd": 148, "NSMP": 146, "p53abn": 163}
if cc_cpe != QC["checks"]["cpe_complete_case_subtype_counts_expected"]:
    stop("cpe_cc_subtype_counts_mismatch", f"{cc_cpe}")

# --- (a) n=506 CPE model design ---
idx506 = np.where(cpe_finite_mask)[0]
X506 = np.column_stack([
    np.ones(506), subt_dum_full[idx506],
    m4_score[idx506], cpe_vals[idx506], composition_score[idx506]])
rank506 = int(np.linalg.matrix_rank(X506))
QC["checks"]["design_n506_ncol"] = X506.shape[1]
QC["checks"]["design_n506_rank"] = rank506
if rank506 < X506.shape[1]:
    stop("rank_deficient_n506", f"rank {rank506} < ncol {X506.shape[1]}")
vif_names506 = ["intercept", "subMMRd", "subNSMP", "subp53abn", "M4", "purity_CPE", "composition"]
vifs506 = {vif_names506[j]: float(variance_inflation_factor(X506, j))
           for j in range(1, X506.shape[1])}
QC["checks"]["VIF_n506_CPE_model"] = vifs506
QC["checks"]["max_VIF_n506"] = max(vifs506.values())

# --- (b) n=507 no-purity model design ---
X507 = np.column_stack([
    np.ones(507), subt_dum_full, m4_score, composition_score])
rank507 = int(np.linalg.matrix_rank(X507))
QC["checks"]["design_n507_ncol"] = X507.shape[1]
QC["checks"]["design_n507_rank"] = rank507
if rank507 < X507.shape[1]:
    stop("rank_deficient_n507", f"rank {rank507} < ncol {X507.shape[1]}")
vif_names507 = ["intercept", "subMMRd", "subNSMP", "subp53abn", "M4", "composition"]
vifs507 = {vif_names507[j]: float(variance_inflation_factor(X507, j))
           for j in range(1, X507.shape[1])}
QC["checks"]["VIF_n507_nopurity_model"] = vifs507
QC["checks"]["max_VIF_n507"] = max(vifs507.values())

# --- ABSOLUTE sensitivity design (n=502) rank+VIF too, for completeness ---
abs_finite_mask = np.isfinite(abs_vals)
QC["checks"]["absolute_complete_case_n"] = int(abs_finite_mask.sum())
abs_dropped = [patient_order[i] for i in range(507) if not abs_finite_mask[i]]
QC["checks"]["absolute_dropped_patients"] = [{"patient": p, "subtype": subt_map[p]} for p in abs_dropped]
cc_abs = {s: 0 for s in SUBTYPES}
for i in range(507):
    if abs_finite_mask[i]:
        cc_abs[cohort["mapped_4way"].iloc[i]] += 1
QC["checks"]["absolute_complete_case_subtype_counts"] = cc_abs
QC["checks"]["absolute_complete_case_subtype_counts_expected"] = {
    "POLE": 48, "MMRd": 148, "NSMP": 144, "p53abn": 162}
if int(abs_finite_mask.sum()) != 502:
    stop("absolute_model_n_not_502", f"got {int(abs_finite_mask.sum())}")
if cc_abs != QC["checks"]["absolute_complete_case_subtype_counts_expected"]:
    stop("absolute_cc_subtype_counts_mismatch", f"{cc_abs}")
idx502 = np.where(abs_finite_mask)[0]
X502 = np.column_stack([
    np.ones(502), subt_dum_full[idx502],
    m4_score[idx502], abs_vals[idx502], composition_score[idx502]])
rank502 = int(np.linalg.matrix_rank(X502))
QC["checks"]["design_n502_ABSOLUTE_rank"] = rank502
QC["checks"]["design_n502_ABSOLUTE_ncol"] = X502.shape[1]
if rank502 < X502.shape[1]:
    stop("rank_deficient_n502", f"rank {rank502} < ncol {X502.shape[1]}")
vifs502 = {vif_names506[j]: float(variance_inflation_factor(X502, j))
           for j in range(1, X502.shape[1])}
QC["checks"]["VIF_n502_ABSOLUTE_model"] = vifs502
QC["checks"]["max_VIF_n502"] = max(vifs502.values())

# global VIF gate: any covariate VIF>10 in any confirmatory/sensitivity design -> STOP
all_max_vif = max(QC["checks"]["max_VIF_n506"], QC["checks"]["max_VIF_n507"],
                  QC["checks"]["max_VIF_n502"])
QC["checks"]["max_VIF_overall"] = all_max_vif
QC["checks"]["VIF_flagged_gt5"] = {
    "n506": {k: v for k, v in vifs506.items() if v > 5},
    "n507": {k: v for k, v in vifs507.items() if v > 5},
    "n502": {k: v for k, v in vifs502.items() if v > 5},
}
if all_max_vif > 10:
    stop("VIF_gt_10", f"max {all_max_vif:.4f} -> model STOP per B2.6")

# -----------------------------------------------------------------
# 7. Confirm fixed C1-C3 matrix + deterministic seeds
# -----------------------------------------------------------------
import itertools
C = {"C1": [-1, -1, -1, 3], "C2": [1, 1, -2, 0], "C3": [1, -1, 0, 0]}
QC["checks"]["contrast_matrix"] = C
cvecs = {k: np.array(v, float) for k, v in C.items()}
orth = {f"{a}.{b}": float(cvecs[a] @ cvecs[b]) for a, b in itertools.combinations(C, 2)}
QC["checks"]["contrast_orthogonality_dotproducts"] = orth
QC["checks"]["contrasts_mutually_orthogonal"] = all(abs(v) < 1e-9 for v in orth.values())
QC["checks"]["seed_substep_examples"] = {
    "bootstrap_M1_C1": substep_seed("bootstrap_M1_C1"),
    "perm_omnibus_M1": substep_seed("perm_omnibus_M1"),
}
QC["checks"]["seeds_deterministic"] = (
    substep_seed("bootstrap_M1_C1") == substep_seed("bootstrap_M1_C1"))

# -----------------------------------------------------------------
# 8. Persist materialized intermediate for Phase 2 (only if PASS)
# -----------------------------------------------------------------
if not QC["stops"]:
    QC["verdict"] = "PASS"
    np.save(os.path.join(INTER, "log2tpm_v3.npy"), LOG)
    with open(os.path.join(INTER, "gene_ids_versioned_v3.json"), "w") as f:
        json.dump(gene_ids_ref, f)
    with open(os.path.join(INTER, "gene_names_v3.json"), "w") as f:
        json.dump(gene_names_ref, f)
    with open(os.path.join(INTER, "patient_order_v3.json"), "w") as f:
        json.dump(patient_order, f)
    cov_df = pd.DataFrame({
        "patient_barcode": patient_order,
        "subtype": cohort["mapped_4way"].tolist(),
        "M4_prolif": m4_score,
        "purity_CPE": cpe_vals,
        "purity_ABSOLUTE": abs_vals,
        "estimate_stromal": strom_score,
        "estimate_immune": imm_score,
        "composition": composition_score,
        "cpe_complete_case": cpe_finite_mask.astype(int),
        "absolute_complete_case": abs_finite_mask.astype(int),
    })
    cov_df.to_csv(os.path.join(INTER, "covariates_v3.tsv"), sep="\t", index=False)
    pd.DataFrame(QC["checks"]["allzero_removals"]).to_csv(
        os.path.join(RESULTS_V3, "allzero_exclusions.tsv"), sep="\t", index=False)
else:
    QC["verdict"] = "STOP"

QC["timestamp_utc"] = utcnow()
with open(os.path.join(RESULTS_V3, "phase1_qc_verdict.json"), "w") as f:
    json.dump(QC, f, indent=2)
print("=== PHASE 1 (v3) verdict:", QC["verdict"], utcnow(), "===", flush=True)
print(json.dumps({k: QC["checks"].get(k) for k in [
    "cohort_n_rows", "subtype_counts_507", "n_genes_gencode_v36",
    "matrix_shape_genes_x_samples", "allzero_removals_n", "unmapped_module_genes_n",
    "cpe_finite_coverage_n", "cpe_complete_case_n", "absolute_finite_coverage_n",
    "absolute_complete_case_n", "design_n506_rank", "design_n506_ncol",
    "design_n507_rank", "design_n507_ncol", "max_VIF_n506", "max_VIF_n507",
    "max_VIF_n502", "max_VIF_overall"]}, indent=2))
