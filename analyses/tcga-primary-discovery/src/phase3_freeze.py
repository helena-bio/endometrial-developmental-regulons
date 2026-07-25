"""
TASK-028 Freeze B v3 -- PHASE 3: RESULTS FREEZE (no biological narrative).
Writes flat result tables (TSV) + ledgers + diagnostics under results_v3/, then the
RESULT_SEAL.sha256 + reproducibility manifest. Claim firewall: result language only.
"""
import os
import sys
import json
import subprocess
import numpy as np
import pandas as pd
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(__file__))
from common_v3 import (MASTER_SEED, GIT_HEAD, RESULTS_V3, INTER, SUBTYPES, sha256_file,
                       ARAN_XLSX, ABSOLUTE_TXT, ESTIMATE_GMT, COHORT_FILE, SUBTYPE_TABLE)


def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


R = json.load(open(os.path.join(RESULTS_V3, "phase2b_models.json")))
QC = json.load(open(os.path.join(RESULTS_V3, "phase1_qc_verdict.json")))
SM = json.load(open(os.path.join(RESULTS_V3, "phase2a_score_manifest.json")))
CORE = SM["CORE_TFS"]
prim = R["primary"]; sens = R["sensitivities"]; dt = R["dt"]
order = ["M1"] + [f"M3_{tf}" for tf in CORE]


def nn(x):
    return "" if x is None else x


# ------------------------------------------------------------------
# (1) PRIMARY result table: per module x {omnibus, C1,C2,C3}
# ------------------------------------------------------------------
rows = []
for mk in order:
    m = prim[mk]
    om = m["omnibus"]
    base = dict(module=mk, n=m["n"], sigma_resid=m["sigma_resid"],
                omnibus_F=om["F"], omnibus_df1=om["df1"], omnibus_df2=om["df2"],
                omnibus_p_F=om["p_F"], omnibus_BH_q_F1=om["BH_q_F1"],
                omnibus_F1_gate_pass=om["F1_gate_pass"],
                kruskal_H=m["kruskal"]["H"], kruskal_p=m["kruskal"]["p"])
    for c in ["C1", "C2", "C3"]:
        cc = m["contrasts"][c]; b = m["bootstrap"][c]
        cat = m.get("b3_category", {}).get(c, {})
        rows.append({**base, "contrast": c,
                     "estimate": cc["estimate"], "d": cc["d"],
                     "est_ci_lo": b["est_ci_lo"], "est_ci_hi": b["est_ci_hi"],
                     "d_ci_lo": b["d_ci_lo"], "d_ci_hi": b["d_ci_hi"],
                     "perm_p": cc.get("perm_p"), "BH_q_F2": cc.get("BH_q_F2"),
                     "b3_category": cat.get("category"),
                     "category_reason": cat.get("reason")})
pd.DataFrame(rows).to_csv(os.path.join(RESULTS_V3, "primary_results.tsv"), sep="\t", index=False)

# equivalence (F4) table -- C1 only, F1-gated modules
eq_rows = []
for mk in order:
    m = prim[mk]
    if "equivalence_C1" in m:
        e = m["equivalence_C1"]
        eq_rows.append(dict(module=mk, contrast="C1", d=e["d"],
                            d_ci_lo=e["d_ci_lo"], d_ci_hi=e["d_ci_hi"],
                            SESOI=e.get("SESOI"), tost_pass=e["tost_pass"]))
pd.DataFrame(eq_rows).to_csv(os.path.join(RESULTS_V3, "equivalence_C1_F4.tsv"), sep="\t", index=False)

# ------------------------------------------------------------------
# (1b) SENSITIVITY result tables (one per config) + concordance table
# ------------------------------------------------------------------
for cfg_name, cfg in sens.items():
    srows = []
    for mk, m in cfg.items():
        if "contrasts" not in m:
            srows.append(dict(module=mk, error=m.get("error")))
            continue
        om = m["omnibus"]
        for c in ["C1", "C2", "C3"]:
            cc = m["contrasts"][c]; b = m["bootstrap"][c]
            cat = m.get("b3_category", {}).get(c, {})
            srows.append(dict(module=mk, n=m["n"], contrast=c,
                              omnibus_F=om["F"], omnibus_p_F=om["p_F"],
                              omnibus_BH_q_F1=om["BH_q_F1"], F1_gate_pass=om["F1_gate_pass"],
                              estimate=cc["estimate"], d=cc["d"],
                              d_ci_lo=b["d_ci_lo"], d_ci_hi=b["d_ci_hi"],
                              perm_p=cc.get("perm_p"), BH_q_F2=cc.get("BH_q_F2"),
                              b3_category=cat.get("category")))
    pd.DataFrame(srows).to_csv(
        os.path.join(RESULTS_V3, f"sensitivity_{cfg_name}.tsv"), sep="\t", index=False)

# concordance vs primary (direction / magnitude-delta / category / gate)
def sign(x):
    return 0 if x is None or abs(x) < 1e-9 else (1 if x > 0 else -1)

conc = []
which = {"SENS_nopurity": order, "SENS_absolute": order, "SENS_compactM1": ["M1"],
         "SENS_signedssgsea": [f"M3_{t}" for t in CORE],
         "SENS_S1": [f"M3_{t}" for t in CORE], "SENS_S2": [f"M3_{t}" for t in CORE]}
for cfg_name, mods in which.items():
    cfg = sens[cfg_name]
    for mk in mods:
        if mk not in cfg or "contrasts" not in cfg[mk]:
            continue
        pg = prim[mk]["omnibus"]["F1_gate_pass"]; sg = cfg[mk]["omnibus"]["F1_gate_pass"]
        for c in ["C1", "C2", "C3"]:
            pdv = prim[mk]["contrasts"][c]["d"]; sdv = cfg[mk]["contrasts"][c]["d"]
            pc = prim[mk].get("b3_category", {}).get(c, {}).get("category")
            scv = cfg[mk].get("b3_category", {}).get(c, {}).get("category")
            conc.append(dict(config=cfg_name, module=mk, contrast=c,
                             primary_d=pdv, sens_d=sdv, delta_d=(sdv - pdv),
                             direction_agree=(sign(pdv) == sign(sdv)),
                             primary_category=pc, sens_category=scv,
                             category_match=(pc == scv),
                             primary_gate=pg, sens_gate=sg, gate_match=(pg == sg)))
pd.DataFrame(conc).to_csv(os.path.join(RESULTS_V3, "sensitivity_concordance.tsv"),
                          sep="\t", index=False)

# ------------------------------------------------------------------
# (1c) D_t (H4) discordance table -- includes negatives
# ------------------------------------------------------------------
drows = []
for tf in CORE:
    d = dt[tf]
    for c in ["C1", "C2", "C3"]:
        D = d["D"][c]; eX = d["e_X"][c]; dec = d["declared"][c]
        drows.append(dict(TF=tf, contrast=c, D=D["D"], D_ci_lo=D["ci_lo"], D_ci_hi=D["ci_hi"],
                          D_ci_excludes_0=D["ci_excludes_0"], perm_p=D["perm_p"],
                          BH_q_F3=dec["BH_q_F3"], e_X=eX["e_X"], e_X_ci_lo=eX["ci_lo"],
                          e_X_ci_hi=eX["ci_hi"], e_X_ci_excludes_0=eX["ci_excludes_0"],
                          cond_i_absD_ge_0p5=dec["cond_i_absD_ge_0.5"],
                          cond_ii_bootCI_excl0=dec["cond_ii_bootCI_excl0"],
                          cond_iii_permF3_q_le_0p05=dec["cond_iii_permF3_q_le_0.05"],
                          cond_iv_direction=dec["cond_iv_direction"],
                          verdict=dec["verdict"]))
pd.DataFrame(drows).to_csv(os.path.join(RESULTS_V3, "discordance_Dt_H4.tsv"), sep="\t", index=False)

# ------------------------------------------------------------------
# (2) sample-flow ledger (per-estimator excluded sets) + gene-flow ledger
# ------------------------------------------------------------------
c = QC["checks"]
sample_flow = {
    "frozen_cohort_n": 507,
    "subtype_counts_507": c["subtype_counts_507"],
    "CPE_primary": {"n": c["cpe_complete_case_n"],
                    "subtype_counts": c["cpe_complete_case_subtype_counts"],
                    "excluded": c["cpe_dropped_patients"],
                    "exclusion_reason": c["cpe_drop_reason"]},
    "no_purity_sensitivity": {"n": 507, "subtype_counts": c["subtype_counts_507"],
                              "excluded": []},
    "ABSOLUTE_sensitivity": {"n": c["absolute_complete_case_n"],
                             "subtype_counts": c["absolute_complete_case_subtype_counts"],
                             "excluded": c["absolute_dropped_patients"]},
}
with open(os.path.join(RESULTS_V3, "sample_flow_ledger.json"), "w") as f:
    json.dump(sample_flow, f, indent=2)

gene_flow = {
    "gene_universe_gencode_v36": c["n_genes_gencode_v36"],
    "module_gene_counts_before_allzero": c["module_gene_counts"],
    "allzero_removals": c["allzero_removals"],
    "unmapped_module_genes": c["unmapped_module_genes"],
    "M1_scored_rows": SM["modules"]["M1"]["n_scored_rows"],
    "compactM1_scored_rows": SM["modules"]["compactM1"]["n_scored_rows"],
    "M3_regulon_edge_accounting": {
        "primary_included_edges": 761, "excluded_no_consensus_edges": 41,
        "both_flagged_edges": 79, "single_flag_edges": 723,
        "orphan_targets": 21, "core_tfs_scorable": len(CORE)},
    "M3_per_regulon": {f"M3_{tf}": SM["modules"][f"M3_{tf}"] for tf in CORE},
}
with open(os.path.join(RESULTS_V3, "gene_flow_ledger.json"), "w") as f:
    json.dump(gene_flow, f, indent=2)

# ------------------------------------------------------------------
# (4) model diagnostics (VIF per model, bootstrap/permutation diagnostics)
# ------------------------------------------------------------------
boot_fail = {}
for mk in order:
    boot_fail[mk] = prim[mk]["bootstrap"].get("_n_boot_failed", 0)
diag = {
    "VIF_n506_CPE": c["VIF_n506_CPE_model"], "max_VIF_n506": c["max_VIF_n506"],
    "VIF_n507_nopurity": c["VIF_n507_nopurity_model"], "max_VIF_n507": c["max_VIF_n507"],
    "VIF_n502_ABSOLUTE": c["VIF_n502_ABSOLUTE_model"], "max_VIF_n502": c["max_VIF_n502"],
    "design_ranks": {"n506": [c["design_n506_rank"], c["design_n506_ncol"]],
                     "n507": [c["design_n507_rank"], c["design_n507_ncol"]],
                     "n502": [c["design_n502_ABSOLUTE_rank"], c["design_n502_ABSOLUTE_ncol"]]},
    "n_boot": R["n_boot"], "n_perm": R["n_perm"],
    "primary_bootstrap_failed_resamples": boot_fail,
    "contrast_orthogonality": c["contrast_orthogonality_dotproducts"],
    "SESOI": R["SESOI"], "H1_floor": R["H1_floor"],
}
with open(os.path.join(RESULTS_V3, "model_diagnostics.json"), "w") as f:
    json.dump(diag, f, indent=2)

# ------------------------------------------------------------------
# (6) software/package versions + runtime env
# ------------------------------------------------------------------
import platform
pipfreeze = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                           capture_output=True, text=True).stdout
env = {
    "python": sys.version, "platform": platform.platform(),
    "numpy": np.__version__, "pandas": pd.__version__,
    "pip_freeze": pipfreeze.strip().splitlines(),
    "timestamp_utc": utcnow(),
}
with open(os.path.join(RESULTS_V3, "software_env.json"), "w") as f:
    json.dump(env, f, indent=2)

# ------------------------------------------------------------------
# (7) reproducibility manifest FIRST, then RESULT_SEAL.sha256.
# The seal covers every RESULT file + intermediate + script. It EXCLUDES two files by
# construction: RESULT_SEAL.sha256 (cannot hash itself) and REPRODUCIBILITY_MANIFEST.json
# (the provenance record that documents the seal; it is self-referential -- it lists the
# sealed set and a timestamp written after hashing). Both are named in the seal header.
# ------------------------------------------------------------------
git_head_now = subprocess.run(["git", "-C", os.environ.get("ANALYSIS_GIT_ROOT", ROOT),
                               "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()

# reproducibility manifest (written BEFORE the seal so the sealed set is stable)
manifest = {
    "task": "TASK-028 Freeze B v3 PRIMARY=C locked execution (resumed)",
    "git_HEAD_start": GIT_HEAD, "git_HEAD_at_seal": git_head_now,
    "git_unchanged": GIT_HEAD == git_head_now,
    "master_seed": MASTER_SEED,
    "sub_seed_rule": "sha256('20260713:<step_id>')[:8] -> int; step ids logged in each config/module",
    "data_sources": {
        "cohort": {"path": COHORT_FILE, "sha256": sha256_file(COHORT_FILE)},
        "subtype_table": {"path": SUBTYPE_TABLE, "sha256": sha256_file(SUBTYPE_TABLE)},
        "aran_cpe_xlsx": {"path": ARAN_XLSX, "sha256": sha256_file(ARAN_XLSX),
                          "column": "CPE", "complete_case_n": 506},
        "pancanatlas_absolute": {"path": ABSOLUTE_TXT, "sha256": sha256_file(ABSOLUTE_TXT),
                                 "column": "purity", "complete_case_n": 502},
        "estimate_gmt": {"path": ESTIMATE_GMT, "sha256": sha256_file(ESTIMATE_GMT)},
        "star_originals_dir": ORIGINALS,
    },
    "sealed_packet": "sealed_v3 (all 23 input SHAs verified vs SEAL_MANIFEST.sha256)",
    "n_boot": R["n_boot"], "n_perm": R["n_perm"],
    "package_versions": {"numpy": np.__version__, "pandas": pd.__version__,
                         "python": sys.version.split()[0]},
    "commands": [
        "python3 src/phase1_qc_v3.py",
        "python3 src/phase2a_score.py",
        "python3 src/phase2b_models.py",
        "python3 src/phase3_freeze.py",
    ],
    "timestamp_utc": utcnow(),
}
# seal EXCLUDES: RESULT_SEAL.sha256 (self) + REPRODUCIBILITY_MANIFEST.json (self-referential)
SEAL_EXCLUDE = {"RESULT_SEAL.sha256", "REPRODUCIBILITY_MANIFEST.json"}
result_files = sorted([f for f in os.listdir(RESULTS_V3)
                       if os.path.isfile(os.path.join(RESULTS_V3, f))
                       and f not in SEAL_EXCLUDE])
manifest["result_files_sealed"] = result_files
manifest["seal_excludes"] = sorted(SEAL_EXCLUDE)
with open(os.path.join(RESULTS_V3, "REPRODUCIBILITY_MANIFEST.json"), "w") as f:
    json.dump(manifest, f, indent=2)

# now compute the seal over the stable result set + intermediates + scripts
seal_lines = []
for f in result_files:
    seal_lines.append(f"{sha256_file(os.path.join(RESULTS_V3, f))}  {f}")
extra = {
    "scores_npz": os.path.join(INTER, "scores_v3.npz"),
    "covariates_tsv": os.path.join(INTER, "covariates_v3.tsv"),
}
for name, path in extra.items():
    seal_lines.append(f"{sha256_file(path)}  ../intermediate/{os.path.basename(path)}")
scripts_dir = os.path.dirname(__file__)
for sf in sorted(os.listdir(scripts_dir)):
    if sf.endswith(".py"):
        seal_lines.append(f"{sha256_file(os.path.join(scripts_dir, sf))}  ../scripts_v3/{sf}")

seal_header = [
    "# TASK-028 FREEZE B v3 -- RESULT_SEAL (SHA-256 of every result + intermediate + script)",
    "# EXCLUDES (self-referential, named not hashed): RESULT_SEAL.sha256, REPRODUCIBILITY_MANIFEST.json",
    f"# git_HEAD_start: {GIT_HEAD}",
    f"# git_HEAD_at_seal: {git_head_now}",
    f"# git_unchanged: {GIT_HEAD == git_head_now}",
    f"# master_seed: {MASTER_SEED}",
    f"# seal_timestamp_utc: {utcnow()}",
    f"# n_boot: {R['n_boot']}  n_perm: {R['n_perm']}  SESOI: {R['SESOI']}  H1_floor: {R['H1_floor']}",
]
with open(os.path.join(RESULTS_V3, "RESULT_SEAL.sha256"), "w") as f:
    f.write("\n".join(seal_header + seal_lines) + "\n")

print("PHASE 3 freeze complete", utcnow())
print("git unchanged:", GIT_HEAD == git_head_now, git_head_now)
print("result files sealed:", len(result_files))
