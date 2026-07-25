#!/usr/bin/env python3
# TASK-028 Freeze B -- FINAL SEALED B1 module builder.
# Tumour-blind, deterministic. Applies set operations + relabeling + edge-ledger
# transforms to ALREADY-DERIVED frozen candidate artifacts. NO tumour value read.
# Master SEED 20260713 (no randomness is used; recorded for provenance only).
#
# Source artifacts are supplied through TCGA_DEFINITION_SOURCE_DIR; outputs go to definitions/.
# ASCII only. Deterministic: every run recomputes identical SHA-256.

import csv
import hashlib
import os

SEED = 20260713  # recorded; no stochastic operation in this build.

DRAFT = os.environ.get(
    "TCGA_DISCOVERY_ROOT",
    os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))
    ),
)
OUT = os.environ.get(
    "TCGA_DEFINITION_OUTPUT",
    os.path.join(DRAFT, "definitions"),
)

# ---- source paths (read-only external inputs) ----
SOURCE_DIR = os.environ.get(
    "TCGA_DEFINITION_SOURCE_DIR",
    os.path.join(DRAFT, "external-inputs", "definition-build"),
)
SRC_M1 = os.path.join(SOURCE_DIR, "module_M1_derived.txt")
SRC_DERIV_LEDGER = os.path.join(SOURCE_DIR, "derivation_ledger.tsv")
SRC_COMPACT = os.path.join(SOURCE_DIR, "module_M1_compact_sensitivity.txt")
SRC_M3P = os.path.join(SOURCE_DIR, "module_M3_primary_collectri_signed.txt")
SRC_M3P_LEDGER = os.path.join(SOURCE_DIR, "m3_primary_ledger.tsv")
SRC_M3S = os.path.join(SOURCE_DIR, "module_M3_sensitivity.txt")
SRC_COMP_REMOVED = os.path.join(
    SOURCE_DIR,
    "composition_removed_from_modules.tsv",
)
SRC_M4 = os.path.join(
    SOURCE_DIR,
    "module_M4_proliferation_progenitor_control.txt",
)

FIREWALL = ["CD74", "SAMHD1", "SPON1"]  # composition-overlap ESTIMATE genes


def load_members(path):
    """Ordered, de-duplicated (keep first occurrence) non-comment, non-blank lines."""
    seen = set()
    out = []
    with open(path) as fh:
        for line in fh:
            g = line.strip()
            if not g or g.startswith("#"):
                continue
            if g not in seen:
                seen.add(g)
                out.append(g)
    return out


def write_lines(path, header_lines, members):
    with open(path, "w") as fh:
        for h in header_lines:
            fh.write(h.rstrip("\n") + "\n")
        for m in members:
            fh.write(m + "\n")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ============================================================
# Load composition-removal map (which genes removed from which module)
# ============================================================
comp_rows = list(csv.DictReader(open(SRC_COMP_REMOVED), delimiter="\t"))
m1_removed = [r for r in comp_rows if r["module"] == "M1"]
m3_removed = [r for r in comp_rows if r["module"] == "M3_primary_collectri_signed"]
m3_removed_genes = set(r["gene"] for r in m3_removed)

produced = []  # (filename) in output order

# ============================================================
# 1 + 2. M1: source (239) and analysis-ready (236)
# ============================================================
m1 = load_members(SRC_M1)
assert len(m1) == 239, "M1 source expected 239, got %d" % len(m1)

M1_RELABEL = "E-MTAB-15475-derived fetal Mullerian epithelial module"

m1_src_header = [
    "# M1_source -- %s" % M1_RELABEL,
    "# RELABEL NOTE: this is an E-MTAB-15475-derived fetal Mullerian epithelial module,",
    "#   NOT a universal fetal-uterus program. Membership is dataset-derived and dataset-scoped.",
    "# STATUS: SEAL CANDIDATE (TASK-028 Freeze B). SHA-256 sealing is Vladimir-gated.",
    "# members = 239 (source, before composition firewall).",
]
p_m1_src = os.path.join(OUT, "M1_source.txt")
write_lines(p_m1_src, m1_src_header, m1)
produced.append("M1_source.txt")

firewall_in_m1 = [g for g in m1 if g in set(FIREWALL)]
m1_ar = [g for g in m1 if g not in set(FIREWALL)]
assert len(m1_ar) == 236, "M1 analysis-ready expected 236, got %d" % len(m1_ar)

m1_ar_header = [
    "# M1_analysis_ready -- %s" % M1_RELABEL,
    "# RELABEL NOTE: E-MTAB-15475-derived fetal Mullerian epithelial module (NOT a universal",
    "#   fetal-uterus program).",
    "# STATUS: SEAL CANDIDATE (TASK-028 Freeze B). SHA-256 sealing is Vladimir-gated.",
    "# COMPOSITION FIREWALL applied: removed {CD74, SAMHD1, SPON1} (composition-overlap ESTIMATE,",
    "#   source SI_geneset). See B1_exclusion_ledger.tsv for per-gene rows.",
    "# members = 236 = 239 minus 3 composition-overlap genes.",
]
p_m1_ar = os.path.join(OUT, "M1_analysis_ready.txt")
write_lines(p_m1_ar, m1_ar_header, m1_ar)
produced.append("M1_analysis_ready.txt")

# ============================================================
# 3. COMPACT-M1 (sensitivity-only): source (188) and analysis-ready (185)
# ============================================================
compact = load_members(SRC_COMPACT)
assert len(compact) == 188, "compact-M1 source expected 188, got %d" % len(compact)

compact_relabel_note = [
    "# RELABEL NOTE: E-MTAB-15475-derived fetal Mullerian epithelial module (compact, SENSITIVITY-only),",
    "#   NOT a universal fetal-uterus program.",
]

compact_src_header = [
    "# compact_M1_source -- E-MTAB-15475-derived fetal Mullerian epithelial module (COMPACT, SENSITIVITY-only)",
] + compact_relabel_note + [
    "# STATUS: SEAL CANDIDATE (TASK-028 Freeze B). SENSITIVITY-only; full 236/239 M1 stays PRIMARY.",
    "# RULE (pre-stated, stricter): donor_fraction_target >= 0.90 AND",
    "#   Mullerian-vs-other-fetal-epithelia log2FC >= 1.0. Size is an OUTPUT, not a target.",
    "# members = 188 (compact source, before composition firewall).",
]
p_compact_src = os.path.join(OUT, "compact_M1_source.txt")
write_lines(p_compact_src, compact_src_header, compact)
produced.append("compact_M1_source.txt")

firewall_in_compact = [g for g in compact if g in set(FIREWALL)]
compact_ar = [g for g in compact if g not in set(FIREWALL)]
# NO gene replacement -- size is simply 188 minus whatever firewall genes are present.

compact_ar_header = [
    "# compact_M1_analysis_ready -- E-MTAB-15475-derived fetal Mullerian epithelial module (COMPACT, SENSITIVITY-only)",
] + compact_relabel_note + [
    "# STATUS: SEAL CANDIDATE (TASK-028 Freeze B). SENSITIVITY-only.",
    "# COMPOSITION FIREWALL applied (SAME firewall as M1): removed the firewall genes",
    "#   {CD74, SAMHD1, SPON1} PRESENT in the 188. Present here: %s (n=%d). NO gene replacement."
    % (",".join(firewall_in_compact) if firewall_in_compact else "none", len(firewall_in_compact)),
    "# members = %d = 188 minus %d composition-overlap genes present in the compact list."
    % (len(compact_ar), len(firewall_in_compact)),
]
p_compact_ar = os.path.join(OUT, "compact_M1_analysis_ready.txt")
write_lines(p_compact_ar, compact_ar_header, compact_ar)
produced.append("compact_M1_analysis_ready.txt")

# ============================================================
# 4. M3 primary: source (603) and analysis-ready (577) + edge ledger with dual-split
# ============================================================
m3p = load_members(SRC_M3P)
assert len(m3p) == 603, "M3 primary source expected 603, got %d" % len(m3p)

m3p_src_header = [
    "# M3_primary_source -- CollecTRI SIGNED regulon (PRIMARY)",
    "# STATUS: SEAL CANDIDATE (TASK-028 Freeze B). SHA-256 sealing is Vladimir-gated.",
    "# members = 603 unique = 20 CORE_TF + 583 admitted CollecTRI-curated signed targets,",
    "#   in GENCODE v36, post-filter. (20 TFs, some also appear as targets.)",
]
p_m3p_src = os.path.join(OUT, "M3_primary_source.txt")
write_lines(p_m3p_src, m3p_src_header, m3p)
produced.append("M3_primary_source.txt")

m3p_ar = [g for g in m3p if g not in m3_removed_genes]
assert len(m3p_ar) == 577, "M3 primary analysis-ready expected 577, got %d" % len(m3p_ar)

m3p_ar_header = [
    "# M3_primary_analysis_ready -- CollecTRI SIGNED regulon (PRIMARY)",
    "# STATUS: SEAL CANDIDATE (TASK-028 Freeze B).",
    "# COMPOSITION FIREWALL applied: removed 26 composition-overlap genes (see B1_exclusion_ledger.tsv).",
    "# members = 577 = 603 minus 26 composition-overlap genes.",
]
p_m3p_ar = os.path.join(OUT, "M3_primary_analysis_ready.txt")
write_lines(p_m3p_ar, m3p_ar_header, m3p_ar)
produced.append("M3_primary_analysis_ready.txt")

# ---- edge ledger with dual-split ----
edge_rows_in = list(csv.DictReader(open(SRC_M3P_LEDGER), delimiter="\t"))
admitted = [r for r in edge_rows_in if r["decision"] == "admit"]
assert len(admitted) == 802, "expected 802 admitted edges, got %d" % len(admitted)


def is_true(x):
    return str(x).strip().lower() == "true"


def confidence_for(sign, r):
    """Per-edge likelihood weight from CollecTRI consensus flags.
    weight=1.0 if the assigned sign is a CollecTRI CONSENSUS call
    (consensus_stimulation for +1, consensus_inhibition for -1);
    weight=0.5 if it is a non-consensus curated call. Deterministic."""
    if sign == 1:
        return 1.0 if is_true(r["collectri_consensus_stimulation"]) else 0.5
    if sign == -1:
        return 1.0 if is_true(r["collectri_consensus_inhibition"]) else 0.5
    raise ValueError("confidence_for called with sign=%r" % sign)


# Build the output edge ledger. Sign +1 -> activation, -1 -> repression.
# Dual edges (source sign == 0) are the 41 dual activation/repression edges;
# each is emitted as TWO signed rows: one +1 and one -1.
edge_header_comment = [
    "# m3_primary_edge_ledger -- M3 PRIMARY CollecTRI-signed edge ledger (SEAL CANDIDATE, TASK-028 Freeze B)",
    "# ONE ROW PER SIGNED EDGE CONTRIBUTION.",
    "# EDGE-WEIGHT / CONFIDENCE RULE:",
    "#   VIPER/aREA mode-of-regulation (MoR) = CollecTRI sign (+1 activation, -1 repression).",
    "#   Per-edge likelihood weight = CollecTRI CONSENSUS confidence, derived deterministically",
    "#   from the pinned CollecTRI columns: weight=1.0 when the assigned sign is a CollecTRI",
    "#   CONSENSUS call (column collectri_consensus_stimulation=True for +1, or",
    "#   collectri_consensus_inhibition=True for -1); weight=0.5 for a non-consensus curated call.",
    "#   evidence_class = CollecTRI-curated (all admitted primary edges).",
    "# DUAL EDGES: the 41 dual activation/repression edges (source sign==0, sign_label in",
    "#   {unknown, dual/unknown}) enter as TWO SEPARATE SIGNED CONTRIBUTIONS -- one +1 row and one",
    "#   -1 row -- NOT collapsed to a single sign. Each half carries its own consensus-derived weight",
    "#   (consensus_stimulation for the +1 half, consensus_inhibition for the -1 half; 0.5 if neither).",
    "#   dual_split=Y flags both halves of a dual edge.",
    "# Columns: TF\ttarget\tsign\tsign_label\tweight_confidence\tevidence_class\tdual_split",
]

out_edges = []  # tuples for writing
n_single = 0
n_dual = 0
for r in admitted:
    tf = r["TF"]
    target = r["target_gencode_symbol"] or r["target"]
    src_sign = int(r["sign"])
    ev = r["evidence_class"]
    if src_sign in (1, -1):
        w = confidence_for(src_sign, r)
        lbl = "activation" if src_sign == 1 else "repression"
        out_edges.append((tf, target, src_sign, lbl, ("%.1f" % w), ev, "N"))
        n_single += 1
    elif src_sign == 0:
        # dual: emit +1 then -1
        w_pos = confidence_for(1, r)
        w_neg = confidence_for(-1, r)
        out_edges.append((tf, target, 1, "activation", ("%.1f" % w_pos), ev, "Y"))
        out_edges.append((tf, target, -1, "repression", ("%.1f" % w_neg), ev, "Y"))
        n_dual += 1
    else:
        raise ValueError("unexpected sign %r" % src_sign)

assert n_dual == 41, "expected 41 dual edges, got %d" % n_dual
total_edges = len(out_edges)
# 761 single (598 act + 163 rep) + 41 dual*2 = 761 + 82 = 843
assert total_edges == n_single + 2 * n_dual, "edge count mismatch"

p_m3p_edge = os.path.join(OUT, "m3_primary_edge_ledger.tsv")
with open(p_m3p_edge, "w") as fh:
    for h in edge_header_comment:
        fh.write(h + "\n")
    fh.write("TF\ttarget\tsign\tsign_label\tweight_confidence\tevidence_class\tdual_split\n")
    for tf, target, sign, lbl, w, ev, dsp in out_edges:
        fh.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\n" % (tf, target, sign, lbl, w, ev, dsp))
produced.append("m3_primary_edge_ledger.tsv")

# ============================================================
# 5. M3 sensitivity (809) -- pre-frozen sensitivity label
# ============================================================
m3s = load_members(SRC_M3S)
assert len(m3s) == 809, "M3 sensitivity expected 809, got %d" % len(m3s)

m3s_header = [
    "# M3_sensitivity -- M3 CollecTRI/DoRothEA regulon SENSITIVITY variant",
    "# LABEL: PRE-FROZEN SENSITIVITY (DoRothEA A/B + CollecTRI ChIP/motif), NOT primary.",
    "# STATUS: SEAL CANDIDATE (TASK-028 Freeze B). Separate from primary; NOT merged into primary.",
    "# members = 809 (CORE_TF + DoRothEA(A/B) targets + CollecTRI ChIP/motif-evidenced targets, GENCODE v36).",
]
p_m3s = os.path.join(OUT, "M3_sensitivity.txt")
write_lines(p_m3s, m3s_header, m3s)
produced.append("M3_sensitivity.txt")

# ============================================================
# 6. M4 covariate (30) -- covariate-only label
# ============================================================
m4 = load_members(SRC_M4)
assert len(m4) == 30, "M4 expected 30, got %d" % len(m4)

m4_header = [
    "# M4_covariate -- proliferation / progenitor control",
    "# LABEL: COVARIATE-ONLY, not a developmental claim module.",
    "# STATUS: SEAL CANDIDATE (TASK-028 Freeze B). Used to partial out proliferation (RISK R-10).",
    "# members = 30 (compact literature-standard proliferation/cell-cycle seed).",
]
p_m4 = os.path.join(OUT, "M4_covariate.txt")
write_lines(p_m4, m4_header, m4)
produced.append("M4_covariate.txt")

# ============================================================
# 7. Consolidated final ledgers
# ============================================================
# 7a. B1_gene_level_ledger.tsv -- one row per (module, gene) with inclusion/exclusion status.
#     Modules covered: M1, compact_M1, M3_primary, M3_sensitivity, M4.
gene_ledger_path = os.path.join(OUT, "B1_gene_level_ledger.tsv")
gl_cols = ["module", "gene", "membership", "status", "reason", "composition_signature", "source"]

# map of firewall signature (from composition_removed file) for reasons
m1_sig = {r["gene"]: r["estimate_signature"] for r in m1_removed}
m3_sig = {r["gene"]: r["estimate_signature"] for r in m3_removed}

gl_rows = []

# M1 (source scope): every source gene, marked included or excluded-by-firewall
for g in m1:
    if g in set(FIREWALL):
        gl_rows.append(["M1", g, "source", "EXCLUDED_analysis_ready",
                        "composition-overlap ESTIMATE (%s)" % m1_sig.get(g, "immune/stromal"),
                        m1_sig.get(g, ""), "SI_geneset"])
    else:
        gl_rows.append(["M1", g, "source", "INCLUDED_analysis_ready",
                        "fetal-Mullerian-epithelial derived member (E-MTAB-15475)", "", "E-MTAB-15475-derived"])

# compact_M1 (source scope): 188, mark firewall-present as excluded
for g in compact:
    if g in set(FIREWALL):
        gl_rows.append(["compact_M1", g, "source", "EXCLUDED_analysis_ready",
                        "composition-overlap ESTIMATE (%s)" % m1_sig.get(g, "immune/stromal"),
                        m1_sig.get(g, ""), "SI_geneset"])
    else:
        gl_rows.append(["compact_M1", g, "source", "INCLUDED_analysis_ready",
                        "compact sensitivity member (donor_fraction>=0.90 & log2FC>=1.0)", "", "E-MTAB-15475-derived"])

# M3_primary (member scope): 603, mark 26 firewall as excluded
for g in m3p:
    if g in m3_removed_genes:
        gl_rows.append(["M3_primary", g, "source", "EXCLUDED_analysis_ready",
                        "composition-overlap ESTIMATE (%s)" % m3_sig.get(g, "immune/stromal"),
                        m3_sig.get(g, ""), "SI_geneset"])
    else:
        gl_rows.append(["M3_primary", g, "source", "INCLUDED_analysis_ready",
                        "CollecTRI-signed regulon member (CORE_TF or admitted target)", "", "CollecTRI + GENCODE v36"])

# M3_sensitivity (member scope): 809, all included (sensitivity variant, no firewall applied per spec)
for g in m3s:
    gl_rows.append(["M3_sensitivity", g, "member", "INCLUDED_sensitivity",
                    "PRE-FROZEN SENSITIVITY (DoRothEA A/B + CollecTRI ChIP/motif), not primary", "",
                    "DoRothEA A/B + CollecTRI ChIP/motif"])

# M4 (member scope): 30, covariate-only
for g in m4:
    gl_rows.append(["M4", g, "member", "INCLUDED_covariate",
                    "COVARIATE-ONLY proliferation/progenitor control, not a developmental claim", "",
                    "literature-standard cell-cycle seed"])

with open(gene_ledger_path, "w") as fh:
    fh.write("\t".join(gl_cols) + "\n")
    for row in gl_rows:
        fh.write("\t".join(str(x) for x in row) + "\n")
produced.append("B1_gene_level_ledger.tsv")

# 7b. B1_exclusion_ledger.tsv -- every removed gene (module, gene, reason, signature, source)
excl_path = os.path.join(OUT, "B1_exclusion_ledger.tsv")
excl_cols = ["module", "gene", "reason", "composition_signature", "source"]
excl_rows = []
# M1 firewall removals (3)
for g in firewall_in_m1:
    excl_rows.append(["M1", g, "composition-overlap ESTIMATE", m1_sig.get(g, ""), "SI_geneset"])
# compact_M1 firewall removals (whatever present)
for g in firewall_in_compact:
    excl_rows.append(["compact_M1", g, "composition-overlap ESTIMATE", m1_sig.get(g, ""), "SI_geneset"])
# M3_primary firewall removals (26) -- preserve source-file order
for r in m3_removed:
    excl_rows.append(["M3_primary", r["gene"], "composition-overlap ESTIMATE",
                      r["estimate_signature"], "SI_geneset"])
with open(excl_path, "w") as fh:
    fh.write("\t".join(excl_cols) + "\n")
    for row in excl_rows:
        fh.write("\t".join(str(x) for x in row) + "\n")
produced.append("B1_exclusion_ledger.tsv")

# also emit this build script's own presence in manifest (self-documenting)
produced_pre_manifest = list(produced)
produced_pre_manifest.append("build_sealed_B1.py")

# ============================================================
# 8. SHA-256 manifest of every produced file (+ this script)
# ============================================================
manifest_path = os.path.join(OUT, "B1_SHA256.manifest")
sha_map = {}
for fn in produced_pre_manifest:
    fp = os.path.join(OUT, fn)
    sha_map[fn] = sha256_file(fp)
with open(manifest_path, "w") as fh:
    for fn in produced_pre_manifest:
        fh.write("%s  %s\n" % (fn, sha_map[fn]))

# print a machine-readable summary to stdout
print("SEED=%d" % SEED)
print("compact_firewall_present=%s n=%d" % (",".join(firewall_in_compact), len(firewall_in_compact)))
print("compact_analysis_ready_size=%d" % len(compact_ar))
print("m1_analysis_ready_size=%d" % len(m1_ar))
print("m3_primary_source_size=%d" % len(m3p))
print("m3_primary_analysis_ready_size=%d" % len(m3p_ar))
print("m3_edges_single=%d m3_edges_dual=%d m3_total_edges_after_split=%d" % (n_single, n_dual, total_edges))
print("gene_level_ledger_rows=%d" % len(gl_rows))
print("exclusion_ledger_rows=%d" % len(excl_rows))
print("---SHA256---")
for fn in produced_pre_manifest:
    print("%s  %s" % (fn, sha_map[fn]))
# manifest's own sha (computed after write, reported separately)
print("B1_SHA256.manifest  %s" % sha256_file(manifest_path))
