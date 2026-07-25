"""
TASK-028 Freeze B v3 PRIMARY=C -- common utilities (LOCKED MECHANICAL RUN, resume).
Uses the public definitions directory. Implements the v3 purity resolution (CPE complete-case n=506;
no-purity n=507; ABSOLUTE complete-case n=502). No tuning. Deterministic. ASCII only.

Master seed: 20260713 (B2.11). Sub-seeds derived DETERMINISTICALLY from (master_seed,
step_id) via sha256; no wall-clock seeding.
git HEAD (frozen, unchanged): 83503bad47b60193598b2b9ebe819c22c83e8ac1
"""
import os
import re
import hashlib
import numpy as np
import pandas as pd

MASTER_SEED = 20260713
GIT_HEAD = "83503bad47b60193598b2b9ebe819c22c83e8ac1"

# ---- public package paths and external inputs ----
ROOT = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
SEALED = ROOT
B1 = os.path.join(ROOT, "definitions")
ORIGINALS = os.environ.get(
    "TCGA_STAR_ORIGINALS",
    os.path.join(ROOT, "external-inputs", "originals"),
)
COHORT_FILE = os.environ.get(
    "TCGA_COHORT_FILE",
    os.path.join(
        ROOT,
        "external-inputs",
        "cohort_selected_primary.tsv",
    ),
)
SUBTYPE_TABLE = os.environ.get(
    "TCGA_SUBTYPE_TABLE",
    os.path.join(
        ROOT,
        "external-inputs",
        "subtype_normalized.tsv",
    ),
)
SRC = os.environ.get(
    "TCGA_DISCOVERY_SOURCE_DIR",
    os.path.join(ROOT, "external-inputs"),
)
ARAN_XLSX = os.path.join(SRC, "aran_cpe.xlsx")
ABSOLUTE_TXT = os.path.join(SRC, "absolute_purity.txt")
ESTIMATE_GMT = os.path.join(SRC, "ESTIMATE_SI_geneset.gmt")

WORK = os.environ.get(
    "TCGA_DISCOVERY_WORK_DIR",
    os.path.join(ROOT, "work"),
)
RESULTS_V3 = os.environ.get(
    "TCGA_DISCOVERY_RESULTS_DIR",
    os.path.join(ROOT, "results"),
)
INTER = os.environ.get(
    "TCGA_DISCOVERY_INTERMEDIATE_DIR",
    os.path.join(WORK, "intermediate"),
)
LOGS = os.environ.get(
    "TCGA_DISCOVERY_LOG_DIR",
    os.path.join(WORK, "logs"),
)

SUBTYPES = ["POLE", "MMRd", "NSMP", "p53abn"]

# The single CPE-non-finite patient dropped from CPE-including models ONLY (v3 amendment).
CPE_NONFINITE_PATIENT = "TCGA-BS-A0TG"
CPE_DROP_REASON = "primary purity covariate non-finite"

PATIENT_RE = re.compile(r"^(TCGA-[0-9A-Za-z]{2}-[0-9A-Za-z]{4})")


def patient_of(barcode):
    m = PATIENT_RE.match(str(barcode).strip())
    return m.group(1) if m else None


def substep_seed(step_id):
    """Deterministic sub-seed from (master_seed, step_id) -- B2.11."""
    h = hashlib.sha256(f"{MASTER_SEED}:{step_id}".encode()).hexdigest()
    return int(h[:8], 16)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_module_genes(path):
    genes = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            genes.append(s)
    return genes


def load_cohort():
    df = pd.read_csv(COHORT_FILE, sep="\t", dtype=str)
    assert list(df.columns) == ["patient_barcode", "kept_uuid",
                                "kept_aliquot_barcode", "mapped_4way"], df.columns
    return df


def load_star_tpm(uuid):
    path = os.path.join(ORIGINALS, f"{uuid}__augmented_star_gene_counts.tsv")
    df = pd.read_csv(path, sep="\t", comment=None, skiprows=1, dtype={
        "gene_id": str, "gene_name": str, "gene_type": str})
    df = df[df["gene_id"].str.startswith("ENSG")].copy()
    return df


def ensembl_unversioned(gid):
    return gid.split(".")[0]


def is_finite_value(v):
    """True iff v is a real finite numeric value (not None, not NaN/NA/NULL/blank/'.')."""
    if v is None:
        return False
    if isinstance(v, (int, float)):
        return v == v
    s = str(v).strip()
    if s == "" or s.upper() in ("NA", "NAN", "NULL", "."):
        return False
    try:
        f = float(s)
        return f == f
    except ValueError:
        return False


def load_cpe_by_patient():
    """Aran CPE (Supp Data 1, header row 4, column 'CPE'), patient-level, finite only.
    Union across aliquots (patient finite if ANY aliquot finite). Returns dict pat->value.
    Matches the sealed ABSOLUTE_PROVENANCE.md CPE re-confirmation (n=506)."""
    import openpyxl
    wb = openpyxl.load_workbook(ARAN_XLSX, read_only=True, data_only=True)
    ws = wb["Supp Data 1"]
    rows = ws.iter_rows(min_row=4, values_only=True)
    header = next(rows)
    header = [str(h).strip() if h is not None else None for h in header]
    idx_sample = header.index("Sample ID")
    idx_cpe = header.index("CPE")
    cpe = {}
    for row in rows:
        sample = row[idx_sample]
        pat = patient_of(sample)
        if pat is None:
            continue
        v = row[idx_cpe]
        if is_finite_value(v):
            # union: keep first finite (aliquots of a patient carry the same CPE)
            if pat not in cpe:
                cpe[pat] = float(str(v).strip())
    wb.close()
    return cpe


def load_absolute_by_patient():
    """PanCanAtlas ABSOLUTE (TCGA_mastercalls...), column 'purity', patient-level finite.
    patient key from 'sample' (aliquot) with 'array' fallback. Matches the sealed
    ABSOLUTE_PROVENANCE.md PanCanAtlas re-derivation (n=502)."""
    absv = {}
    with open(ABSOLUTE_TXT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx_array = header.index("array")
        idx_sample = header.index("sample")
        idx_purity = header.index("purity")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_array, idx_sample, idx_purity):
                continue
            pat = patient_of(parts[idx_sample]) or patient_of(parts[idx_array])
            if pat is None:
                continue
            if is_finite_value(parts[idx_purity]):
                if pat not in absv:
                    absv[pat] = float(parts[idx_purity].strip())
    return absv


def load_estimate_signatures():
    """Return dict {StromalSignature:[genes], ImmuneSignature:[genes]} from the GMT."""
    sig = {}
    with open(ESTIMATE_GMT) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            sig[p[0]] = p[2:]
    return sig
