"""
TASK-028 Freeze B v2 -- scoring methods (LOCKED; pinned params).

ssGSEA: Barbie et al. 2009 rank-based single-sample GSEA.
  - per sample, rank genes by expression (ascending), rank^alpha weighting alpha=0.25
  - ES = sum over gene-set positions of the running-sum deviation (integral form),
    normalized by (n_genes) range as in the standard ssGSEA (GSVA) implementation.
  - normalized across samples by range (max-min) of the ES vector (ssgsea.norm=TRUE
    style). We report BOTH raw and range-normalized; primary uses range-normalized
    (the standard ssGSEA output). Deterministic; no randomness.

VIPER/aREA: Alvarez et al. 2016 analytic rank-based enrichment (aREA, 1-tail + 2-tail
  signature). Signed mode-of-regulation (MoR = consensus sign) with per-edge likelihood
  weights (consensus confidence). Regulon NES computed per sample from the gene
  expression signature (z-scored across samples). Deterministic.

All params are FIXED a priori (B2.3). No tuning.
"""
import numpy as np
from scipy.stats import rankdata, norm

SSGSEA_ALPHA = 0.25  # pinned


def ssgsea_matrix(expr, genesets, alpha=SSGSEA_ALPHA, norm_range=True):
    """
    expr : ndarray genes x samples (log2 TPM+1)
    genesets : dict name -> list of ROW indices into expr
    returns : dict name -> ndarray (samples,) ssGSEA enrichment score
    Standard ssGSEA (Barbie/GSVA ssgsea): for each sample, rank genes ascending;
    weight in-set genes by |rank_expr|^alpha; running sum of (weighted in-set minus
    uniform out-set); ES = sum of running-sum deviations (integral), then range-norm.
    """
    n_genes, n_samples = expr.shape
    # ranks per sample: ascending, ties averaged then argsort ordering
    # We follow the GSVA ssgsea: order genes by expression; use rank as the statistic.
    out = {name: np.zeros(n_samples) for name in genesets}
    # Precompute per-sample ordering
    for s in range(n_samples):
        x = expr[:, s]
        order = np.argsort(x, kind="mergesort")  # ascending, stable/deterministic
        # rank_stat: position rank 1..n in ascending order -> larger = higher expr
        ranks = np.empty(n_genes, dtype=np.float64)
        ranks[order] = np.arange(1, n_genes + 1, dtype=np.float64)
        # weight = rank^alpha (Barbie ssGSEA uses |expr rank|^alpha)
        for name, gset in genesets.items():
            if len(gset) == 0:
                out[name][s] = np.nan
                continue
            in_set = np.zeros(n_genes, dtype=bool)
            in_set[gset] = True
            # process genes from highest expr to lowest (descending order of ranks)
            desc = order[::-1]
            w = np.power(ranks[desc], alpha)
            hits = in_set[desc]
            # P_hit running: cumulative weighted hits / total weighted hits
            wsum_hit = w[hits].sum()
            if wsum_hit == 0:
                out[name][s] = np.nan
                continue
            step_up = np.where(hits, w / wsum_hit, 0.0)
            n_miss = n_genes - hits.sum()
            step_dn = np.where(hits, 0.0, 1.0 / n_miss)
            running = np.cumsum(step_up - step_dn)
            # ssGSEA ES = sum of running (integral form, GSVA ssgsea)
            out[name][s] = running.sum()
    if norm_range:
        for name in out:
            v = out[name]
            finite = np.isfinite(v)
            if finite.sum() > 0:
                rng = np.nanmax(v) - np.nanmin(v)
                if rng > 0:
                    out[name] = v / rng
    return out


def _make_signature_z(expr):
    """Gene expression signature = z-score across samples (per gene). genes x samples."""
    mu = expr.mean(axis=1, keepdims=True)
    sd = expr.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = np.nan
    z = (expr - mu) / sd
    return z  # genes x samples


def area_regulon(expr, gene_index, regulon, weights, mor):
    """
    aREA (Alvarez et al. 2016, viper::aREA) SIGNED single-sample regulon activity.
    expr : genes x samples (log2 TPM+1)
    gene_index : dict gene_symbol -> single best row index into expr
    regulon : list of target gene symbols
    weights : list of per-edge likelihood weights (consensus confidence 1.0/0.5)
    mor : list of per-edge mode-of-regulation sign (+1 activation / -1 repression)
    Returns: (nes ndarray (samples,), n_targets_used)

    Canonical aREA:
      per sample, signature = z-scored expression across samples.
      t2 = qnorm( rank(sig)/(N+1) )                 (directional / two-tail position)
      t1 = qnorm( (rank(|sig|)/(N+1) + 1)/2 )        (magnitude / one-tail)
      per-edge footprint weight  w_i (likelihood), MoR sign m_i.
      Normalize footprint: wts = w / max(w); wtss = wts / sum(wts).
      Two-tail enrichment (directional):
          s2 = sum_i ( wtss_i * m_i * t2_i )
      One-tail enrichment (undirected magnitude):
          s1 = sum_i ( wtss_i * t1_i )
      Combined enrichment (viper aREA integration):
          es = (|s2| + s1*(s1>0)) * sign(s2)
      NES = es * sqrt( sum_i wts_i^2 )   (likelihood-based normalization to ~N(0,1))
    Deterministic (rank-based; no randomness). Missing/non-finite targets dropped.
    """
    z = _make_signature_z(expr)  # genes x samples (zero-variance genes -> NaN rows)
    n_samples = expr.shape[1]
    # The gene expression SIGNATURE universe = genes with a finite z (measured / non-
    # zero-variance). Zero-variance (all-zero / undetectable) genes carry no signal and
    # are EXCLUDED from the rank universe -- they must not poison the ranking (canonical
    # viper::aREA ranks the measured signature). This is the same "undetectable" exclusion
    # already applied to module genes (B2.2 all-zero rule), applied to the ranking universe.
    finite_gene = np.all(np.isfinite(z), axis=1)   # genes present in every sample's signature
    # position of each ORIGINAL row index within the reduced (finite) universe
    reduced_pos = -np.ones(z.shape[0], dtype=np.int64)
    reduced_pos[finite_gene] = np.arange(int(finite_gene.sum()))
    z_red = z[finite_gene]                          # (Nfin x samples)
    N = z_red.shape[0]
    rows, w, m = [], [], []
    for t, wt, mo in zip(regulon, weights, mor):
        r = gene_index.get(t)
        if r is None:
            continue
        if not finite_gene[r]:
            continue          # target itself undetectable/zero-variance -> contributes nothing
        rows.append(reduced_pos[r]); w.append(float(wt)); m.append(float(mo))
    if len(rows) == 0:
        return np.full(n_samples, np.nan), 0
    rows = np.array(rows); w = np.array(w); m = np.array(m)
    # footprint weight normalization (viper aREA)
    wmax = w.max()
    if wmax <= 0:
        return np.full(n_samples, np.nan), 0
    wts = w / wmax
    swts = wts.sum()
    wtss = wts / swts
    nes_norm = np.sqrt(np.sum(wts ** 2))
    nes = np.zeros(n_samples)
    for s in range(n_samples):
        sig = z_red[:, s]
        r_signed = rankdata(sig, method="average") / (N + 1.0)
        t2 = norm.ppf(r_signed)[rows]
        r_abs = rankdata(np.abs(sig), method="average") / (N + 1.0)
        t1 = norm.ppf((r_abs + 1.0) / 2.0)[rows]
        s2 = np.sum(wtss * m * t2)
        s1 = np.sum(wtss * t1)
        es = (abs(s2) + (s1 if s1 > 0 else 0.0)) * np.sign(s2)
        nes[s] = es * nes_norm
    return nes, len(rows)


def signed_ssgsea_regulon(expr, gene_index, regulon, mor):
    """
    SIGNED-ssGSEA = the pre-frozen SENSITIVITY scoring method for M3 regulons (B2.3).
    Standard ssGSEA computed separately over the ACTIVATION target set (mor=+1) and the
    REPRESSION target set (mor=-1) using the SAME frozen ssGSEA params (alpha=0.25,
    range-normalized), then combined as a SIGNED difference:
        signed_score = ES(activation targets) - ES(repression targets).
    A regulon with only one sign-set uses that set alone (the other ES is 0). This is the
    ssGSEA analogue of the signed VIPER/aREA primary; deterministic, no randomness, no tuning.
    Returns (score ndarray (samples,), n_pos_used, n_neg_used).
    """
    pos_rows, neg_rows = [], []
    used = set()
    for t, mo in zip(regulon, mor):
        r = gene_index.get(t)
        if r is None or r in used:
            continue
        used.add(r)
        if mo > 0:
            pos_rows.append(r)
        elif mo < 0:
            neg_rows.append(r)
    sets = {}
    if pos_rows:
        sets["POS"] = pos_rows
    if neg_rows:
        sets["NEG"] = neg_rows
    if not sets:
        return np.full(expr.shape[1], np.nan), 0, 0
    es = ssgsea_matrix(expr, sets, alpha=SSGSEA_ALPHA, norm_range=True)
    pos = es["POS"] if "POS" in es else np.zeros(expr.shape[1])
    neg = es["NEG"] if "NEG" in es else np.zeros(expr.shape[1])
    signed = pos - neg
    return signed, len(pos_rows), len(neg_rows)
