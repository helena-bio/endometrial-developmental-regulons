# Molecular-subtype sensitivity

This package contains the post-hoc decomposition of the frozen C2
molecular-subtype contrast for GATA2 and SOX9.

The analysis asks whether the pooled negative C2 signal is contributed by
both POLE and MMRd relative to NSMP, or whether it is driven by one class,
sample-size composition, or opposing class-specific effects.

## Frozen contrast

C2 is defined as:

    0.5 * (POLE - NSMP) + 0.5 * (MMRd - NSMP)

The weights are fixed biological contrast weights. They are not
sample-size weights.

Observed POLE and MMRd class sizes therefore cannot mechanically reweight
the pooled contrast.

## Analysis design

The package estimates three contrasts from each unchanged full model:

- POLE versus NSMP;
- MMRd versus NSMP;
- POLE versus MMRd.

Pairwise subsets are not refitted.

The decomposition is performed for:

- TCGA primary;
- TCGA no-purity;
- CPTAC Discovery;
- CPTAC Confirmatory;
- CPTAC fixed-effect meta-analysis.

Separate descriptive BH-18 families are retained for documentation. They
do not replace or modify the frozen multiplicity families.

## Main result

Both class-versus-NSMP estimates are negative for GATA2 and SOX9 across
the TCGA models, both CPTAC strata, and the CPTAC meta-analysis.

For GATA2:

- TCGA primary POLE versus NSMP: d = -0.682;
- TCGA primary MMRd versus NSMP: d = -0.424;
- CPTAC meta POLE versus NSMP: d = -1.071;
- CPTAC meta MMRd versus NSMP: d = -0.620.

For SOX9:

- TCGA primary POLE versus NSMP: d = -0.690;
- TCGA primary MMRd versus NSMP: d = -0.364;
- CPTAC meta POLE versus NSMP: d = -1.050;
- CPTAC meta MMRd versus NSMP: d = -0.607.

The direct CPTAC meta POLE-versus-MMRd intervals include zero for both
targets.

The larger POLE point estimates therefore do not establish a definitive
POLE-versus-MMRd difference.

The CPTAC POLE strata are small, which limits precision.

## Package contents

The package includes:

- the hypothesis and analysis plan;
- a machine-readable frozen specification;
- executable analysis and release-building code;
- per-class contrast tables in JSON, CSV, and TSV formats;
- C2 reconstruction checks;
- CPTAC fixed-effect meta-analysis tables;
- heterogeneity and interpretation-taxonomy tables;
- GATA2 and SOX9 forest plots in PNG, PDF, and SVG formats;
- two deterministic reproducibility runs;
- independent validation and presentation-audit records.

Key directories include:

    protocol/
    src/
    results/
    figures/
    reproducibility/
    validation/

## Interpretation boundary

This analysis is post-hoc and descriptive.

It may explain the composition of C2, but it cannot:

- change a frozen TCGA category;
- grant deletion-robustness credit;
- rescue a failed external target;
- replace a frozen CPTAC verdict;
- establish equality when an interval spans zero;
- establish causal subtype-specific regulation.

## Reproduction

Read the package protocol and repository-level
[REPRODUCIBILITY.md](../../REPRODUCIBILITY.md) before execution.

A rerun must preserve the documented full-model design, contrast
construction, residual standardisation, seeds, source mappings, and
meta-analysis procedure.
