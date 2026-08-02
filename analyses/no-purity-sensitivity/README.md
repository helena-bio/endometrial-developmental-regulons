# No-purity sensitivity

This package contains the prespecified TCGA-UCEC sensitivity analysis
that repeats the six-target evaluation without tumour purity in the
model.

Its purpose is to test whether omission of the primary purity covariate
materially changes the direction, magnitude, pointwise category, or
single-target-deletion behaviour of the reported signals.

## Scientific scope

The package evaluates six frozen targets:

- GATA2;
- SOX9;
- HOXA9;
- WT1;
- PAX8;
- LHX1.

The same biological contrast definitions used in the primary TCGA
analysis are retained.

The no-purity model is a sensitivity analysis. It does not replace the
primary purity-adjusted model.

## Main result

All six targets retain their primary direction and pointwise category in
the no-purity model.

For the two externally concordant C2 targets:

- GATA2 changes from d = -0.553 to d = -0.539;
- SOX9 changes from d = -0.527 to d = -0.509.

The small changes do not explain the larger CPTAC effect estimates.

This stability does not establish purity independence. It shows only that
removing the documented TCGA purity covariate produces limited change
within this model.

## Universal single-target-deletion analysis

The complete mapped-target leave-one-out analysis is repeated under the
no-purity specification.

Universal robustness credit requires every deletion to retain:

- the locked direction;
- absolute standardized effect at least 0.50;
- a confidence interval excluding zero;
- the omnibus gate;
- the full gated multiplicity criterion;
- the locked scoring and adjustment condition.

No-purity deletion-floor failures are:

- GATA2: 2 of 75;
- SOX9: 28 of 70;
- HOXA9: 2 of 19;
- WT1: 1 of 92;
- PAX8: 0 of 23;
- LHX1: 0 of 6.

PAX8 and LHX1 retain universal deletion credit. GATA2, SOX9, HOXA9,
and WT1 do not receive universal credit.

## Package contents

The package includes:

- protocol revisions and seed records;
- executable no-purity and leave-one-out code;
- machine-readable target results;
- effect and purity decomposition tables;
- analytical reports;
- deterministic rerun records;
- independent validation and targeted verdict records.

Key directories include:

    protocol/
    results/
    validation/

## Interpretation boundary

The package does not establish:

- causal transcription-factor activity;
- direct DNA binding;
- purity independence;
- biomarker validity;
- clinical utility;
- treatment prediction;
- therapeutic dependence.

The no-purity result is one model-sensitivity check within the frozen
evidence framework.

## Reproduction

Read the package protocol and repository-level
[REPRODUCIBILITY.md](../../REPRODUCIBILITY.md) before execution.

A rerun must use the documented TCGA cohort, source mappings, target
definitions, contrast weights, seeds, thresholds, and dependency
environment.
