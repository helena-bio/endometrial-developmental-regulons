# TCGA primary discovery

This package contains the frozen TCGA-UCEC discovery analysis for the
study:

> Robustness and cross-cohort concordance of developmental regulons in
> endometrial carcinoma

The analysis evaluates one fetal Mullerian epithelial module and a
prespecified panel of 20 developmental transcription-factor regulons in
507 primary TCGA-UCEC tumours.

## Scientific objective

The discovery analysis asks whether molecular subtype is associated with:

- a broad fetal Mullerian epithelial programme;
- selected developmental transcription-factor target-expression patterns;
- effects that remain stable when individual mapped targets are removed.

Pointwise subtype association and universal single-target-deletion
robustness are treated as separate evidential properties.

## Cohort

The frozen cohort contains:

- 49 POLE-ultramutated tumours;
- 148 mismatch-repair-deficient tumours, abbreviated MMRd;
- 147 no specific molecular profile tumours, abbreviated NSMP;
- 163 p53-abnormal tumours.

The primary purity-adjusted complete-case model contains 506 tumours.

Prespecified sensitivity configurations include:

- a no-purity model with 507 tumours;
- a PanCanAtlas ABSOLUTE-purity sensitivity with 502 tumours.

One primary-tumour analytical unit is retained per patient.

## Molecular-subtype contrasts

Three contrasts are defined:

- C1 compares p53-abnormal tumours with equal biological weights across
  POLE, MMRd, and NSMP;
- C2 compares `0.5 * POLE + 0.5 * MMRd` with NSMP;
- C3 compares POLE with MMRd and is descriptive.

Observed group sizes do not alter the contrast weights.

A three-degree-of-freedom molecular-subtype omnibus test is applied before
a gated contrast claim is considered.

## Modules and regulons

The package includes:

- the broad fetal Mullerian epithelial module;
- the frozen 20-factor developmental transcription-factor inventory;
- signed CollecTRI consensus regulons;
- module and factor provenance;
- mapped-target and edge records;
- frozen scoring and exclusion rules.

The broad fetal module and the individual regulons represent different
biological hypotheses and are interpreted separately.

## Primary model

The primary model includes:

- molecular subtype;
- the frozen M4 proliferation covariate;
- ESTIMATE stromal score;
- ESTIMATE immune score;
- Aran consensus tumour purity.

The standardized effect is the contrast coefficient divided by the model
residual standard deviation.

Bootstrap confidence intervals and subtype-label permutations use the
frozen seed scheme recorded in the package.

## Pointwise category

A positive pointwise category requires:

- the locked direction;
- absolute standardized effect at least 0.50;
- a confidence interval excluding zero;
- the omnibus gate;
- the full gated multiplicity criterion.

C1 equivalence uses the frozen plus-or-minus 0.30 region and the
confidence-interval-in-region rule.

## Universal single-target-deletion robustness

The complete analysis is repeated after removing each unique mapped target
gene one at a time.

Universal robustness credit requires every deletion to retain:

- the locked direction;
- absolute standardized effect at least 0.50;
- a confidence interval excluding zero;
- the omnibus gate;
- the full gated multiplicity criterion;
- the locked scoring and adjustment condition.

One failed deletion is sufficient to withhold universal credit.

## Main result

For C1:

- LHX1 shows d = 0.962 and retains universal deletion credit;
- PAX8 shows d = 0.693 and retains universal deletion credit.

For C2:

- GATA2 shows d = -0.553;
- SOX9 shows d = -0.527;
- HOXA9 shows d = -0.646;
- WT1 shows d = -0.591.

All four C2 targets meet the pointwise criteria, but none receives
universal single-target-deletion credit.

Primary deletion-floor failures are:

- GATA2: 1 of 75;
- SOX9: 15 of 70;
- HOXA9: 2 of 19;
- WT1: 1 of 92.

The broad fetal Mullerian epithelial module is non-confirmatory under the
prespecified contrast and equivalence rules.

## Package contents

The package includes:

- frozen model and claim-boundary specifications;
- amendments and provenance records;
- cohort and module definitions;
- executable discovery code;
- primary and sensitivity results;
- universal deletion analyses;
- machine-readable tables;
- model diagnostics;
- validation and checksum records.

Key directories include:

    protocol/
    src/
    results/
    validation/

## Running the code

Install the packages listed in the repository-level `requirements.txt`.

Large source and intermediate files are not stored in Git.

Their locations can be supplied through the `TCGA_*` environment
variables defined in `src/common_v3.py`.

Generated intermediates default to:

    work/intermediate

Frozen public outputs default to:

    results/

Read the frozen protocol, provenance records, and repository-level
[REPRODUCIBILITY.md](../../REPRODUCIBILITY.md) before execution.

## Interpretation boundary

The analysis does not establish:

- causal transcription-factor activity;
- direct DNA binding;
- uniform fetal reactivation;
- diagnostic, prognostic, or predictive biomarker validity;
- clinical utility;
- treatment prediction;
- therapeutic dependence.

Bulk regulon scores summarise signed target-expression patterns at cohort
level.

Pointwise association, deletion robustness, and later external
transportability must remain separately reported.
