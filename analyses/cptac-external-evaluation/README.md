# CPTAC external evaluation

This package contains the frozen target-level external evaluation for the
study:

> Robustness and cross-cohort concordance of developmental regulons in
> endometrial carcinoma

Six TCGA-derived target-level signals are evaluated in one independent
CPTAC-UCEC cohort represented by two donor-independent strata.

## Scientific objective

The package asks whether the prespecified target-level direction and
external criteria are retained in CPTAC-UCEC without retuning targets,
contrast definitions, thresholds, multiplicity rules, or verdict logic
after the outcome-bearing expression read.

The external analysis is model-adapted to CPTAC rather than an identical
covariate replication of TCGA.

## Cohort structure

The external cohort contains:

- Discovery, n = 95;
- Confirmatory, n = 135.

The public linkage record maps:

- patient-level molecular-subtype labels;
- GDC cases;
- primary-tumour samples;
- aliquots;
- RNA files;
- stratum membership.

No patient, sample, or aliquot overlaps between the two strata.

Raw expression is analysed separately within each stratum and is not
pooled.

## Frozen target family

The external target family contains:

- GATA2, C2;
- SOX9, C2;
- HOXA9, C2;
- WT1, C2;
- PAX8, C1;
- LHX1, C1.

The target family and verdict rules were frozen before the
outcome-bearing CPTAC expression read.

## External model

The CPTAC model includes:

- molecular subtype;
- the frozen M4 proliferation covariate;
- ESTIMATE stromal score;
- ESTIMATE immune score.

Tumour purity is omitted because no methodologically comparable
RNA-independent open purity estimate is available for the primary
external model, and an ESTIMATE-derived surrogate would be redundant
with the composition terms.

The model uses the same C1 and C2 biological contrast definitions as the
TCGA analysis.

## Meta-analysis

Stratum-specific standardized effects are combined by inverse-variance
fixed-effect meta-analysis.

A frozen CPTAC target-level pass requires:

- the prespecified direction;
- absolute meta-analytic standardized effect at least 0.50;
- a meta-analytic confidence interval excluding zero;
- BH-adjusted q at most 0.05 across the six-target family;
- consistent direction in both strata.

An opposite stratum direction triggers a veto.

## Evaluability and power

C2 was classified as confirmatory-evaluable before the CPTAC expression
read.

C1 and external M1 equivalence testing were classified as underpowered.

Underpowered questions remain sensitivity or unresolved evidence even when
their point estimates appear favourable.

## Main result

GATA2 retains the prespecified negative C2 direction:

- Discovery d = -0.689;
- Confirmatory d = -0.999;
- fixed-effect meta d = -0.877;
- 95 percent confidence interval -1.31 to -0.45;
- BH q = 0.00036.

SOX9 also retains the prespecified negative C2 direction:

- Discovery d = -0.993;
- Confirmatory d = -0.641;
- fixed-effect meta d = -0.789;
- 95 percent confidence interval -1.25 to -0.32;
- BH q = 0.0018.

Both meet the frozen CPTAC target-level conditions.

Because the source TCGA signals did not receive universal
single-target-deletion credit, the supported manuscript interpretation is
cross-cohort directional concordance rather than confirmatory external
replication of a robustness-credited discovery.

Other target outcomes are:

- HOXA9 is evaluable but does not confirm;
- WT1 is evaluable but does not confirm;
- PAX8 remains sensitivity-only because C1 was underpowered;
- LHX1 triggers the opposite-direction veto.

## Package contents

The package includes:

- cohort-independence and subtype-provenance records;
- GDC query and download manifests;
- patient, sample, aliquot, and RNA-file linkage records;
- frozen analysis and power plans;
- acquisition and matrix quality-control code;
- scoring, modelling, sensitivity, and table-freezing code;
- stratum-specific and fixed-effect result tables;
- model diagnostics and sample-flow ledgers;
- result seals and reproducibility manifests.

Key directories include:

    data-linkage/
    protocol/
    src/
    results/
    reproducibility/

## Running the code

Install the packages listed in the repository-level `requirements.txt`.

Raw STAR-count files and generated matrices are not stored in Git.

Default working locations are under:

    work/

They may be relocated with:

- `CPTAC_WORK_DIR`;
- `CPTAC_INTERMEDIATE_DIR`;
- `CPTAC_STAR_COUNTS_DIR`.

Public linkage records are read from:

    data-linkage/

Generated result tables default to:

    results/

The external ESTIMATE gene set and the documented Dou et al.
supplementary workbook must be supplied separately through the package
environment variables.

Read the frozen protocol and repository-level
[REPRODUCIBILITY.md](../../REPRODUCIBILITY.md) before execution.

## Interpretation boundary

The package does not establish:

- equality of TCGA and CPTAC estimands;
- purity-adjusted external replication;
- causal transcription-factor activity;
- direct DNA binding;
- biomarker validity;
- clinical utility;
- treatment prediction;
- therapeutic dependence.

Directional concordance is a target-level cross-cohort result within the
documented model and verdict framework.
