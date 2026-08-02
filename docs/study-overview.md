# Study release overview

This repository is the public computational and reproducibility record for
the study:

> Robustness and cross-cohort concordance of developmental regulons in
> endometrial carcinoma

The release contains five analysis packages covering discovery, external
evaluation, and explanatory sensitivity analyses.

## Scientific objective

The study asks whether developmental transcription-factor regulon signals
associated with endometrial carcinoma molecular classes remain interpretable
after three distinct evidential questions are separated:

1. whether a prespecified contrast shows a pointwise association;
2. whether the signal survives deletion of every mapped target gene;
3. whether the signal retains its prespecified direction in an independent
   external cohort.

These properties are reported separately. Statistical significance does not
automatically establish target-set robustness, and an external target-level
pass does not automatically establish confirmatory replication.

## Cohorts

### TCGA-UCEC

The discovery cohort contains 507 primary tumours:

- 49 POLE-ultramutated;
- 148 mismatch-repair-deficient, abbreviated MMRd;
- 147 no specific molecular profile, abbreviated NSMP;
- 163 p53-abnormal.

The primary purity-adjusted complete-case model contains 506 tumours. The
prespecified no-purity sensitivity contains 507 tumours.

### CPTAC-UCEC

The external-evaluation cohort contains 230 cases divided into two
donor-independent strata:

- Discovery, n = 95;
- Confirmatory, n = 135.

The strata are analysed separately and combined by inverse-variance
fixed-effect meta-analysis. Raw expression is not pooled.

## Analysis packages

### TCGA primary discovery

`analyses/tcga-primary-discovery/` contains the frozen discovery analysis of
the fetal Mullerian epithelial module and the prespecified 20-factor
developmental regulon panel.

It includes:

- cohort and model specifications;
- module and regulon provenance;
- molecular-subtype contrasts;
- omnibus and gated multiple-testing families;
- pointwise result categories;
- universal single-target-deletion analyses;
- result and validation records.

### CPTAC external evaluation

`analyses/cptac-external-evaluation/` contains the frozen target-level
external evaluation of six TCGA signals.

It includes:

- cohort and subtype linkage records;
- donor-independence checks;
- power and evaluability records;
- stratum-specific effects;
- fixed-effect meta-analysis;
- target-level verdict and provenance records.

### No-purity sensitivity

`analyses/no-purity-sensitivity/` evaluates the six TCGA targets after
omitting tumour purity from the model.

It includes:

- no-purity point estimates;
- complete mapped-target leave-one-out analyses;
- comparison with the primary model;
- purity and coefficient decomposition;
- deterministic and independent verification.

### Molecular-subtype sensitivity

`analyses/molecular-subtype-sensitivity/` provides the post-hoc
decomposition of the frozen C2 contrast into:

- POLE versus NSMP;
- MMRd versus NSMP;
- direct POLE versus MMRd.

The decomposition uses unchanged full models and fixed biological contrast
weights. It does not refit pairwise subsets.

### Grade and histology sensitivity

`analyses/grade-histology-sensitivity/` contains the outcome-blind
clinical-composition sensitivity analysis.

It includes:

- TCGA binary-histology adjustment;
- endometrioid-restricted binary-grade adjustment;
- matched base and adjusted models;
- attenuation decomposition;
- influence diagnostics;
- independent reconstruction and validation.

CPTAC clinical modelling is not included because the available annotations
did not satisfy the frozen harmonisation rules.

## Main evidence map

- GATA2 and SOX9 retain the prespecified negative C2 direction in TCGA and
  both CPTAC strata. Their supported interpretation is cross-cohort
  directional concordance, not confirmatory external replication, because
  neither source TCGA signal passed the universal deletion gate.
- HOXA9 and WT1 are pointwise TCGA signals that were externally evaluable
  but did not confirm.
- PAX8 and LHX1 passed universal single-target-deletion robustness in TCGA.
  Their C1 external evaluation was underpowered before outcome read; PAX8
  remains sensitivity-only and LHX1 triggered the opposite-direction veto.
- The broad fetal Mullerian epithelial module was non-confirmatory.

## Release contents

Each package may contain:

- frozen protocols and amendments;
- executable analysis code;
- source acquisition and linkage records;
- machine-readable results;
- figures and manuscript-facing tables;
- deterministic rerun records;
- independent validation records;
- checksums and provenance seals.

The repository-level `release/` directory contains the public file manifest,
SHA-256 checksums, and source mapping for the published repository state.

## What is not included

Large upstream expression datasets are not redistributed as a general-purpose
copy.

Unfinished experiments, abandoned execution attempts, private workspaces,
runtime environments, caches, and superseded intermediate cycles are outside
the public release unless they are required to explain a retained deviation
or validation result.

## Interpretation boundary

The release does not establish:

- causal transcription-factor activity;
- direct DNA binding;
- uniform fetal-program reactivation;
- biomarker validity or clinical utility;
- treatment response;
- therapeutic dependence.

Bulk regulon scores summarise signed target-expression patterns at cohort
level. The manuscript, supplementary material, frozen protocols, and
machine-readable verdict records govern the supported scientific wording.
