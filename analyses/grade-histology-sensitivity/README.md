# Grade and histology sensitivity

This package contains the outcome-blind clinical-composition sensitivity
analysis for GATA2 and SOX9 in TCGA-UCEC.

It evaluates whether the reported negative C2 effects are materially
attenuated after permitted adjustment for histology or grade.

## Outcome-blind design

Clinical-variable availability, linkage, missingness, and category support
were audited before molecular outcomes were read.

The frozen feasibility rules permitted:

- TCGA binary-histology adjustment;
- endometrioid-restricted binary-grade adjustment.

The rules excluded:

- CPTAC clinical modelling;
- all-histology grade adjustment;
- joint grade-plus-histology models;
- grade-stratified models;
- interaction models;
- rescue models.

These exclusions were fixed before the outcome-bearing analysis.

## Cohort counts

The TCGA clinical record contains:

- 382 endometrioid tumours;
- 124 non-endometrioid tumours;
- 1 tumour with missing histology;
- 207 low-grade tumours;
- 291 high-grade tumours;
- 9 tumours with missing grade.

The matched primary analyses contain:

- 505 tumours for histology sensitivity;
- 380 endometrioid tumours for grade sensitivity.

The corresponding no-purity counts are 506 and 381.

## Main result

Binary histology adjustment produces only small same-direction increases
in the GATA2 and SOX9 C2 magnitudes.

The paired delta-d confidence intervals span zero. Histology does not
attenuate the TCGA effects in the permitted matched analysis.

Within endometrioid carcinomas:

- GATA2 is attenuated by approximately 17 to 18 percent;
- SOX9 is attenuated by approximately 9 to 10 percent.

Both retain the negative direction.

Neither meets the frozen joint material-attenuation rule, which requires
both:

- absolute delta-d at least 0.10;
- attenuation of at least 20 percent.

These analyses do not change any frozen TCGA or CPTAC verdict.

## Package contents

The package records both data preparation and final analysis.

It includes:

- acquisition and canonicalisation records;
- source-field dictionaries;
- category and missingness summaries;
- linkage and harmonisation ledgers;
- outcome-firewall attestations;
- the frozen Phase 2 analysis charter;
- matched model results;
- attenuation decompositions;
- influence diagnostics;
- manuscript-ready sensitivity tables;
- two deterministic execution runs;
- independent reconstruction;
- provenance, preservation, and adversarial audits;
- validation verdicts and checksum manifests.

Key directories include:

    data-preparation/
    protocol/
    final-analysis/
    validation/

## Interpretation boundary

The supported statements are descriptive TCGA sensitivity statements.

The package does not establish:

- causal adjustment;
- that grade is purely a confounder;
- that histology or grade has no biological role;
- CPTAC clinical independence;
- absence of residual clinical-composition confounding;
- biomarker validity or clinical utility.

Grade may be a correlate of differentiation state, a confounder, or partly
downstream of molecular subtype. The adjusted model is not a uniquely
correct causal estimand.

## Reproduction

Read the outcome-blind protocol, frozen Phase 2 charter, and
repository-level [REPRODUCIBILITY.md](../../REPRODUCIBILITY.md) before
execution.

A rerun must preserve the documented source acquisition, harmonisation
rules, matched-sample construction, model specification, seeds, and
material-attenuation criterion.
