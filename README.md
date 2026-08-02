# Endometrial developmental regulons

This repository is the version-controlled computational and reproducibility
record for the study:

> Robustness and cross-cohort concordance of developmental regulons in
> endometrial carcinoma

The study evaluates a fetal Mullerian epithelial module and a prespecified
panel of 20 developmental transcription-factor regulons in TCGA-UCEC,
followed by frozen target-level external evaluation in CPTAC-UCEC.

The repository contains frozen analysis specifications, executable code,
acquisition and provenance records, derived result tables, figures,
deterministic reruns, independent validation records, and release checksums.

## Scientific question

The study separates three evidential properties that are often merged:

1. **Pointwise association** - whether a prespecified molecular-subtype
   contrast meets its locked direction, effect-size, confidence-interval,
   omnibus, and multiplicity criteria.
2. **Universal single-target-deletion robustness** - whether the complete
   result survives deletion of every mapped target gene, one at a time.
3. **Cross-cohort directional concordance** - whether a target retains its
   prespecified direction and external-evaluation criteria in CPTAC-UCEC.

A statistically significant association is not automatically described as a
robust regulon. An external target-level pass is not described as
confirmatory replication when its source TCGA signal lacks universal
single-target-deletion credit.

## Study design

### TCGA-UCEC discovery

The discovery cohort contains 507 primary tumours assigned to four molecular
classes:

- POLE-ultramutated;
- mismatch-repair-deficient, abbreviated MMRd;
- no specific molecular profile, abbreviated NSMP;
- p53-abnormal.

The primary complete-case model contains 506 tumours and includes molecular
subtype, a frozen proliferation covariate, broad stromal and immune
composition terms, and tumour purity. The no-purity sensitivity contains
507 tumours.

### CPTAC-UCEC external evaluation

Six frozen target-level signals are evaluated in one independent CPTAC-UCEC
cohort represented by two donor-independent strata:

- Discovery, n = 95;
- Confirmatory, n = 135.

Stratum-specific standardized effects are combined by inverse-variance
fixed-effect meta-analysis. Raw expression is not pooled across strata.

The CPTAC analysis is model-adapted to the external cohort. It is not an
identical covariate replication of the TCGA model, and the standardized
effects are not interchangeable estimands.

### Post-hoc explanatory analyses

Separately frozen descriptive analyses examine:

- the POLE-versus-NSMP and MMRd-versus-NSMP components of C2;
- the direct POLE-versus-MMRd contrast;
- TCGA binary-histology sensitivity;
- endometrioid-restricted binary-grade sensitivity;
- matched-sample attenuation and influence diagnostics.

These analyses may explain or delimit a signal. They cannot promote, rescue,
or replace a frozen verdict.

## Molecular-subtype contrasts

The contrast weights are fixed biological weights rather than sample-size
weights.

- **C1** compares p53-abnormal tumours with equal weights across the three
  non-p53-abnormal classes.
- **C2** compares `0.5 * POLE + 0.5 * MMRd` with NSMP.
- **C3** compares POLE with MMRd and is descriptive.

Observed class proportions do not alter these weights.

## Main evidence map

The public record retains favourable, non-confirming, and underpowered
results.

- **GATA2** and **SOX9** show negative C2 effects in TCGA and retain the
  same direction in both CPTAC strata. Their supported interpretation is
  cross-cohort directional concordance, not confirmatory external
  replication, because neither TCGA source signal passed the universal
  single-target-deletion gate.
- **HOXA9** and **WT1** are pointwise TCGA signals that were externally
  evaluable but did not confirm.
- **PAX8** and **LHX1** passed universal single-target-deletion robustness
  in TCGA. Their C1 external evaluation was underpowered before the CPTAC
  expression read. PAX8 remains sensitivity-only, while LHX1 triggered the
  opposite-direction veto.
- The broad fetal Mullerian epithelial module was non-confirmatory under
  the prespecified rules.

The supported result is a differentiated evidence map of developmental
regulons. It is not a clinically validated biomarker panel or a demonstrated
fetal-reactivation mechanism.

## Repository structure

    .
    |-- analyses/
    |   |-- tcga-primary-discovery/
    |   |-- cptac-external-evaluation/
    |   |-- no-purity-sensitivity/
    |   |-- molecular-subtype-sensitivity/
    |   `-- grade-histology-sensitivity/
    |-- data/
    |-- docs/
    |-- release/
    |-- scripts/
    |-- CHANGELOG.md
    |-- CODE_OF_CONDUCT.md
    |-- CONTRIBUTING.md
    |-- README.md
    |-- REPRODUCIBILITY.md
    `-- requirements.txt

## Analysis packages

| Package | Contents |
| --- | --- |
| `analyses/tcga-primary-discovery/` | Frozen TCGA-UCEC discovery workflow, module and regulon definitions, model specification, primary results, deletion analyses, and validation records |
| `analyses/cptac-external-evaluation/` | Frozen CPTAC-UCEC target-level evaluation, linkage records, model-adapted analyses, fixed-effect meta-analysis, result tables, and provenance seals |
| `analyses/no-purity-sensitivity/` | TCGA no-purity model, complete mapped-target leave-one-out results, purity decomposition, deterministic reruns, and independent checks |
| `analyses/molecular-subtype-sensitivity/` | Post-hoc decomposition of C2 into POLE-versus-NSMP, MMRd-versus-NSMP, and direct POLE-versus-MMRd estimates |
| `analyses/grade-histology-sensitivity/` | Outcome-blind TCGA histology and endometrioid-grade sensitivity, matched analyses, influence diagnostics, and independent verification |

Each package contains its own README and package-local protocol, source,
result, or validation materials.

## Supporting directories

| Directory | Purpose |
| --- | --- |
| `data/` | Source acquisition records, identifiers, mappings, checksums, and descriptions of the public data boundary |
| `docs/` | Study overview, provenance, limitations, path mapping, and reproducibility guidance |
| `release/` | Public file manifest, SHA-256 checksums, and source mapping |
| `scripts/` | Repository-level release verification |

## Verify the published files

From the repository root:

    python3 scripts/verify-release.py

The verifier checks every path listed in `release/file-manifest.tsv` against
its recorded size and SHA-256 digest.

A successful run confirms that the local files are complete and
byte-consistent with the published manifest. It does not download source
data, reconstruct cohorts, rerun models, or establish scientific validity.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the full reproducibility
boundary.

## Source data

Large upstream expression datasets are not redistributed as a
general-purpose copy in this repository.

A complete rerun requires acquisition of the identified TCGA-UCEC,
CPTAC-UCEC, and reference inputs from their original repositories, subject
to the source terms and access controls.

The public acquisition records and package-local provenance files identify
the releases, files, mappings, and cohort-construction rules used in the
reported analysis.

Restricted or controlled-access material must be obtained independently by
an authorised user.

## Interpretation boundary

The repository does not establish:

- causal transcription-factor activity;
- direct DNA binding from bulk-expression associations;
- uniform fetal-program reactivation;
- diagnostic, prognostic, or predictive biomarker validity;
- clinical utility;
- treatment response;
- therapeutic dependence;
- equality or equivalence merely because an interval spans zero;
- absence of confounding merely because one sensitivity analysis is stable.

Bulk regulon scores summarise signed target-expression patterns at cohort
level. They are not direct measurements of transcription-factor protein
activity or binding.

## Contributing

Corrections, reproduction reports, documentation improvements, validation
improvements, and clearly separated extensions are welcome.

Changes to frozen specifications, inputs, numerical results, categories,
figures, or verdict records require explicit scientific and provenance
review.

Read [CONTRIBUTING.md](CONTRIBUTING.md) before opening a pull request.

## Citation and release status

The manuscript and immutable citation metadata are being prepared.

Until a versioned release and persistent identifier are published, identify
the repository by organisation, repository name, exact commit, and access
date. Do not cite the moving `main` branch as though it were an immutable
scientific release.

## Maintainer

Maintained by [Helena Bioinformatics](https://helena.bio).

Scientific and reproducibility questions may be submitted through GitHub
Issues. Privacy, restricted-data, and security concerns should be reported
privately to `contact@helena.bio`.
