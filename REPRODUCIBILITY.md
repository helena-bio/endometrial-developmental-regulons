# Reproducibility

This repository accompanies the study:

> Robustness and cross-cohort concordance of developmental regulons in
> endometrial carcinoma

It contains frozen analysis specifications, acquisition records, executable
code, derived results, figures, sensitivity analyses, independent
validation records, and release checksums.

## Reproducibility has two distinct meanings here

The repository supports two different operations:

1. verification of the published repository files;
2. re-execution of an analysis from its documented source inputs.

These operations are not equivalent.

A checksum establishes that a recorded file has not changed. It does not
reconstruct a study cohort, rerun a statistical model, establish
scientific correctness, or prove that a different computing environment
will produce byte-identical numerical output.

## Verify the published release

From the repository root, run:

    python3 scripts/verify-release.py

The verifier reads `release/file-manifest.tsv` and checks every listed path
against its recorded file size and SHA-256 digest.

A successful run confirms that the local checkout is complete and
byte-consistent with the release manifest.

It does not:

- download upstream source data;
- reconstruct TCGA-UCEC or CPTAC-UCEC cohorts;
- calculate module or regulon scores;
- fit statistical models;
- repeat bootstrap or permutation procedures;
- reassess the scientific interpretation;
- establish causal or clinical validity.

## Python environment

Create a local environment with:

    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install --upgrade pip
    python -m pip install -r requirements.txt

Package-local environment and provenance records identify the software
conditions used for the reported analyses.

Installing the current repository requirements supports inspection and
execution of the public Python code. It is not a guarantee of
byte-identical output across operating systems, processors, BLAS
implementations, Python builds, or future dependency releases.

## Analysis packages

The principal workflows are organised under `analyses/`.

### TCGA primary discovery

`analyses/tcga-primary-discovery/` contains:

- the frozen TCGA-UCEC cohort and model specification;
- the fetal Mullerian epithelial module definition;
- the prespecified 20-factor developmental regulon panel;
- molecular-subtype contrasts and covariates;
- omnibus and gated multiplicity families;
- pointwise association results;
- universal single-target-deletion analyses;
- primary result and validation records.

Pointwise subtype association and universal deletion robustness are
reported separately.

### CPTAC external evaluation

`analyses/cptac-external-evaluation/` contains:

- Discovery and Confirmatory cohort-linkage records;
- patient, sample, and RNA-file mappings;
- the frozen six-target external-evaluation family;
- model-adapted stratum analyses;
- fixed-effect meta-analysis outputs;
- power and evaluability records;
- target-level verdict and provenance records.

The CPTAC analysis uses the same biological contrast definitions but is not
an identical covariate replication of the TCGA model.

### No-purity sensitivity

`analyses/no-purity-sensitivity/` contains:

- the TCGA no-purity model;
- comparison with the primary purity-adjusted model;
- complete mapped-target leave-one-out analyses;
- effect decomposition;
- deterministic rerun records;
- independent numerical checks.

A stable no-purity result does not establish purity independence.

### Molecular-subtype sensitivity

`analyses/molecular-subtype-sensitivity/` contains the post-hoc
decomposition of the frozen C2 contrast into:

- POLE versus NSMP;
- MMRd versus NSMP;
- POLE versus MMRd.

The estimates are derived from unchanged full models. Pairwise subsets are
not refitted.

These analyses are descriptive and cannot alter the frozen TCGA or CPTAC
verdicts.

### Grade and histology sensitivity

`analyses/grade-histology-sensitivity/` contains:

- the outcome-blind clinical-variable feasibility audit;
- TCGA binary-histology adjustment;
- endometrioid-restricted binary-grade adjustment;
- matched base and adjusted models;
- attenuation decomposition;
- influence diagnostics;
- deterministic and independent verification records.

CPTAC clinical modelling is not included because the available annotations
did not satisfy the frozen harmonisation rules.

Each package must be read together with its own README, protocol files,
source code, result records, and validation material.

## Source data boundary

Large upstream TCGA-UCEC and CPTAC-UCEC expression datasets are not
redistributed as a general-purpose copy in this repository.

Where redistribution permits, the repository retains:

- acquisition records;
- source-release identifiers;
- file and case identifiers;
- patient-to-sample mappings;
- checksums;
- cohort-construction records;
- permitted derived inputs and outputs.

A complete rerun requires the user to acquire the relevant source files
from the Genomic Data Commons or the named upstream repository, subject to
the original access terms and repository policies.

Restricted or controlled-access material must be obtained by an authorised
user. This repository does not bypass access requirements.

## Re-executing an analysis

Before executing a package:

1. record the repository commit or immutable release identifier;
2. run the repository release verifier;
3. read the package README and frozen protocol;
4. acquire the exact source release named in the acquisition record;
5. confirm input checksums and sample mappings;
6. create the documented software environment;
7. execute package-local scripts in the documented order;
8. preserve generated logs and environment information;
9. compare generated outputs with the published result records;
10. record every deviation from the frozen protocol.

The published machine-readable results and verdict records remain the
scientific record for the manuscript release.

A rerun using any of the following is a new analysis:

- a newer or different source-data release;
- changed clinical annotations;
- a modified cohort definition;
- substituted sample mappings;
- a different target network;
- changed contrast weights;
- changed covariates;
- changed thresholds or multiplicity families;
- a different outcome definition;
- an undocumented dependency substitution.

Such a run must not be represented as reproduction of the frozen study
without a complete account of the deviation.

## Determinism and numerical comparison

The reported workflows use prespecified seeds and version-controlled code
for bootstrap, permutation, reconstruction, and validation procedures.

Deterministic execution does not imply that every future platform will
produce byte-identical floating-point output.

A reproduction report must distinguish among:

- exact file identity;
- exact parsed-value identity;
- agreement within a declared numerical tolerance;
- agreement of a result category or verdict.

The comparison method and tolerance must be stated.

## Independent verification

The repository retains separate reconstruction and comparison records for
load-bearing analyses.

Their documented scope includes:

- TCGA omnibus and gated multiple-testing families;
- the six frozen external-evaluation targets;
- complete mapped-target deletion categories;
- no-purity effect decomposition;
- per-class C2 contrasts;
- CPTAC fixed-effect meta-analysis;
- grade and histology sensitivities;
- scientific-byte preservation;
- release-level checksums.

These records establish numerical and provenance integrity within their
documented scope. They do not establish causal validity.

## Scientific-state hierarchy

Repository materials have different roles and authority:

1. frozen protocols define the prespecified analysis;
2. machine-readable result and verdict records define the stored numerical
   outcome;
3. independent validation records document reconstruction or comparison;
4. the manuscript and supplementary material define the supported public
   interpretation;
5. post-hoc explanatory analyses may delimit or explain a result but cannot
   promote or replace a frozen verdict.

Presentation files must not silently override machine-readable scientific
records.

## Interpretation boundary

The repository supports pointwise associations, prespecified robustness
checks, clinical-composition sensitivities, and target-level cross-cohort
directional comparisons within the documented models.

It does not establish:

- causal transcription-factor activity;
- direct DNA binding from bulk-expression associations;
- uniform fetal-program reactivation;
- diagnostic, prognostic, or predictive biomarker validity;
- clinical utility;
- treatment response;
- therapeutic dependence;
- equality or equivalence merely because an interval spans zero;
- interchangeability of TCGA and CPTAC standardized effects;
- absence of confounding because one sensitivity analysis is stable.

Bulk regulon scores summarise signed target-expression patterns at cohort
level. They are not direct measurements of transcription-factor protein
activity or binding.

GATA2 and SOX9 are interpreted as showing cross-cohort directional
concordance rather than confirmatory replication because their source TCGA
signals did not receive universal single-target-deletion credit.

## Reporting a reproduction

A reproduction report should identify:

- the repository commit or release;
- source-data releases;
- input checksums;
- patient and sample mappings;
- dependency versions;
- operating system and hardware context;
- the command sequence;
- generated output checksums;
- the numerical comparison method and tolerance;
- every deviation from the frozen protocol.

A successful integrity check may be reported as verification of the
published repository files.

It must not be described as a complete rerun unless the source data were
reacquired, the cohort was reconstructed, and the analysis was executed
from its documented inputs.
