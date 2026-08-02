# Reproducibility

This document provides a compact map of the repository reproducibility
model. The complete operational guidance is in the repository-level
[REPRODUCIBILITY.md](../REPRODUCIBILITY.md).

## Layers of reproducibility

The public release distinguishes five layers:

1. **Frozen specification** - the prespecified cohort, model, contrast,
   threshold, multiplicity, seed, and verdict rules.
2. **Executable implementation** - the version-controlled code that applies
   those rules.
3. **Recorded output** - machine-readable results, tables, figures, and
   diagnostics produced by the reported run.
4. **Deterministic rerun comparison** - repeated execution under the
   documented inputs and seeds.
5. **Independent reconstruction** - a separately instantiated comparison or
   rebuild of load-bearing calculations.

These layers answer different questions. None should be used as a substitute
for another.

## Release verification

Run from the repository root:

    python3 scripts/verify-release.py

The verifier checks the paths in `release/file-manifest.tsv` against their
recorded sizes and SHA-256 values.

A successful result establishes that the checked-out files match the
published release record.

It does not:

- acquire source data;
- reconstruct cohorts;
- execute analysis code;
- compare numerical outputs;
- establish biological correctness;
- establish clinical validity.

## Full re-execution

A complete rerun requires:

- the exact repository commit or immutable release;
- the named upstream source releases;
- source-file checksums;
- patient, sample, aliquot, and file mappings;
- the documented software environment;
- the package-local execution order;
- the frozen seeds and analysis rules;
- comparison with the published result objects;
- an explicit record of every deviation.

The repository does not redistribute every large upstream data object.
Users must acquire those files independently under the applicable source
terms and access controls.

## Package-local authority

Each analysis package defines its own execution boundary.

Before running a package, read:

- its `README.md`;
- its frozen protocol or analysis plan;
- its acquisition and provenance records;
- its source code;
- its result schema;
- its validation records.

Scripts should not be moved between packages or run against substituted
inputs without recording the change.

## Numerical comparison

A reproduction report must state what kind of agreement was assessed:

- exact byte identity;
- exact parsed-value identity;
- agreement within a declared numerical tolerance;
- agreement of the scientific category or verdict.

Floating-point results may depend on the operating system, hardware, Python
build, numerical libraries, and dependency versions.

A numerically close result is not the same as a byte-identical result, and a
matching category is not the same as exact numerical reconstruction.

## Independent checks

The retained validation records cover load-bearing parts of the study,
including:

- TCGA omnibus and gated multiple-testing families;
- universal mapped-target deletion categories;
- the six frozen CPTAC targets;
- no-purity effect decomposition;
- per-class C2 contrasts;
- CPTAC fixed-effect meta-analysis;
- grade and histology sensitivities;
- scientific-byte preservation;
- release-level checksum integrity.

These checks document numerical and provenance integrity within their stated
scope.

They do not establish:

- causal transcription-factor activity;
- direct DNA binding;
- biomarker validity;
- clinical utility;
- therapeutic prediction;
- absence of unmeasured confounding.

## Frozen and post-hoc analyses

Frozen verdicts and post-hoc explanatory analyses have different roles.

A post-hoc analysis may:

- decompose a pooled contrast;
- test sensitivity to a clinical composition variable;
- describe magnitude heterogeneity;
- identify a plausible explanation that needs further study.

It may not:

- promote a pointwise signal to universal robustness;
- rescue a failed external result;
- replace a frozen target-level verdict;
- establish equality because an interval spans zero;
- establish causality.

## New analysis boundary

A run becomes a new analysis when it changes any load-bearing element,
including:

- source-data release;
- cohort construction;
- sample mapping;
- target network;
- scoring method;
- contrast weight;
- covariate set;
- effect threshold;
- multiplicity family;
- endpoint definition;
- clinical-variable harmonisation rule.

Such a run may be scientifically valuable, but it must not be presented as
an exact reproduction of the frozen study.

## Reporting standard

A reproduction report should identify:

- repository commit or release;
- source releases and checksums;
- patient and sample mappings;
- operating system and hardware;
- dependency versions;
- commands executed;
- output checksums;
- numerical comparison method;
- tolerances;
- deviations from the frozen protocol.

A release-integrity PASS should be described as verification of published
files, not as a complete analysis rerun.

## Interpretation boundary

Matching code, inputs, checksums, and numerical outputs establishes
reproducibility within the documented computational scope.

It does not by itself establish:

- biological truth;
- causal validity;
- clinical validity;
- fitness for diagnosis, prognosis, or treatment selection;
- generalisability to another cohort, platform, or patient population.

The manuscript, supplementary material, frozen protocols, and
machine-readable verdict records govern the scientific interpretation.
