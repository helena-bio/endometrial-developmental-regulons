# Reproducibility

This repository accompanies the endometrial developmental-regulon study. It contains the analysis code, frozen model and module definitions, acquisition records, result tables, figures, sensitivity analyses, independent audit records, and release checksums used for the manuscript.

The repository supports two distinct operations:

1. verification of the published release and its recorded results;
2. re-execution of analyses after the required source data have been acquired.

These operations should not be treated as equivalent. A checksum verifies that a published file has not changed. It does not establish scientific correctness or recreate an analysis from raw data.

## Verify the release

From the repository root, run:

```bash
python3 scripts/verify-release.py
```

The verifier checks the files listed in `release/file-manifest.tsv` against their recorded sizes and SHA-256 values in `release/SHA256SUMS`. A successful run confirms that the checked-out release is complete and byte-consistent with the published manifest.

This check does not download data, fit models, or reassess the scientific interpretation.

## Python environment

The public scripts use Python and the packages declared in `requirements.txt`. A local environment can be created with:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The frozen analysis records also identify the software versions and execution conditions used for the reported runs. Those records are retained within the corresponding analysis packages. Installing the current dependency file is sufficient for inspecting and running the public Python code, but it is not a claim of byte-identical numerical output across all operating systems, hardware, BLAS implementations, or future package releases.

## Analysis packages

The principal workflows are organised under `analyses/`:

- `tcga-primary-discovery/` contains the frozen TCGA-UCEC discovery analysis, module and regulon definitions, model specification, source code, results, and validation records.
- `cptac-external-evaluation/` contains the target-level CPTAC-UCEC external evaluation, cohort records, analysis rules, results, and provenance seals.
- `no-purity-sensitivity/` contains the TCGA no-purity sensitivity analysis, gene leave-one-out assessment, deterministic rerun records, and independent numerical audits.
- `molecular-subtype-sensitivity/` contains the post-hoc per-class decomposition.
- `grade-histology-sensitivity/` contains the TCGA grade and histology sensitivity analysis.

Each package should be read from its own `README.md`, protocol files, and `src/` directory. The package-local documentation defines the expected inputs, analysis boundary, output files, and interpretation limits. Scripts should not be moved between packages or run against substituted inputs without recording that change.

## Source data

Raw TCGA-UCEC and CPTAC-UCEC files are not redistributed as a general-purpose copy inside this repository. Their acquisition records, identifiers, release information, mappings, checksums, and cohort-construction records are published where redistribution permits.

A full rerun therefore requires the user to acquire the relevant files from the Genomic Data Commons or the named upstream source, subject to the source terms and access controls. The acquisition records under `data/acquisition/` and the package-local provenance files define the inputs used in the study.

Restricted or controlled-access material must be obtained by an authorised user. This repository does not bypass source access requirements.

## Re-executing an analysis

Before running a workflow:

1. verify the repository release;
2. read the package README and frozen protocol;
3. acquire the exact source release identified by the acquisition record;
4. confirm input checksums and sample mappings;
5. install the declared Python dependencies;
6. run the package-local script in the documented order;
7. compare generated outputs with the published result and audit tables.

The public result tables remain the manuscript record. A rerun using a newer GDC release, changed clinical annotations, altered dependencies, or a different cohort definition is a new analysis and should not be presented as a reproduction of the frozen study.

## Independent checks

The repository includes independent reconstruction and comparison records for the load-bearing analyses. These records cover family-wise multiple-testing calculations, the six prespecified targets, complete gene leave-one-out credit, purity decomposition, deterministic reruns, and scientific-byte preservation.

They document the checks performed for the reported release. Their presence does not convert target-level directional concordance into causal evidence, clinical validation, therapeutic prediction, or proof of direct transcription-factor binding.

## Interpretation boundary

The reported analyses support pointwise associations, prespecified robustness checks, and target-level cross-cohort directional comparisons within the documented models.

The repository does not claim:

- causal transcription-factor activity;
- direct DNA binding from bulk-expression associations;
- fetal-program reactivation as an established mechanism;
- biomarker validity or clinical utility;
- treatment response or therapeutic dependence;
- equivalence between molecular classes when an interval spans zero.

The manuscript, supplementary material, and frozen verdict records govern the scientific wording. Repository files should be interpreted within those boundaries.

## Reporting a reproduction

A reproduction report should identify the repository commit, data releases, input checksums, dependency versions, operating environment, command sequence, and every deviation from the frozen protocol.

A successful release-integrity check may be reported as verification of the published files. It should not be described as a complete rerun unless the source data were reacquired and the analysis was executed from its documented inputs.
