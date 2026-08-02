# Data boundary

The `data/` directory documents the boundary between upstream source data,
public acquisition and linkage records, permitted derived inputs, and
machine-readable outputs retained in the analysis packages.

It is not a general-purpose mirror of TCGA-UCEC, CPTAC-UCEC, or other
upstream repositories.

## Directory structure

    data/
    |-- acquisition/
    |   `-- tcga/
    |-- external/
    `-- README.md

Package-specific data and provenance records may also be located inside the
relevant directory under `analyses/`.

## `acquisition/`

The `acquisition/` directory contains public records needed to identify and
reconstruct source acquisition and cohort construction.

Depending on the source and redistribution terms, these records may include:

- source repository and release identifiers;
- query and download manifests;
- case, sample, aliquot, and file identifiers;
- patient-to-sample mappings;
- cohort inclusion and exclusion rules;
- source checksums;
- frozen cohort definitions;
- provenance notes and reconciliation records.

The TCGA-UCEC acquisition record is located under:

    data/acquisition/tcga/

CPTAC-UCEC acquisition and linkage materials are retained within:

    analyses/cptac-external-evaluation/data-linkage/

## `external/`

The `external/` directory documents external inputs that are referenced by
the study but are not redistributed as complete source datasets.

An external input may require independent acquisition because:

- the file is large;
- redistribution is restricted;
- access is controlled by the source repository;
- the source terms require retrieval from the original provider;
- the input is supplied through a documented environment variable during
  execution.

The package README and provenance files define how each external input is
used.

## Data retained in Git

Small files may be retained when redistribution is appropriate and they are
needed to inspect, verify, or reproduce the public record.

Examples include:

- frozen cohort rosters;
- identifier mappings;
- acquisition manifests;
- public clinical-variable extracts;
- harmonisation maps;
- source-field dictionaries;
- model-ready derived inputs where permitted;
- machine-readable result tables;
- figure source data;
- validation and checksum records.

Retention in Git does not by itself make a file authoritative. Its role is
defined by the package protocol, provenance record, and scientific-state
hierarchy.

## Data not retained in Git

The repository does not publish:

- large upstream expression matrices as a general-purpose copy;
- controlled-access or restricted source data;
- identifiable or potentially re-identifiable participant information;
- credentials, authentication tokens, or signed download URLs;
- private reviewer links;
- private workspace copies;
- temporary downloads and runtime caches;
- files whose redistribution is not permitted.

Users must obtain such material independently from the responsible source
under the applicable access and use conditions.

## Source identity

A valid source record should distinguish:

- the upstream repository;
- the release or retrieval date;
- the source identifier;
- the biological unit represented;
- the selected file;
- the transformation applied;
- the analysis package that consumes the result.

A filename alone is not sufficient provenance.

Patient, case, sample, aliquot, and file identifiers must not be treated as
interchangeable.

## Cohort construction

The frozen cohort definitions and package-local linkage records govern which
source records enter each analysis.

For TCGA-UCEC, the public record identifies the primary-tumour analytical
unit and the model-specific complete-case counts.

For CPTAC-UCEC, the public record preserves separate Discovery and
Confirmatory strata and checks for patient, sample, and aliquot overlap.

Changing source releases, subtype annotations, sample selection, or linkage
rules creates a new analysis unless the deviation is documented explicitly.

## Derived data and outputs

Derived data are generated from documented source inputs through
version-controlled transformations.

A derived file should retain enough context to identify:

- its source inputs;
- the script or workflow that produced it;
- its schema;
- the package that owns it;
- whether it is an input, intermediate, result, or validation artifact;
- whether it has frozen verdict authority or is descriptive support.

Post-hoc derived outputs cannot replace or promote a frozen verdict.

## Checksums

Repository-level published files are recorded in:

- `release/file-manifest.tsv`;
- `release/SHA256SUMS`.

Additional package-local checksum and seal files may protect source inputs,
results, or reproducibility records within a narrower scope.

A matching checksum establishes byte identity. It does not establish
biological correctness, causal validity, or clinical utility.

## Full rerun requirements

A user attempting a complete rerun should record:

- the exact repository commit or release;
- upstream source releases;
- input checksums;
- patient and sample mappings;
- dependency versions;
- transformation commands;
- generated output checksums;
- all deviations from the frozen protocol.

See [../REPRODUCIBILITY.md](../REPRODUCIBILITY.md) and
[../docs/data-provenance.md](../docs/data-provenance.md) for the full
reproducibility and provenance guidance.

## Privacy and access

Do not open a public issue or pull request containing restricted data,
credentials, private reviewer links, or potentially identifiable
participant information.

Sensitive data, privacy, or access concerns should be reported privately to
`contact@helena.bio`.
