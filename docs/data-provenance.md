# Data provenance

This document describes how source identity, cohort construction, file
selection, transformations, and derived outputs are recorded in the public
repository.

The repository does not treat a filename or a downloaded archive as sufficient
provenance. A reproducible analytical input requires a traceable relationship
between the upstream source, the selected file, the biological unit, the
transformation applied, and the result package that consumes it.

## Provenance principles

The public record follows these principles:

- source repositories and releases are named explicitly;
- source identifiers are retained where redistribution permits;
- patient, case, sample, aliquot, and file mappings are kept distinct;
- one primary-tumour analytical unit is retained per patient where specified;
- frozen cohort and linkage rules are recorded before outcome-bearing reads;
- file hashes and release checksums are retained for load-bearing inputs;
- transformations are represented by version-controlled code or documented
  derivation records;
- generated outputs are linked to the package and protocol that produced them;
- substituted or newer inputs are treated as a new analysis unless the
  deviation is documented explicitly.

## TCGA-UCEC source record

The TCGA discovery analysis uses primary tumour transcriptomic and clinical
records acquired through the Genomic Data Commons.

The public repository retains, where applicable:

- the frozen cohort definition;
- source release and query records;
- case, sample, and file identifiers;
- primary-tumour selection rules;
- patient-to-sample mappings;
- molecular-subtype assignments used by the analysis;
- tumour-purity and composition-source provenance;
- input checksums and reconciliation records;
- the gene-universe and zero-gene handling rules.

The discovery cohort contains 507 primary tumours. The primary
purity-adjusted model contains 506 complete cases, while the no-purity model
contains 507.

Repository-level TCGA acquisition material is located under:

    data/acquisition/tcga/

Additional package-local provenance records are retained under:

    analyses/tcga-primary-discovery/

## CPTAC-UCEC source record

The external evaluation uses one non-TCGA CPTAC-UCEC cohort represented by
two donor-independent strata:

- Discovery, n = 95;
- Confirmatory, n = 135.

The public linkage record distinguishes:

- patient-level molecular-subtype labels;
- GDC cases;
- primary-tumour samples;
- aliquots;
- RNA files;
- stratum membership.

The record includes overlap checks confirming that no patient, sample, or
aliquot is shared between the two strata.

CPTAC linkage and acquisition materials are retained under:

    analyses/cptac-external-evaluation/data-linkage/

The source transcriptomic files are not pooled across strata. Frozen scoring
and model rules are applied within each stratum before fixed-effect
meta-analysis.

## Developmental reference and regulon provenance

The broad fetal Mullerian epithelial module and the developmental regulon
panel have separate provenance chains.

The fetal reference is tied to its named upstream developmental
reproductive-tract dataset and the frozen module-definition record.

The developmental regulons use the frozen 20-factor inventory and signed
CollecTRI consensus network recorded before tumour outcome analysis.

The repository preserves:

- module provenance and source annotation;
- factor inventory provenance;
- mapped target and edge records where included in the frozen package;
- scoring rules;
- exclusions and composition-overlap rules;
- amendments to the frozen specification.

No missing edge is fabricated, and no factor is added or reprioritised from
the observed TCGA or CPTAC effects.

## Clinical-variable provenance

Grade and histology variables were subjected to an outcome-blind feasibility
and harmonisation audit before sensitivity modelling.

The repository retains:

- raw acquisition responses;
- canonicalised copies;
- source-field dictionaries;
- category counts;
- linkage ledgers;
- missingness summaries;
- harmonisation maps;
- outcome-firewall attestations;
- independent source and semantic audits.

The frozen rules permitted TCGA binary-histology adjustment and
endometrioid-restricted binary-grade adjustment.

CPTAC clinical variables were not modelled because their available
annotations did not satisfy the frozen harmonisation criteria.

## Derived data

Permitted derived data may include:

- frozen cohort rosters;
- identifier mappings;
- module and regulon scores;
- model-ready matrices where redistribution is appropriate;
- machine-readable result tables;
- model diagnostics;
- meta-analysis tables;
- sensitivity summaries;
- figure source data;
- validation and checksum records.

Derived data are not treated as self-explanatory. Their package location,
schema, source relationship, and interpretation boundary must be documented.

## Data not redistributed

The repository does not redistribute large upstream datasets as a
general-purpose mirror.

It also does not publish:

- controlled-access or restricted source material;
- identifiable or potentially re-identifiable participant information;
- credentials or authenticated download material;
- private reviewer links;
- private workspace copies;
- runtime caches and temporary downloads;
- files whose redistribution terms do not permit publication.

Users must obtain externally hosted or restricted data from the responsible
source under the applicable access conditions.

## Checksums and release records

The repository-level public release record contains:

- `release/file-manifest.tsv` - path, size, and SHA-256 for each published
  file;
- `release/SHA256SUMS` - checksum list used by the release verifier;
- `release/source-mapping.tsv` - mapping between public paths and retained
  source-package records.

Package-local checksum, seal, and reproducibility files provide additional
scope-specific records.

Matching a checksum establishes byte identity for that file. It does not
establish biological correctness, causal validity, or clinical utility.

## Path normalisation

Host-specific private paths are not part of the scientific identity of an
input.

Public text records use normalised repository-relative or documented external
paths. The source mapping preserves the relationship to retained package
material without exposing private workspace locations.

See [path-mapping.md](path-mapping.md) for the repository path policy.

## Reuse and deviation

A user re-executing the analysis should record:

- the exact source release;
- source-file checksums;
- patient and sample mappings;
- any missing or substituted input;
- dependency versions;
- transformation commands;
- output checksums;
- every deviation from the frozen protocol.

Using a newer GDC release, changed annotation, alternative sample selection,
different target network, or substituted clinical variable creates a new
analysis. It must not be described as an exact reproduction of the frozen
study without a complete deviation record.
