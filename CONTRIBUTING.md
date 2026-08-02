# Contributing

Thank you for helping improve the computational and reproducibility record
for the endometrial developmental-regulon study.

This repository is a version-controlled scientific record. It is not a
general-purpose analysis framework. Every contribution must preserve a clear
distinction between:

- correction of the frozen study record;
- improvement of documentation or presentation;
- independent reproduction of a reported analysis;
- development of a new analysis or extension.

## Suitable contributions

Appropriate contributions include:

- correction of a documented typographical or formatting error;
- clarification of repository navigation or execution instructions;
- correction of a broken internal link;
- a release-verification or reproduction report;
- identification of a provenance inconsistency;
- identification of a numerical discrepancy;
- a validation test that strengthens integrity checks;
- accessibility improvements to documentation, tables, or figures;
- a clearly separated extension that does not overwrite the frozen study.

Scientific questions and challenges are welcome when they identify the
relevant target, contrast, model, file, and evidential boundary.

## Open an issue first

Open an issue before preparing a substantial change.

The issue should identify:

- the affected analysis package;
- the affected file, result object, or figure;
- the repository commit used;
- the expected and observed behaviour;
- whether the matter concerns code, data provenance, execution, numerical
  output, presentation, or interpretation;
- whether the proposed change could alter a reported value, category, or
  verdict.

Do not post restricted data, credentials, tokens, private reviewer links, or
controlled-access identifiers in a public issue.

Privacy and security concerns should be reported privately to
`contact@helena.bio`.

## Contribution categories

Every pull request must identify one of the following categories.

### Documentation-only

Changes wording, navigation, formatting, accessibility, or explanatory
material without altering scientific inputs, executable code paths,
numerical outputs, categories, or verdicts.

### Reproducibility tooling

Changes verification, validation, packaging, or environment support without
altering the frozen scientific specification.

The pull request must identify every file added to or removed from a release
manifest.

### Scientific correction

Corrects an error in a frozen specification, input mapping, executable
analysis, numerical output, figure, category, verdict, or interpretation.

Scientific corrections require explicit scientific and provenance review.
They must not be presented as documentation-only changes.

### Independent reproduction

Adds a report from an independently executed workflow.

A reproduction report does not replace the frozen result. It must identify
the source data, environment, commands, comparison method, tolerances, and
every deviation from the frozen protocol.

### New analysis or extension

Introduces a different dataset, updated source release, alternative model,
new covariate, changed target network, new endpoint, or modified threshold.

A new analysis must be placed in a clearly separated package. It must not
overwrite frozen inputs, results, or verdict records.

## Frozen scientific record

Do not silently modify:

- frozen protocols;
- module or regulon definitions;
- transcription-factor or target inventories;
- molecular-subtype contrasts or weights;
- covariate definitions;
- cohort-construction rules;
- patient, sample, or file mappings;
- statistical thresholds;
- multiplicity families;
- deterministic seeds;
- numerical result tables;
- verdict records;
- independent validation records;
- release manifests or checksums.

When a scientific correction is necessary, preserve the prior state in Git
history and describe the consequence explicitly in the pull request and
`CHANGELOG.md`.

A post-hoc analysis must remain labelled post-hoc and cannot promote, rescue,
or replace a frozen verdict.

## Data and privacy requirements

Contributions must not include:

- controlled-access source files;
- identifiable or potentially re-identifiable participant data;
- credentials or access tokens;
- private workspace material;
- cached upstream datasets;
- files whose redistribution is not permitted;
- private reviewer links.

Permitted derived data must retain sufficient provenance to identify their
source, transformation, and analysis role.

## Requirements for new analysis code

New analysis code must include:

- a clearly stated scientific question;
- an expected input schema;
- input provenance and acquisition requirements;
- environment or dependency information;
- deterministic seed handling where relevant;
- machine-readable results;
- validation checks;
- documented output paths;
- an interpretation boundary;
- a reproducibility record.

Runtime environments, downloaded upstream datasets, caches, and private
workspace files are not accepted.

## Pull-request requirements

A pull request must state:

1. the problem being addressed;
2. the affected package and files;
3. the contribution category;
4. the scientific impact:
   - none;
   - presentation only;
   - numerical values changed;
   - category or verdict changed;
5. the validation commands and results;
6. manifest or checksum consequences;
7. source-data consequences;
8. interpretation consequences;
9. related issue references.

Do not mix unrelated changes in one pull request.

## Validation

Run the repository-level verifier:

    python3 scripts/verify-release.py

Also run the package-local checks documented by the affected analysis.

For documentation-only changes, inspect:

    git diff --check
    git diff --stat
    git status --short

A release-verification failure must be explained. Do not regenerate
checksums merely to conceal an unexplained file change.

## Commit messages

Use a scoped and factual commit message, for example:

    docs: clarify repository evidence boundaries
    fix(cptac): correct target mapping provenance
    test(release): validate manifest path uniqueness

Avoid wording that implies scientific confirmation beyond the stored
verdict.

## Authorship and attribution

A code, documentation, or review contribution does not automatically
establish manuscript authorship.

Authorship decisions remain with the manuscript authors and must follow the
applicable journal criteria. Repository contributions may be credited
through Git history, release notes, or acknowledgements as appropriate.

## Review standard

Review evaluates:

- scientific scope and claim discipline;
- provenance and numerical traceability;
- privacy and redistribution constraints;
- reproducibility;
- consistency with frozen records.

A technically correct patch may still be declined if it obscures the frozen
analysis boundary or changes the meaning of the scientific record.
