# Support

This repository is a scientific computational and reproducibility record.
Support is organised by the type of question.

## Scientific and methodological questions

Open a GitHub issue when asking about:

- a frozen analysis specification;
- a molecular-subtype contrast;
- a target-level result or verdict;
- interpretation boundaries;
- a possible numerical discrepancy;
- an independent reproduction.

Identify the repository commit, analysis package, file, target, contrast,
and model involved.

Use the reproduction-report template for an independently executed run.

## Scientific corrections

Use the scientific-correction template when reporting a possible error in:

- a frozen protocol;
- an input mapping;
- executable analysis code;
- a numerical result;
- a figure or table;
- a category or verdict;
- manuscript-facing interpretation.

A correction report should explain whether the issue could alter numerical
values, categories, verdicts, figures, or supported wording.

## Repository and execution questions

Open a standard GitHub issue for:

- repository navigation;
- environment setup;
- release verification;
- missing public files;
- broken internal links;
- documentation problems.

Before opening an issue, run:

    python3 scripts/verify-release.py

Include the exact command and complete error message when reporting an
execution problem.

## Sensitive matters

Do not use a public issue for credentials, restricted data, participant
information, private reviewer links, or security vulnerabilities.

Report those privately to:

`contact@helena.bio`

See [SECURITY.md](SECURITY.md) for the security-reporting process.

## What support does not provide

Repository support does not include:

- access to controlled upstream datasets;
- patient-specific interpretation;
- clinical diagnosis or treatment advice;
- validation of the materials for clinical use;
- guarantees that modified inputs reproduce the frozen study;
- support for undocumented forks or altered analysis specifications.

## Maintainer

The repository is maintained by Helena Bioinformatics.
