# Endometrial developmental regulons

This repository contains the code, frozen analysis specifications,
derived results, figures, and independent validation records for three
sensitivity analyses linked to the endometrial developmental regulon study.

The release is organised by scientific question.

## Analyses

### No-purity sensitivity

The package contains the final no-purity comparison, gene leave-one-out
analysis, machine-readable results, and independent validation record.

### Molecular-subtype sensitivity

The package contains the frozen analysis specification, executable code,
target-level results, deterministic reruns, figures, and independent
validation material.

### Grade and histology sensitivity

The package contains the outcome-blind data protocol, frozen analysis plan,
harmonised inputs, final scientific outputs, deterministic run comparison,
and independent validation material.

## Repository structure

- analyses contains the three completed analysis packages.
- data describes the data boundary.
- docs describes provenance, reproducibility, and limitations.
- release contains the public file manifest and checksums.
- scripts contains repository-level verification code.

Large upstream datasets are not redistributed. Source records and checksums
are retained where they form part of the completed analysis packages.

## Verification

Run from the repository root:

    python3 scripts/verify-release.py
