## Summary

Describe the problem and the change.

## Contribution category

- [ ] Documentation-only
- [ ] Reproducibility tooling
- [ ] Scientific correction
- [ ] Independent reproduction
- [ ] New analysis or extension

## Affected packages and files

List every analysis package and file in scope.

## Scientific impact

- [ ] No scientific impact
- [ ] Presentation only
- [ ] Numerical values changed
- [ ] Result category changed
- [ ] Frozen verdict changed
- [ ] Figure or manuscript wording changed

Explain any selected scientific impact:

## Validation

List the commands run and their results.

    python3 scripts/verify-release.py

## Data and provenance

State whether source data, identifiers, mappings, checksums, or acquisition
records changed.

## Release metadata

State whether `release/file-manifest.tsv` or `release/SHA256SUMS` changed
and why.

## Interpretation boundary

Confirm that the change does not silently promote, rescue, or replace a
frozen verdict.

## Checklist

- [ ] I read `CONTRIBUTING.md`.
- [ ] I identified the exact repository commit used.
- [ ] I did not include restricted or identifiable data.
- [ ] I documented all deviations from the frozen protocol.
- [ ] I ran the repository verifier.
- [ ] I inspected the complete diff.
- [ ] I updated release metadata for every published file changed or added.
