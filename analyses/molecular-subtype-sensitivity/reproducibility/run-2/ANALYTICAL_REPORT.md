# TASK A analytical report

Status: COMPLETE producer execution; independent reviewer verdict pending.

This is a post-hoc explanatory sensitivity. It is descriptive and non-causal. It cannot change any frozen TCGA/CPTAC category or replication verdict, and no manuscript byte was edited.

## Sample guards

All four exact analytic counts passed: TCGA primary 506, TCGA no-purity 507, CPTAC Discovery 95, and CPTAC Confirmatory 135. Exact subtype counts are in `results/SAMPLE_COUNTS.tsv`.

## GATA2 and SOX9

The mechanical taxonomy below uses the direct same-fit POLE-minus-MMRd interval. A compatible or unresolved result is not equality or equivalence.

### GATA2

- TCGA_PRIMARY_CPE_N506: `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`
- TCGA_NOPURITY_N507: `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`
- CPTAC_DISCOVERY_N95: `SAME_DIRECTION_POINT_ESTIMATES_UNRESOLVED`
- CPTAC_CONFIRMATORY_N135: `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`
- CPTAC fixed-effect meta: `CPTAC_META_SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`

### SOX9

- TCGA_PRIMARY_CPE_N506: `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`
- TCGA_NOPURITY_N507: `SAME_DIRECTION_DISTINGUISHABLE_BELOW_MATERIALITY_FLOOR`
- CPTAC_DISCOVERY_N95: `SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`
- CPTAC_CONFIRMATORY_N135: `SAME_DIRECTION_POINT_ESTIMATES_UNRESOLVED`
- CPTAC fixed-effect meta: `CPTAC_META_SAME_DIRECTION_COMPATIBLE_BOTH_RESOLVED`

The equal-weight C2 is not sample-size weighted; POLE:MMRd sample proportions alone are not a class-mix artifact.

## Algebra and multiplicity

All 24 direct target/model coefficient and d reconstructions passed scaled tolerance 1e-12, including replicate-level bootstrap identity. Maximum observed replicate errors are recorded in `results/C2_RECONSTRUCTION_CHECKS.tsv`.
CPTAC fixed-effect C2 identity is not required because residual scales and contrast-specific inverse-bootstrap-variance weights differ. Direct FE(C2) discrepancies and every weight are reported in `results/CPTAC_FIXED_EFFECT_META.tsv`.
Exactly five separate descriptive BH-18 families were computed. These q values do not confer confirmatory credit.

## Interpretation boundary

Cross-cohort differences can reflect biology, composition, acquisition, platform, classifier, or scoring context; this analysis cannot identify cause. Bulk subtype-level results are not individual-patient biomarkers. Frozen verdicts are preserved and no manuscript edit was made.
