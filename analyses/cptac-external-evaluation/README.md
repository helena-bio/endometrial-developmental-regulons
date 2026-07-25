# CPTAC external evaluation

This package contains the frozen target-level evaluation in one independent CPTAC-UCEC cohort represented by two donor-independent strata.

GATA2 and SOX9 retained the prespecified direction in both strata and met the stored target-level conditions.

Because the source TCGA effects did not receive universal single-target-deletion credit, the supported interpretation is directional concordance rather than confirmatory replication of a robustness-credited discovery.

HOXA9 and WT1 were evaluable but did not confirm. PAX8 remained sensitivity-only, and LHX1 triggered the opposite-direction veto.

## Running the code

Install the packages listed in the repository-level `requirements.txt`. Raw STAR-count files and generated matrices are not stored in Git. They default to `work/` and may be relocated with `CPTAC_WORK_DIR`, `CPTAC_INTERMEDIATE_DIR`, and `CPTAC_STAR_COUNTS_DIR`. Public linkage records are read from `data-linkage/`, and generated result tables default to `results/`. The external ESTIMATE gene set and Dou et al. supplementary workbook must be supplied separately through the documented environment variables.
