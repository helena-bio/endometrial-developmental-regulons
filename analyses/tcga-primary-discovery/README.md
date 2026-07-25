# TCGA primary discovery

This package contains the frozen TCGA-UCEC discovery analysis of the fetal Mullerian epithelial module and 20 developmental transcription-factor regulons.

Pointwise subtype association and universal single-target-deletion robustness are reported separately.

LHX1 and PAX8 retain universal deletion credit for the C1 contrast. GATA2, HOXA9, SOX9, and WT1 retain directional C2 point estimates but do not receive universal deletion credit.

The analysis does not establish causal transcription-factor activity, direct DNA binding, fetal reactivation, biomarker validity, treatment prediction, or therapeutic dependence.

## Running the code

Install the packages listed in the repository-level `requirements.txt`. Large source and intermediate files are not stored in Git. Their locations can be supplied with the `TCGA_*` environment variables defined in `src/common_v3.py`. Generated intermediates default to `work/intermediate`; frozen public outputs default to `results/`.
