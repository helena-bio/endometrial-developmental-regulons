# Cycle 6 raw analytical report

Status: COMPLETE PRODUCER EXECUTION. This is a post-hoc explanatory TCGA-UCEC sensitivity only. It does not change any frozen TCGA/CPTAC verdict, target category, manuscript, or claim.

Two post-seal runs completed with 60 model rows, 24 matched decompositions, 60 diagnostic rows, and 5,231 influence records each. All 13 scientific files are byte-identical. All 60 fitted nodes passed the frozen gates; all 24 decompositions are interpretable. Taxonomy counts: 16 largely_retained and 8 amplified.

Full-base primary CPE C2 d values: GATA2 -0.553031, SOX9 -0.527111, HOXA9 -0.646019, WT1 -0.590952, PAX8 -0.246147, LHX1 -0.277588.

GATA2/SOX9 matched results:
- primary_cpe / histology_matched / GATA2: n=505, d_base=-0.555879, d_adjusted=-0.559771, delta_d=-0.003892, percent_attenuation=-0.700%, taxonomy=amplified, status=PASS.
- primary_cpe / histology_matched / SOX9: n=505, d_base=-0.529131, d_adjusted=-0.529781, delta_d=-0.000650, percent_attenuation=-0.123%, taxonomy=amplified, status=PASS.
- primary_cpe / endometrioid_grade_matched / GATA2: n=380, d_base=-0.522793, d_adjusted=-0.433510, delta_d=+0.089283, percent_attenuation=+17.078%, taxonomy=largely_retained, status=PASS.
- primary_cpe / endometrioid_grade_matched / SOX9: n=380, d_base=-0.616573, d_adjusted=-0.560640, delta_d=+0.055932, percent_attenuation=+9.071%, taxonomy=largely_retained, status=PASS.
- no_purity / histology_matched / GATA2: n=506, d_base=-0.542232, d_adjusted=-0.546938, delta_d=-0.004705, percent_attenuation=-0.868%, taxonomy=amplified, status=PASS.
- no_purity / histology_matched / SOX9: n=506, d_base=-0.511319, d_adjusted=-0.512671, delta_d=-0.001352, percent_attenuation=-0.264%, taxonomy=amplified, status=PASS.
- no_purity / endometrioid_grade_matched / GATA2: n=381, d_base=-0.510411, d_adjusted=-0.420402, delta_d=+0.090009, percent_attenuation=+17.635%, taxonomy=largely_retained, status=PASS.
- no_purity / endometrioid_grade_matched / SOX9: n=381, d_base=-0.596497, d_adjusted=-0.537973, delta_d=+0.058524, percent_attenuation=+9.811%, taxonomy=largely_retained, status=PASS.

Diagnostics maxima: condition=4.400066, VIF=3.716361, adjusted GVIF=1.154460, hat=0.202652, single-case absolute d change=0.042131.

Upstream point-estimate reconciliation: MATCH; 54 checks, 0 mismatches, maximum absolute delta 4.441e-15.

No q values or new categories were computed. Raw Student-t p values and diagnostic permutation p values are descriptive only. No CPTAC adjusted model was run.
