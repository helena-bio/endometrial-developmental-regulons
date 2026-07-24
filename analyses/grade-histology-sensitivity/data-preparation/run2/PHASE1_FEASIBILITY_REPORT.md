# TASK B Phase-1 clinical feasibility report

Mechanical state: RESTRICT_TCGA_ONLY. This is an outcome-blind feasibility state, not a scientific verdict.

TCGA-UCEC has exact frozen patient and analytic-sample linkage. Grade comes from cBioPortal PanCanAtlas sample-level GRADE tokens; because the locked source package does not document one FIGO construct across all histologies, grade is permitted only inside the endometrioid subset. Histology comes from patient-level ICD-O-3 morphology and supports the binary endometrioid versus non-endometrioid node. The three-level histology option fails the frozen support rule for other non-endometrioid cases. Supported later nodes are base, base plus binary histology, and endometrioid-only with grade. Optional grade stratification is forbidden because the low-grade endometrioid stratum fails the frozen leverage rule, even though the high-grade stratum passes. All-histology grade and grade-plus-histology are forbidden.

CPTAC Discovery uses cBioPortal HISTOLOGIC_GRADE_FIGO and HISTOLOGIC_TYPE plus pinned PDC000125 case-level morphology/tumor_grade and the pinned biospecimen roster. One analytic case has two distinct Primary Tumor samples with no frozen exact selector. Discovery also has cBio/PDC clinical disagreements; conflicts remain missing. Grade high-category and subtype support fail, and reconciled histology lacks a supported varying clinical factor. The stratum is BLOCKED.

CPTAC Confirmatory uses pinned PDC000439 case-level morphology, primary diagnosis, tumor_grade, and the pinned biospecimen roster. Two analytic cases have one primary sample but two aliquots with no frozen exact selector. The generic tumor_grade field is not locally documented as FIGO, and retained histology is single-level. The stratum is BLOCKED.

Definitions are compatible-but-not-identical for TCGA versus CPTAC histology and for the explicitly FIGO Discovery grade field; the Confirmatory generic grade is not harmonizable to the frozen FIGO mapping without inference. No clinical parameter may be pooled across cohorts. Formal Phase 2 is not currently permitted. After fresh reviewer verification and separate authorization, only the supported TCGA nodes may proceed; CPTAC adjusted models remain forbidden.

The exact full molecular design was not tested in Phase 1. Rank, condition, VIF/GVIF, leverage, and influence diagnostics must be repeated on the exact Phase-2 design. No molecular conclusion and no attenuation statement was made.
