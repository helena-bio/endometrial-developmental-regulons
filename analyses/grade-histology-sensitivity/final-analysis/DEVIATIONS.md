# Cycle 6 deviations

## D1 - DISCLOSED_PRESEAL_PLUMBING

Two pre-seal synthetic computations passed, but strace followed live tee helpers and deadlocked. Both attempts and exit-137 tracer terminations are preserved. The audit wrapper was corrected before sealing; the third traced synthetic test passed and exited cleanly.

Scientific effect: none; no outcome was opened and no fit used real outcomes.

## D2 - DISCLOSED_AUDIT_BOOTSTRAP

One post-patch attempt to invoke audit_exec.sh directly failed with Permission denied before the wrapper could log itself; the exact command and error are disclosed in COMMAND_COVERAGE_NOTE.md and the tool transcript. The immediately following bash-invoked command restored executable mode and is logged.

Scientific effect: none; pre-seal and before outcome access.

## D3 - DISCLOSED_RECONCILIATION_SCOPE

While resolving Task-028 module identifiers after output freeze, a filtered command printed GATA2/SOX9 C1 and C3 rows in addition to the authorized C2 rows. This occurred only after both Cycle-6 outputs were read-only and byte-frozen; no analysis code, taxonomy, output, or verdict changed.

Scientific effect: none on Cycle-6 computation; broader-than-needed upstream read disclosed.

## D4 - DISCLOSED_UPSTREAM_COVERAGE

The pinned Task-030 point table contains GATA2, SOX9, HOXA9, and WT1 only. The pinned Task-028 primary and no-purity tables supplied all six targets, including PAX8 and LHX1.

Scientific effect: none; all six point estimates reconciled to Task-028.

## D5 - DISCLOSED_PROVENANCE_CLOSEOUT

The first provenance-only finalizer attempt stopped before writing the manifest because the trace auditor treated bwrap ancestor-directory traversal as an unexpected workspace access. The raw traceback and exit status are preserved. The auditor was corrected to allow only lexical ancestors of the sealed Cycle-6 root and then rerun; frozen run1/run2 files were never changed.

Scientific effect: none; both scientific executions had already succeeded and been made read-only.
