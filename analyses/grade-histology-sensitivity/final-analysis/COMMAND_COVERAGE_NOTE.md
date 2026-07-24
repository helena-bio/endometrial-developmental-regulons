# Command coverage note

The first shell command created the Cycle-6 directory and is recorded verbatim in
audit_bootstrap.txt. The next action created audit_exec.sh and the bootstrap note
through apply_patch; apply_patch is a tool action rather than a shell command.

Every subsequently executed shell command was routed through audit_exec.sh and has
per-command metadata, stdout, stderr, exit status, timestamps, and hashes under
audit/commands/, except one invocation that could not start the wrapper:

    experiments/taskB_grade_histology/phase2_execution_cycle6/audit_exec.sh -- chmod 0555 experiments/taskB_grade_histology/phase2_execution_cycle6/audit_exec.sh

It failed immediately with:

    /bin/bash: line 1: experiments/taskB_grade_histology/phase2_execution_cycle6/audit_exec.sh: Permission denied

This happened after apply_patch reset executable mode and before any outcome access.
The immediately following bash-invoked wrapper command restored mode and is logged.
Pre-seal tracer deadlocks and forced tracer terminations are retained as exit-137
commands. COMMANDS_CHRONOLOGICAL.tsv is sorted by command start ID; raw immutable
command evidence remains under audit/.
