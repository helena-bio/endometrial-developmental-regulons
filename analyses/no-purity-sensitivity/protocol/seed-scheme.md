# TASK-030 cycle-2 frozen stochastic step namespace

Recorded before any cycle-2 leave-one-gene-out stochastic calculation.

- Master seed: `20260713`.
- Hash rule: `int(sha256("20260713:" + step_id).hexdigest()[:8], 16)`.
- Bootstrap step ID: `loo_boot__{config}__M3_{target}__drop_{dropped_gene}`.
- Permutation step ID: `loo_perm__{config}__M3_{target}__drop_{dropped_gene}`.
- `config` is exactly `PRIMARY` or `SENS_nopurity`.
- `target` and `dropped_gene` are the case-sensitive symbols in the sealed edge ledger.
- One unique mapped target gene is deleted at a time; every admitted edge row for that symbol is deleted.

This namespace will not be revised based on output.
