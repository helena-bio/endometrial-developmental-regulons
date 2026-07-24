# Seed scheme

Master seed: `20260722`.

Subseed: `int(sha256('20260722:' + step_id).hexdigest()[:8], 16)`, with exact frozen step IDs `boot__<model>__<target>` and `perm__<model>__<target>`.
