#!/usr/bin/env bash
set -euo pipefail

PKG_DIR=data/external/original-workspace/revgate-tcga-no-purity-verify/experiments/task030_verify_current_main
VENV_PY=data/external/original-workspace/task028-freeze-b-draft/verifier_v3/venv_match/bin/python

cd "$PKG_DIR/gene_loo"
"$VENV_PY" -m py_compile run_task030_cycle2_gene_loo.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$VENV_PY" run_task030_cycle2_gene_loo.py > run.stdout.log 2> run.stderr.log

cd "$PKG_DIR/base_audit"
"$VENV_PY" -m py_compile run_task030_audit.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$VENV_PY" run_task030_audit.py > run.stdout.log 2> run.stderr.log

cd "$PKG_DIR"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$VENV_PY" build_delivery.py --capture-first

cd "$PKG_DIR/gene_loo"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$VENV_PY" run_task030_cycle2_gene_loo.py > ../logs/gene_loo_second.stdout.log 2> ../logs/gene_loo_second.stderr.log

cd "$PKG_DIR/base_audit"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$VENV_PY" run_task030_audit.py > ../logs/base_second.stdout.log 2> ../logs/base_second.stderr.log

cd "$PKG_DIR"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$VENV_PY" build_delivery.py --final
sha256sum -c SHA256SUMS.txt
