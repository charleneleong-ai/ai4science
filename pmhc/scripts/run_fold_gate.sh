#!/usr/bin/env bash
# Fold novelty over the whole cohort, including the designs the sequence gate rejected: those
# are the known-bad set the gate is validated against.
set -euo pipefail
cd "$HOME/pmhc"
export PYTHONPATH="$HOME/pmhc:$HOME/pmhc/scripts"
echo 1000 > /proc/self/oom_score_adj 2>/dev/null || true

POSES="$HOME/pmhc/cohort_poses"
rm -rf "$POSES" && mkdir -p "$POSES"
for v in bulge_core bulge_wide mut_proximal hub; do
  for f in gen_out/$v/intermediate_designs_inverse_folded/refold_cif/*.cif; do
    ln -sf "$(readlink -f "$f")" "$POSES/$(basename "$f")"
  done
done
echo "staged $(find "$POSES" -name '*.cif' | wc -l) poses"

~/.local/bin/autoresearch-wait --label pmhc-fold --ram-free-gb 8 --stable 60 --timeout 1800 || true

~/tnf/.venv-bg/bin/python scripts/fold_novelty.py "$POSES" \
  --out fold_novelty_cohort.csv --work fold_work_cohort
echo "PMHC_FOLD_DONE $(date -u +%Y-%m-%dT%H:%M:%SZ)"
