#!/usr/bin/env bash
# BoltzGen against the intact KRAS G12D pMHC, four hotspot variants over the exposed peptide
# bulge. 64 designs each; the smoke test ran at about 130 s per design on a free A100, so the
# whole sweep is roughly nine hours.
#
# The venv bin goes on PATH first because boltzgen spawns a subprocess calling bare `python`,
# which does not exist on a detached PATH -- the same failure that killed the first GLP1R run.
#
# Expect a third of the output to be beta-2 microglobulin regurgitated from the target; the
# novelty gate in scripts/novelty.py rejects it, and nothing downstream should read a score
# before that gate has run, because the plagiarised poses outscore every genuine design.
set -u
cd ~/tnf
export PATH=$HOME/tnf/.venv-bg/bin:$HOME/.local/bin:$PATH PYTHONNOUSERSITE=1

NUM=${NUM_DESIGNS:-64}
OUT=${OUT_DIR:-$HOME/pmhc/gen_out}

autoresearch-wait --label pmhc-gen --gpu-free-gb 25 --ram-free-gb 16 --stable 120 --timeout 10800 \
  || { echo "GATE_TIMEOUT no room"; exit 1; }
echo 1000 > /proc/self/oom_score_adj 2>/dev/null || true

for S in bulge_core bulge_wide mut_proximal hub; do
  echo "=== $S $(date -u +%FT%TZ) ==="
  ./.venv-bg/bin/boltzgen run "$HOME/pmhc/specs/$S.yaml" \
    --protocol protein-anything --num_designs "$NUM" \
    --output "$OUT/$S" --devices 1 || echo "$S FAILED rc=$?"
done
echo "PMHC_GEN_DONE $(date -u +%FT%TZ)"
