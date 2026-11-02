#!/usr/bin/env bash
# Three-way co-fold screen. See cofold_specs.py for why only two of the three conditions are
# screened on every design and the third is reserved for the finalists.
#
# Invocation gotchas inherited from the TNF campaign, both already paid for once:
#   --use_msa_server is a boolean FLAG, so passing it a value makes click read a stray
#   positional argument and exit 2. The specs set `msa: empty`, so it is simply omitted.
#   This venv lacks cuequivariance_torch, so --no_kernels is required; same numerics, slower.
set -uo pipefail
cd "$HOME/pmhc"
export PYTHONNOUSERSITE=1
echo 1000 > /proc/self/oom_score_adj 2>/dev/null || true

# One prediction on this complex (~470 residues, 1 diffusion sample) measured at ~70 s
# amortised and peaked well under 20 GB, so the gate asks for room rather than the whole card.
~/.local/bin/autoresearch-wait --label pmhc-cofold --gpu-free-gb 24 --ram-free-gb 16 \
  --stable 120 --timeout 7200 || echo "WARN gate timed out, proceeding"

echo "=== cofold start $(date -u +%FT%TZ): $(ls cofold_specs/*.yaml | wc -l) specs"
time ~/tnf/.venv-boltz/bin/boltz predict cofold_specs \
  --out_dir cofold_out \
  --override \
  --output_format mmcif \
  --diffusion_samples 1 \
  --write_full_pae \
  --no_kernels
rc=$?
echo "=== PMHC_COFOLD_DONE rc=$rc $(date -u +%FT%TZ)"
find cofold_out -name 'confidence*.json' | wc -l
