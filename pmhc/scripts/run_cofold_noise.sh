#!/usr/bin/env bash
# How much of a design's peptide-dependence drop is sampling noise?
#
# The screen runs one diffusion sample per condition, so `drop = complex_iptm - apo_iptm` is a
# difference of two single draws and carries the variance of both. Ranking 78 designs on it is
# only defensible once the spread of that difference is known, so this re-runs a spanning subset
# with five samples per condition and reports the per-design spread of the drop. Differences
# between designs smaller than that spread are not differences.
#
# Sized to the same measured envelope as run_cofold.sh: ~70 s and <20 GB per single-sample
# prediction on this ~470-residue complex.
set -uo pipefail
cd "$HOME/pmhc"
export PYTHONNOUSERSITE=1
echo 1000 > /proc/self/oom_score_adj 2>/dev/null || true

SPECS="${1:-cofold_noise_specs}"
OUT="${2:-cofold_noise_out}"
SAMPLES="${3:-5}"

n=$(ls "$SPECS"/*.yaml 2>/dev/null | wc -l)
if [ "$n" -eq 0 ]; then
  echo "no specs in $SPECS -- run: python scripts/cofold_noise.py --specs $SPECS" >&2
  exit 2
fi

~/.local/bin/autoresearch-wait --label pmhc-noise --gpu-free-gb 24 --ram-free-gb 16 \
  --stable 120 --timeout 7200 || echo "WARN gate timed out, proceeding"

echo "=== noise start $(date -u +%FT%TZ): $n specs x $SAMPLES samples"
time ~/tnf/.venv-boltz/bin/boltz predict "$SPECS" \
  --out_dir "$OUT" \
  --override \
  --output_format mmcif \
  --diffusion_samples "$SAMPLES" \
  --write_full_pae \
  --no_kernels
rc=$?
echo "=== PMHC_NOISE_DONE rc=$rc $(date -u +%FT%TZ)"
find "$OUT" -name 'confidence*.json' | wc -l
