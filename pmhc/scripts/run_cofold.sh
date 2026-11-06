#!/usr/bin/env bash
# Boltz co-fold runner. One script for every pass, because they differ only in which specs go
# in, where the output lands, and how many diffusion samples each prediction draws:
#
#   pilot:     run_cofold.sh                                     # every spec, 1 sample
#   probe:     run_cofold.sh cofold_noise_specs cofold_noise_out 5
#   cohort:    run_cofold.sh cofold_specs       cofold_full5_out 5
#   wild type: run_cofold.sh cofold_wt_specs    cofold_wt_out    5 2
#   replicate: run_cofold.sh cofold_control_specs cofold_ctrl_out 5 8
#
# A fourth argument runs the pass that many times as SEPARATE invocations, into OUT/rep1..repN.
# That is the only way to measure a between-run component: Boltz seeds once per process, so
# more samples inside one invocation share that state. Written as a replicate count on the one
# runner rather than a second script, because a fork of this file immediately lost the venv
# isolation, the spec-count precheck, the gate's timeout fallback and the evidence count below.
# `cofold_followup.py select-control` writes the control-only spec directory it takes.
#
# After two or more passes over the same specs, `cofold_followup.py pool-replicates` turns them
# into the one table the submission ranks on; a single pass is read by `report-spread` instead.
#
# The cohort pass is the measurement of record: the probe over a 26-spec subset showed that a
# single diffusion sample reverses the sign of a design's peptide-dependence drop, so the whole
# cohort is re-run at five and the one-sample pass is retrospectively a pilot.
#
# See cofold_specs.py for why only two of the three conditions are screened on every design,
# and cofold_followup.py for what reads each pass.
#
# Invocation gotchas inherited from the TNF campaign, both already paid for once:
#   --use_msa_server is a boolean FLAG, so passing it a value makes click read a stray
#   positional argument and exit 2. The specs set `msa: empty`, so it is simply omitted.
#   This venv lacks cuequivariance_torch, so --no_kernels is required; same numerics, slower.
set -uo pipefail
cd "$HOME/pmhc"
export PYTHONNOUSERSITE=1
echo 1000 > /proc/self/oom_score_adj 2>/dev/null || true

SPECS="${1:-cofold_specs}"
OUT="${2:-cofold_out}"
SAMPLES="${3:-1}"
REPLICATES="${4:-1}"
LABEL="$(basename "$SPECS")"

n=$(find "$SPECS" -name '*.yaml' 2>/dev/null | wc -l)
if [ "$n" -eq 0 ]; then
  echo "no specs in $SPECS" >&2
  exit 2
fi

# One prediction on this complex (~470 residues, 1 diffusion sample) measured at ~21 s and
# peaked well under 20 GB, so the gate asks for room rather than the whole card.
~/.local/bin/autoresearch-wait --label "$LABEL" --gpu-free-gb 24 --ram-free-gb 16 \
  --stable 120 --timeout 7200 || echo "WARN gate timed out, proceeding"

echo "=== $LABEL start $(date -u +%FT%TZ): $n specs x $SAMPLES sample(s) x $REPLICATES run(s)"
for i in $(seq 1 "$REPLICATES"); do
  # One out_dir per invocation only when there is more than one, so a single pass keeps the
  # layout every existing reader was written against.
  if [ "$REPLICATES" -gt 1 ]; then
    dest="$OUT/rep$i"
    # Re-gate before each one: the card is shared, and the previous invocation has just freed it.
    ~/.local/bin/autoresearch-wait --label "$LABEL-rep$i" --gpu-free-gb 24 --ram-free-gb 16 \
      --stable 120 --timeout 7200 || echo "WARN gate timed out, proceeding"
  else
    dest="$OUT"
  fi
  time ~/tnf/.venv-boltz/bin/boltz predict "$SPECS" \
    --out_dir "$dest" \
    --override \
    --output_format mmcif \
    --diffusion_samples "$SAMPLES" \
    --write_full_pae \
    --no_kernels
  rc=$?
  echo "=== PMHC_COFOLD_REP_DONE $LABEL $i/$REPLICATES rc=$rc $(date -u +%FT%TZ)"
done
echo "=== PMHC_COFOLD_DONE $LABEL rc=$rc $(date -u +%FT%TZ)"
find "$OUT" -name 'confidence*.json' | wc -l
