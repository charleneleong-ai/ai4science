#!/usr/bin/env bash
# Step 0 of reports/2026-11-14-scorer-options.md: build the structure-disjoint holdout, train the calibration
# variants on it, and judge the calibration gate. Run detached; exits non-zero when the gate fails.
#   kit     the shipped schedule, 16 x 500 steps (the reference)
#   long    40 x 2000 steps (public LB: lower than kit, 0.453 against 0.472)
#   s1..s3  the kit schedule with three new seeds; kit+s1+s2+s3 is the seed blend (public LB: -0.016)
# Gate (fixed before any result): kit beats long, and kit beats kit+s1+s2+s3 (eval_holdout.py gate).
set -euo pipefail
cd "${CASMI_HOLDOUT_ROOT:-$HOME/casmi_holdout}"
PY=${PY:-$HOME/casmi_seeds/.venv/bin/python}
echo 1000 > /proc/self/oom_score_adj
$PY -u build_holdout.py
$PY -u eval_holdout.py lists > logs/lists.log 2>&1 &
pids=($!)
specs=("kit 20261023 16 500" "long 20261023 40 2000" "s1 20261101 16 500" "s2 20261102 16 500" "s3 20261103 16 500")
names=()
for v in "${specs[@]}"; do
  set -- $v
  names+=("$1")
  $PY -u train_variant.py "$@" > "logs/train_$1.log" 2>&1 &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
$PY -u eval_holdout.py score "${names[@]}"
$PY -u eval_holdout.py compare kit long kit+s1+s2+s3 s1 s2 s3
$PY -u eval_holdout.py gate
