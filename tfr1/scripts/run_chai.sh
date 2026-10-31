#!/usr/bin/env bash
#   run_chai.sh [in_dir] [out_dir] [scores.jsonl] [after_pattern]
#
# Directories are relative to $CHAI_WORK, which is bound to /work in the container. The driver is
# chai_independent.py, kept beside the FASTAs in $CHAI_WORK rather than vendored into this repo:
# the host's repo clone is stale and does not carry it, so a bind would pin an older revision than
# the one that ran. Both it and this launcher belong in autoresearch beside autoresearch-wait and
# autoresearch-heartbeat once the apptainer invocation stops being host-specific.
#
# `after_pattern` holds the job behind a predecessor still on the GPU, which the memory gate below
# does not cover on its own: Chai between designs can momentarily leave room for a second pass
# that then cannot fit. It is a pgrep pattern and must be bracketed -- '[c]hai_independent', not
# 'chai_independent' -- or it matches this script's own command line and waits forever.
# autoresearch-wait growing an --after-pattern of its own would retire both the loop and the
# bracket trick.
set -u
IN_DIR=${1:-chai_in}
OUT_DIR=${2:-chai_out}
SCORES=${3:-chai_scores.jsonl}
AFTER=${4:-}

WORK=${CHAI_WORK:-$HOME/tfr1}
SIF=${CHAI_SIF:-$HOME/RFdiffusion2/sifs/chai.sif}
DOWNLOADS=${CHAI_DOWNLOADS:-$HOME/tnf/chai_downloads}
AFTER_TIMEOUT=${CHAI_AFTER_TIMEOUT:-21600}

export PATH=$HOME/.local/bin:$PATH

# Checks in while it holds, so a stuck predecessor is distinguishable from a long one.
waited=0
while [ -n "$AFTER" ] && pgrep -f "$AFTER" >/dev/null 2>&1; do
  if [ "$waited" -ge "$AFTER_TIMEOUT" ]; then
    echo "AFTER_TIMEOUT $AFTER still running after ${waited}s $(date -u +%FT%TZ)"
    exit 1
  fi
  [ $((waited % 900)) -eq 0 ] && echo "HOLDING for $AFTER ${waited}s $(date -u +%FT%TZ)"
  sleep 60
  waited=$((waited + 60))
done
[ -n "$AFTER" ] && echo "PREDECESSOR_CLEAR after ${waited}s $(date -u +%FT%TZ)"

autoresearch-wait --label tfr1-chai --gpu-free-gb 25 --ram-free-gb 16 --stable 120 --timeout 10800
gate=$?
if [ "$gate" -eq 127 ]; then
  echo "GATE_MISSING autoresearch-wait is not on PATH"
  exit 127
elif [ "$gate" -ne 0 ]; then
  echo "GATE_TIMEOUT no room after 10800s"
  exit "$gate"
fi

apptainer exec --nv \
  -B "$WORK":/work \
  -B "$DOWNLOADS":/usr/envs/chai/lib/python3.11/site-packages/downloads \
  "$SIF" python /work/chai_independent.py "/work/$IN_DIR" "/work/$OUT_DIR" "/work/$SCORES"
rc=$?
echo "CHAI_TFR1_DONE $IN_DIR rc=$rc $(date -u +%FT%TZ)"
exit $rc
