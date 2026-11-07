#!/usr/bin/env bash
# Launch the seed runs detached (PPID=1), each with a heartbeat. Prep prints nothing until it finishes,
# so the heartbeat's state line reports peaks.bin growth and the stall window covers the prep.
set -euo pipefail
export PATH=$HOME/.local/bin:$PATH
cd ~/casmi_seeds
echo "b18b53c84c3c6a5cc2604322e54b19f0  data/train.parquet" | md5sum -c -
TS=$(date -u +%Y%m%dT%H%M%SZ)
for s in "$@"; do
  log=logs/seed_${s}_${TS}.log
  setsid nohup bash -c "echo 1000 > /proc/self/oom_score_adj; exec .venv/bin/python -u run_seed.py $s" \
    </dev/null >>"$log" 2>&1 & disown
  setsid nohup autoresearch-heartbeat --label "seed$s" --pattern "[r]un_seed.py $s" --log "$log" \
    --interval 1200 --stall-after 7200 \
    --state-cmd "echo peaks=\$(du -sh runs/$s/work/peaks.bin 2>/dev/null | cut -f1) epochs=\$(grep -c -E '^(HOLDOUT|REFIT)' $log)" \
    </dev/null >>logs/heartbeat_seed${s}_${TS}.log 2>&1 & disown
done
sleep 5
ps -eo pid,ppid,etime,cmd | grep -E '[r]un_seed.py|[a]utoresearch-heartbeat'
