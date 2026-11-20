#!/usr/bin/env bash
# Embed train molecules with the frozen DreaMS model, train the molecule tower, then replay the holdout's reach.
# Run detached from ~/dreams on the A100; gated on room because the host is shared.
set -euo pipefail
cd "$HOME/dreams"
PY=${PY:-$HOME/casmi_seeds/.venv/bin/python}
echo 1000 > /proc/self/oom_score_adj
"$HOME/.local/bin/autoresearch-wait" --label dreams --gpu-free-gb 12 --ram-free-gb 45 --stable 120
for step in embed train reach; do $PY -u dreams_tower.py "$step"; done
