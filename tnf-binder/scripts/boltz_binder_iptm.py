"""Extract the BINDER-to-target ipTM from Boltz confidence json, not the global ipTM.

Global ipTM averages every inter-chain pair, so a multi-chain target inflates it with its own
internal interfaces. For a binder against the TNF-alpha trimer the three TNF-TNF pairs score
~0.88 and swamp the one number that matters. This is the same artifact that made the BoltzGen
dimer->trimer switch look like a win (global 0.20 -> 0.52) when the binder interface was flat.

Chain 0 is the designed binder (listed first in the spec); chains 1..n are target protomers.
"""
import json
import statistics as st
import sys
from pathlib import Path

for d in [Path(x) for x in sys.argv[1:]]:
    files = sorted(d.rglob("confidence*.json")) if d.is_dir() else [d]
    if not files:
        print(f"no confidence json under {d}")
        continue
    print(f"\n===== {d.name}")
    per_model = []
    for f in files:
        c = json.loads(f.read_text())
        pc = c.get("pair_chains_iptm", {})
        if not pc:
            continue
        ids = sorted(pc.keys(), key=lambda x: int(x))
        binder, targets = ids[0], ids[1:]
        bt = [float(pc[binder][t]) for t in targets] + \
             [float(pc[t][binder]) for t in targets]
        tt = [float(pc[a][b]) for a in targets for b in targets if a != b]
        per_model.append((f.name, float(c["iptm"]), max(bt), st.mean(bt),
                          st.mean(tt) if tt else float("nan"),
                          float(c.get("complex_ipde", float("nan"))),
                          float(pc[binder][binder])))
    if not per_model:
        continue
    print(f"  {'model':34s} {'global':>7s} {'bind2tgt_best':>14s} {'bind2tgt_mean':>14s} "
          f"{'tgt2tgt_mean':>13s} {'i_pde':>7s} {'binder_ptm':>11s}")
    for n, g, bmax, bmean, tmean, ipde, bptm in per_model:
        print(f"  {n:34s} {g:7.3f} {bmax:14.3f} {bmean:14.3f} {tmean:13.3f} "
              f"{ipde:7.2f} {bptm:11.3f}")
    print(f"\n  VERDICT: binder-to-target ipTM best={max(p[2] for p in per_model):.3f} "
          f"mean={st.mean([p[3] for p in per_model]):.3f}")
    print(f"           global ipTM mean={st.mean([p[1] for p in per_model]):.3f} "
          f"-- inflated by target-internal interfaces at {st.mean([p[4] for p in per_model]):.3f}")
