"""Reject-sampling selection step for BoltzGen RLVR — verifier-as-reward, rung 1.

Score a BoltzGen-generated pool with touchstone (the CSD-calibrated reward), keep the
high-reward / TRUST designs, and write them out as the fine-tuning set. BoltzGen has no DPO,
only a standard supervised trainer (`boltzgen.task.train`, checkpoint-resume), so the loop is
RAFT (reward-ranked fine-tuning): train on the verifier's winners, repeat.

    BoltzGen generate  →  rlvr_select (this)  →  BoltzGen resume-train on the kept set  →  repeat

CSD metal knowledge enters through touchstone's reward (its prior is the CSD pull), not as raw
training data — which is why this works for a *protein* generator. See docs/specs for the loop.

    uv run python scripts/rlvr_select.py \
        --npz-dir boltzgen_out/.../fold_out_npz --cif-dir boltzgen_out/.../refold_cif \
        --out rlvr_round1 --metal Ni2+ --keep trust      # strict: only TRUST designs
    # or --keep 8  → top-8 by reward (softer; use to bootstrap when the pool has few/no TRUST)
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import typer

from touchstone.reward import reward_from_result
from touchstone.service import verify_structure


def main(
    npz_dir: Path = typer.Option(..., help="BoltzGen fold_out_npz dir (confidence per design)"),
    cif_dir: Path = typer.Option(..., help="matching refold_cif dir (structures to score)"),
    out: Path = typer.Option(..., help="output dir for the selected fine-tuning set + rewards"),
    metal: str = typer.Option("Ni2+", help="target metal"),
    keep: str = typer.Option("trust", help="'trust' (only TRUST designs) or an int N (top-N by reward)"),
) -> None:
    out.mkdir(parents=True, exist_ok=True)
    dataset = out / "dataset"
    dataset.mkdir(exist_ok=True)

    scored = []
    for cif in sorted(cif_dir.glob("*.cif")):
        try:
            r = verify_structure(cif, metal)
        except Exception as e:
            scored.append({"design": cif.stem, "reward": 0.0, "consensus": "error", "error": str(e)})
            continue
        npz_path = npz_dir / f"{cif.stem}.npz"
        iptm = None
        if npz_path.exists():
            npz = np.load(npz_path, allow_pickle=True)
            for k in ("design_to_target_iptm", "ligand_iptm", "iptm"):
                if k in npz:
                    iptm = float(np.atleast_1d(npz[k]).astype(float).max())
                    break
        scored.append({"design": cif.stem, "cif": str(cif), "reward": reward_from_result(r),
                       "consensus": r["consensus"], "boltzgen_iptm": iptm})
    scored.sort(key=lambda s: s["reward"], reverse=True)

    if keep == "trust":
        winners = [s for s in scored if s["consensus"] == "trust"]
    else:
        winners = scored[: int(keep)]

    for s in winners:
        if s.get("cif"):
            shutil.copy(s["cif"], dataset / Path(s["cif"]).name)
    (out / "rewards.jsonl").write_text("".join(json.dumps(s) + "\n" for s in scored))

    n_trust = sum(s["consensus"] == "trust" for s in scored)
    print(f"scored {len(scored)} designs · {n_trust} TRUST · kept {len(winners)} → {dataset}")
    print(f"rewards → {out / 'rewards.jsonl'}")
    if not winners:
        print("⚠ no designs selected — lower the bar (--keep N) or improve generation before fine-tuning")


if __name__ == "__main__":
    typer.run(main)
