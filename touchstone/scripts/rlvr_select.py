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
    # --deep  → also run the MLIP relax+MD tiers, so the reward selects for dynamically-stable
    #           sites, not just geometrically-clean ones (needs a GPU; only spends MLIP on
    #           geometry-plausible designs). Run in the mace env on the GPU box.
    # --selectivity Ni2+,Cu2+,Co2+  → fold the metal-swap ΔE tier into the reward, selecting designs
    #           whose *target* metal binds most favourably — the lever for *selective* binders (e.g.
    #           Cu2+ over Ni2+/Co2+ in a mixed leachate). Needs --deep.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import typer

from touchstone.boltzgen import boltzgen_confidence
from touchstone.reward import rank_structures
from touchstone.service import mlip_backbone


def main(
    npz_dir: Path = typer.Option(..., help="BoltzGen fold_out_npz dir (confidence per design)"),
    cif_dir: Path = typer.Option(..., help="matching refold_cif dir (structures to score)"),
    out: Path = typer.Option(..., help="output dir for the selected fine-tuning set + rewards"),
    metal: str = typer.Option("Ni2+", help="target metal"),
    keep: str = typer.Option("trust", help="'trust' (only TRUST designs) or an int N (top-N by reward)"),
    deep: bool = typer.Option(False, "--deep", help="fold the MLIP relax+MD tiers into the reward (needs a GPU + touchstone[mace])"),
    selectivity: str = typer.Option("", "--selectivity", help="comma-separated competing metals for the selectivity tiers, e.g. Ni2+,Cu2+,Co2+ — motif enrichment (MetalPDB occupancy) runs on CPU; the MLIP metal-swap ΔE tier is added by --deep (and currently defers: no backbone passes the Irving-Williams gate)"),
) -> None:
    # validate before touching the filesystem — a bad flag combination shouldn't leave an empty out/
    calc = mlip_backbone() if deep else None  # build the MACE backbone once, share across the batch
    if deep and calc is None:  # honor an explicit --deep: never silently ship a geometry-only reward as "deep"
        raise typer.Exit("--deep needs a GPU + touchstone[mace]; refusing to score geometry-only")

    out.mkdir(parents=True, exist_ok=True)
    dataset = out / "dataset"
    dataset.mkdir(exist_ok=True)

    # rank_structures does the batch verify+reward+sort; gate_defer spends MLIP only on geometry-plausible designs
    ranked = rank_structures(sorted(cif_dir.glob("*.cif")), metal, deep=deep, gate_defer=deep, calc=calc,
                             selectivity_metals=selectivity.split(",") if selectivity else None)
    scored = []
    for r in ranked:
        cif = Path(r["structure"])
        v = r.get("verifiers", {})
        scored.append({"design": cif.stem, "cif": str(cif), "reward": r["reward"],
                       "consensus": r["consensus"],
                       "boltzgen_iptm": (boltzgen_confidence(npz_dir / f"{cif.stem}.npz") or {}).get("iptm"),
                       "mlip": v.get("mlip", {}).get("label"), "mlip_md": v.get("mlip_md", {}).get("label"),
                       "selectivity": v.get("selectivity", {}).get("label")})

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
