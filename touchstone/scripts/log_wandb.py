"""Verify a pool of designed metal binders and log the filtering to Weights & Biases.

Demonstrates the verifier's value prop: generate many candidates (RFdiffusionAA →
LigandMPNN packs), let the geometry verifier filter to the trustworthy few. Logs a
per-design table, summary counts, three interactive charts, and a 3D view of the best
design — all W&B-native (no matplotlib).

    uv run --extra viz python scripts/log_wandb.py <design_pdb_dir>

Reads WANDB_API_KEY / WANDB_PROJECT from touchstone/.env (gitignored).
"""

from __future__ import annotations

import glob
import math
from collections import Counter
from pathlib import Path

import typer
import wandb
from dotenv import load_dotenv

from touchstone import BinderDesign, GeometryVerifier, PDBReference, coordination_site_from_pdb, rank

load_dotenv(Path(__file__).parent.parent / ".env")


def verdict_of(v) -> str:
    return "defer" if v.ood else ("trust" if v.trust else "weak")


def load(design_dir: str) -> list[tuple[BinderDesign, str]]:
    out = []
    for pdb in sorted(glob.glob(f"{design_dir}/**/*.pdb", recursive=True)):
        try:
            site = coordination_site_from_pdb(pdb, "NI", "Ni2+")
        except ValueError:
            continue
        out.append((BinderDesign(Path(pdb).stem, site, "rfaa+ligmpnn", float("nan")), pdb))
    return out


def main(design_dir: str = "ligmpnn_out") -> None:
    ref = PDBReference()
    g = ref.geometry("Ni2+")
    verifier = GeometryVerifier(ref)

    pool = load(design_dir)
    paths = {d.sequence: p for d, p in pool}
    ranked = rank([d for d, _ in pool], verifier)

    run = wandb.init(
        project="touchstone", name="ligmpnn-nickel-filter",
        config={"n_designs": len(ranked), "metal": "Ni2+",
                "ref_mean": g.bond_length_mean, "ref_std": g.bond_length_std, "ref_cn_range": list(g.cn_range)},
    )

    cols = ["design", "CN", "score", "log_score", "verdict", "donors", "min_bond", "max_bond"]
    designs = wandb.Table(columns=cols)
    bonds = wandb.Table(columns=["bond_length"])
    counts: Counter = Counter()
    for d, v in ranked:
        bl = d.site.bond_lengths()
        counts[verdict_of(v)] += 1
        designs.add_data(d.sequence, d.site.coordination_number, v.score,
                         math.log10(max(v.score, 1e-30)), verdict_of(v), "".join(d.site.ligand_elems),
                         round(float(bl.min()), 2) if len(bl) else 0.0,
                         round(float(bl.max()), 2) if len(bl) else 0.0)
        for b in bl:
            bonds.add_data(round(float(b), 2))

    verdict_counts = wandb.Table(data=[[k, counts[k]] for k in ("trust", "weak", "defer")],
                                 columns=["verdict", "count"])
    run.log({
        "designs": designs,
        "filtering": wandb.plot.bar(verdict_counts, "verdict", "count",
                                    title="Verifier filtering of the design pool"),
        "score_vs_CN": wandb.plot.scatter(designs, "CN", "log_score",
                                          title="Score landscape (log) vs coordination number"),
        "design_bond_lengths": wandb.plot.histogram(bonds, "bond_length",
                                                    title=f"Designed donor bonds vs real Ni²⁺ {g.bond_length_mean}±{g.bond_length_std} Å"),
    })

    for k in ("trust", "weak", "defer"):
        run.summary[f"n_{k}"] = counts[k]
    best, bv = ranked[0]
    run.summary["best_score"] = bv.score
    run.summary["best_design"] = f"{best.sequence} ({verdict_of(bv)}, score {bv.score:.2g})"

    # 3D structures: the best design overall + the best of each verdict class, so you can
    # rotate them and compare the coordination (wandb.Molecule(open(pdb)) per the W&B report).
    mols, seen = {"best_design_3d": wandb.Molecule(open(paths[best.sequence]))}, set()
    for d, v in ranked:
        cls = verdict_of(v)
        if cls not in seen:
            seen.add(cls)
            mols[f"3d_best_{cls}"] = wandb.Molecule(open(paths[d.sequence]))
    try:
        run.log(mols)
    except Exception as e:  # 3D view is a nice-to-have, not worth failing the run
        print(f"(skipped 3D molecules: {e})")

    print(f"logged {len(ranked)} designs → {run.url}")
    print(f"  trust={counts['trust']}  weak={counts['weak']}  defer={counts['defer']}  best={bv.score:.2g}")
    run.finish()


if __name__ == "__main__":
    typer.run(main)
