"""Optional Weights & Biases wiring for touchstone experiments.

The verifier core stays dependency-light (numpy + typer); W&B and dotenv are imported
lazily here so `import touchstone` never pulls them. The key is read from
touchstone/.env (gitignored). Typical use from a script:

    from touchstone import tracking

    run = tracking.init("ligmpnn-nickel-filter", config={"n": len(ranked)})
    tracking.log_ranked(run, ranked, reference)          # table + bar/scatter/histogram
    tracking.log_molecules(run, {"best": best_pdb})      # 3D wandb.Molecule views
    run.finish()
"""

from __future__ import annotations

import math
import os
from collections import Counter
from pathlib import Path

_ENV = Path(__file__).resolve().parents[2] / ".env"  # touchstone/.env


def init(name: str, config: dict | None = None, project: str | None = None):
    """Start a W&B run, loading WANDB_API_KEY / WANDB_PROJECT from touchstone/.env."""
    import wandb
    from dotenv import load_dotenv

    load_dotenv(_ENV)
    return wandb.init(
        project=project or os.environ.get("WANDB_PROJECT", "touchstone"),
        name=name,
        config=config or {},
    )


def log_ranked(run, ranked, reference, metal: str = "Ni2+") -> Counter:
    """Log a ranked design pool: per-design table + filtering bar + score/CN scatter +
    a donor bond-length histogram. Returns the trust/weak/defer counts."""
    import wandb

    g = reference.geometry(metal)
    designs = wandb.Table(
        columns=["design", "CN", "score", "log_score", "verdict", "donors", "min_bond", "max_bond"]
    )
    bonds = wandb.Table(columns=["bond_length"])
    counts: Counter = Counter()
    for d, v in ranked:
        counts[v.label] += 1
        bl = d.site.bond_lengths()
        designs.add_data(
            d.sequence, d.site.coordination_number, v.score, math.log10(max(v.score, 1e-30)),
            v.label, "".join(d.site.ligand_elems),
            round(float(bl.min()), 2) if len(bl) else 0.0,
            round(float(bl.max()), 2) if len(bl) else 0.0,
        )
        for b in bl:
            bonds.add_data(round(float(b), 2))

    verdicts = wandb.Table(data=[[k, counts[k]] for k in ("trust", "weak", "defer")],
                           columns=["verdict", "count"])
    run.log({
        "designs": designs,
        "filtering": wandb.plot.bar(verdicts, "verdict", "count", title="Verifier filtering of the pool"),
        "score_vs_CN": wandb.plot.scatter(designs, "CN", "log_score", title="Score (log) vs coordination number"),
        "design_bond_lengths": wandb.plot.histogram(
            bonds, "bond_length", title=f"Donor bonds vs real {metal} {g.bond_length_mean}±{g.bond_length_std} Å"),
    })
    for k in ("trust", "weak", "defer"):
        run.summary[f"n_{k}"] = counts[k]
    return counts


def log_molecules(run, named_pdbs: dict[str, str]) -> None:
    """Log each PDB as an interactive 3D wandb.Molecule under its given key."""
    import wandb

    run.log({key: wandb.Molecule(open(pdb)) for key, pdb in named_pdbs.items()})
