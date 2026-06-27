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


def _verdict_bar(counts: Counter, title: str):
    """A trust/weak/defer count bar — the verdict-distribution plot shared by the loggers."""
    import wandb

    table = wandb.Table(data=[[k, counts[k]] for k in ("trust", "weak", "defer")], columns=["verdict", "count"])
    return wandb.plot.bar(table, "verdict", "count", title=title)


def _set_verdict_summary(run, counts: Counter) -> None:
    for k in ("trust", "weak", "defer"):
        run.summary[f"n_{k}"] = counts[k]


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

    run.log({
        "designs": designs,
        "filtering": _verdict_bar(counts, "Verifier filtering of the pool"),
        "score_vs_CN": wandb.plot.scatter(designs, "CN", "log_score", title="Score (log) vs coordination number"),
        "design_bond_lengths": wandb.plot.histogram(
            bonds, "bond_length", title=f"Donor bonds vs real {metal} {g.bond_length_mean}±{g.bond_length_std} Å"),
    })
    _set_verdict_summary(run, counts)
    return counts


def log_molecules(run, named_pdbs: dict[str, str]) -> None:
    """Log each PDB as an interactive 3D wandb.Molecule under its given key."""
    import wandb

    run.log({key: wandb.Molecule(open(pdb)) for key, pdb in named_pdbs.items()})


_CONSENSUS_WEIGHT = {"trust": 1.0, "weak": 0.5, "defer": 0.0}


def _candidate_metrics(result: dict) -> dict:
    """Per-candidate scalars for a stepped log: reward, the consensus weight (trust=1 /
    weak=0.5 / defer=0), and each model's score under `score/<stage>` — so every
    verifier's verdict on the candidate plots as its own series across the run."""
    m = {"reward": result.get("reward", 0.0), "consensus_weight": _CONSENSUS_WEIGHT.get(result.get("consensus"), 0.0)}
    for e in result.get("stack", []):
        if e.get("status") == "ran":
            m[f"score/{e['stage']}"] = e["score"]
    return m


def log_candidates(run, results: list[dict]) -> None:
    """Step through candidates: log each at its own W&B step as a captioned 3D molecule
    plus every model's score for it. The media panel becomes a per-candidate slider, and
    each verifier's `score/<stage>` plots across the series — the cross-model comparison."""
    import wandb

    for i, r in enumerate(results):
        metrics = _candidate_metrics(r)
        src = r.get("structure")
        if src and Path(src).exists():
            metrics["candidate"] = wandb.Molecule(open(src), caption=Path(src).stem)
        run.log(metrics, step=i)


def _stack_cell(entry: dict) -> str:
    """One verifier's cell in the stack table: 'label score' if it ran, else its status."""
    return f"{entry['label']} {entry['score']:.2f}" if entry.get("status") == "ran" else entry["status"]


def _stack_rows(results: list[dict]) -> tuple[list[str], list[list], Counter]:
    """Pure shaping for `log_stack`: (columns, rows, consensus counts) from a batch of
    verify_structure results (each carrying a `reward`). One row per design; one column
    per stack tier (driven by the result's own `stack` view, so new tiers appear for
    free), plus a stress-robustness summary when present."""
    stages = [e["stage"] for e in results[0]["stack"]] if results and results[0].get("stack") else []
    has_stress = any(r.get("stress") for r in results)
    cols = ["design", "CN", "donors", *stages, *(["stress_holds"] if has_stress else []), "consensus", "reward"]
    rows, counts = [], Counter()
    for r in results:
        counts[r.get("consensus", "defer")] += 1
        by_stage = {e["stage"]: e for e in r.get("stack", [])}
        row = [Path(r.get("structure", "?")).name, r.get("coordination_number", 0), "".join(r.get("donors", []))]
        row += [_stack_cell(by_stage[s]) if s in by_stage else "—" for s in stages]
        if has_stress:
            st = r.get("stress", {})
            row.append(f"{sum(v.get('label') != 'defer' for v in st.values())}/{len(st)}" if st else "—")
        row += [r.get("consensus", "defer"), r.get("reward", 0.0)]
        rows.append(row)
    return cols, rows, counts


def log_stack(run, results: list[dict]) -> Counter:
    """Log a full-stack verification batch: a per-design × per-tier verdict table, a
    consensus-distribution bar, and a per-design reward bar. `results` are
    verify_structure dicts, each with a `reward` (see reward.reward_from_result)."""
    import wandb

    cols, rows, counts = _stack_rows(results)
    rewards = wandb.Table(data=[[row[0], row[-1]] for row in rows], columns=["design", "reward"])
    run.log({
        "stack": wandb.Table(columns=cols, data=rows),
        "consensus": _verdict_bar(counts, "Consensus across the stack"),
        "reward_by_design": wandb.plot.bar(rewards, "design", "reward", title="Verifier reward per design"),
    })
    _set_verdict_summary(run, counts)
    return counts
