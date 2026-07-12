"""Validate the full verifier stack on a batch of structures and log it to W&B.

Runs verify_structure(deep, stress) over every matched structure, attaches the scalar
reward, and logs a per-design × per-tier verdict table + consensus/reward bars to a
single W&B run. The auditable end-to-end picture of the stack on real generator output.

    conda run -n mlip python scripts/wandb_stack_eval.py \
        'bg_motif_pdbs/*.pdb' 'boltzgen_out16/intermediate_designs/*.cif' --deep --stress
"""

from __future__ import annotations

import glob
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import BarColumn, MofNCompleteColumn, Progress, TextColumn, TimeElapsedColumn

from touchstone import tracking
from touchstone.reward import reward_from_result
from touchstone.cofold import cif_provider
from touchstone.service import mlip_backbone, verify_structure

console = Console()


def cofold_provider(paths: list[str], cofold_dir: str, metal: str):
    """Map each design to its predicted structure in `cofold_dir` (matched by filename
    stem — e.g. ni_motif_02.pdb → cofold_dir/ni_motif_02_m3d.pdb from AllMetal3D/Chai-1)
    and build the CofoldCrossCheck provider over it."""
    predictions = {}
    for p in paths:
        hits = sorted(glob.glob(f"{cofold_dir}/{Path(p).stem}*"))
        if hits:
            predictions[p] = hits[0]
    metal_atom = "".join(c for c in metal if c.isalpha()).upper()
    return cif_provider(predictions, metal_atom=metal_atom, metal=metal), len(predictions)


def main(
    patterns: list[str] = typer.Argument(..., help="one or more glob patterns of structures (.pdb/.cif)"),
    metal: str = "Ni2+",
    deep: bool = True,
    stress: bool = True,
    cofold_dir: str = typer.Option("", help="dir of independent predictions (AllMetal3D/Chai-1) for the co-fold tier"),
    name: str = "stack-eval",
    project: str = "touchstone",
) -> None:
    paths = sorted({p for pat in patterns for p in glob.glob(pat)})
    if not paths:
        raise typer.BadParameter(f"no structures match {patterns!r}")

    cofold = None
    if cofold_dir:
        cofold, n = cofold_provider(paths, cofold_dir, metal)
        console.log(f"co-fold tier on: {n}/{len(paths)} designs have a prediction in {cofold_dir}")

    calc = None
    if deep:  # build the MLIP backbone once (a ~minute model load), share across the batch
        with console.status("loading MLIP backbone (MACE-MP)…"):
            calc = mlip_backbone()
        console.log("MLIP backbone ready" if calc else "no MLIP backend — geometry tiers only")

    results = []
    with Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(), MofNCompleteColumn(), TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("verifying", total=len(paths))
        for p in paths:
            progress.update(task, description=f"verifying {Path(p).name}")
            try:
                r = verify_structure(p, metal, deep=deep, stress=stress, cofold_provider=cofold, calc=calc)
            except Exception as e:  # unparseable / no metal ⇒ worst reward, recorded, batch continues
                r = {"structure": p, "consensus": "defer", "verifiers": {}, "stack": [], "error": f"{type(e).__name__}: {e}"}
                progress.console.log(f"[red]{Path(p).name}: {r['error']}")
            r["reward"] = reward_from_result(r)
            results.append(r)
            progress.console.log(f"{Path(p).name}: {r['consensus']} (reward {r['reward']})")
            progress.advance(task)

    with console.status("logging to W&B…"):
        cfg = {"n": len(paths), "metal": metal, "deep": deep, "stress": stress, "patterns": patterns}
        run = tracking.init(name, config=cfg, project=project)
        tracking.log_candidates(run, results)  # per-candidate step: 3D molecule + each model's score
        counts = tracking.log_stack(run, results)  # summary table + consensus/reward bars
        run.finish()
    console.log(f"logged {dict(counts)} → {run.url}")


if __name__ == "__main__":
    typer.run(main)
