"""Put BoltzGen's own confidence next to touchstone's independent verdict, per design.

BoltzGen scores its predictions (iPTM / pLDDT / pTM) — the generator grading itself.
touchstone re-judges the *chemistry* of the same structure, independent of that confidence.
This pairs a design's BoltzGen fold-output `.npz` (confidence arrays, one per refold sample)
with its refold CIF and prints both side by side — so you can spot a confident design with
a bad metal site (high iPTM, touchstone `defer`) or the reverse. With `--wandb` it logs the
comparison table + a consensus bar + an iPTM-vs-touchstone-reward scatter to Weights & Biases.

    uv run python scripts/boltzgen_scores.py \
        --npz-dir  boltzgen_out/intermediate_designs_inverse_folded/fold_out_npz \
        --cif-dir  boltzgen_out/intermediate_designs_inverse_folded/refold_cif \
        --metal Ni2+ --wandb

(run from the touchstone repo so `touchstone` is importable; needs the BoltzGen output dirs.)
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from touchstone.boltzgen import boltzgen_confidence
from touchstone.reward import reward_from_result
from touchstone.service import verify_structure

COLOR = {"trust": "green", "weak": "yellow", "defer": "red"}


def score_rows(npz_dir: Path, cif_dir: Path, metal: str) -> list[dict]:
    """One row per design: BoltzGen confidence (best refold sample) + touchstone verdict."""
    rows = []
    for npz_path in sorted(npz_dir.glob("*.npz")):
        cif = cif_dir / f"{npz_path.stem}.cif"
        if not cif.exists():
            continue
        row = {"design": npz_path.stem, **boltzgen_confidence(npz_path)}
        try:
            r = verify_structure(cif, metal)
            row["consensus"] = r["consensus"]
            row["reward"] = reward_from_result(r)
            # why it isn't trusted: the reason from each tier that flagged it
            row["why"] = "; ".join(
                f"{n} {v['reason'].replace(' — defer', '')}"
                for n, v in r["verifiers"].items() if v.get("label") not in (None, "trust")
            ) or "—"
        except Exception as e:  # unparseable / no metal ⇒ record, don't crash the batch
            row["consensus"], row["reward"], row["why"] = "error", 0.0, type(e).__name__
        rows.append(row)
    return rows


def print(rows: list[dict], metal: str) -> None:
    console = Console()
    table = Table(title=f"BoltzGen confidence  vs  touchstone verdict · {metal}")
    for col in ("design", "iPTM↑", "pLDDT↑", "pTM↑", "touchstone", "reward", "why (touchstone)"):
        table.add_column(col)
    fmt = lambda x: f"{x:.2f}" if isinstance(x, float) else "—"
    for r in rows:
        c = r["consensus"]
        verdict = f"[{COLOR[c]}]{c.upper()}[/]" if c in COLOR else "[dim]error[/]"
        table.add_row(r["design"], fmt(r["iptm"]), fmt(r["plddt"]), fmt(r["ptm"]), verdict, fmt(r["reward"]), r["why"])
    console.print(table)
    console.print("[dim]BoltzGen iPTM/pLDDT/pTM = generator self-confidence (best refold sample); "
                  "touchstone = independent chemistry verdict on the same structure.[/]")


def to_wandb(rows: list[dict], metal: str, project: str | None) -> None:
    import wandb

    from touchstone import tracking

    run = tracking.init(f"boltzgen-vs-touchstone-{metal}", config={"metal": metal, "n": len(rows)}, project=project)
    tbl = wandb.Table(columns=["design", "iptm", "plddt", "ptm", "consensus", "reward", "why"])
    counts = {"trust": 0, "weak": 0, "defer": 0}
    for r in rows:
        counts[r["consensus"]] = counts.get(r["consensus"], 0) + 1
        tbl.add_data(r["design"], r["iptm"], r["plddt"], r["ptm"], r["consensus"], r["reward"], r["why"])
    bar = wandb.Table(data=[[k, counts.get(k, 0)] for k in ("trust", "weak", "defer")], columns=["verdict", "count"])
    run.log({
        "boltzgen_vs_touchstone": tbl,
        "consensus": wandb.plot.bar(bar, "verdict", "count", title="touchstone consensus over the BoltzGen pool"),
        "iptm_vs_reward": wandb.plot.scatter(tbl, "iptm", "reward", title="BoltzGen iPTM vs touchstone reward"),
    })
    for k, n in counts.items():
        run.summary[f"n_{k}"] = n
    run.finish()


def main(
    npz_dir: Path = typer.Option(..., help="BoltzGen fold_out_npz dir (confidence arrays per design)"),
    cif_dir: Path = typer.Option(..., help="matching refold_cif dir (structures to verify)"),
    metal: str = typer.Option("Ni2+", help="target metal, e.g. Ni2+ / Cu2+ / Co2+"),
    wandb: bool = typer.Option(False, "--wandb", help="also log the table + plots to Weights & Biases"),
    project: str = typer.Option(None, help="W&B project (default: WANDB_PROJECT or 'touchstone')"),
) -> None:
    rows = score_rows(npz_dir, cif_dir, metal)
    print(rows, metal)
    if wandb:
        to_wandb(rows, metal, project)


if __name__ == "__main__":
    typer.run(main)
