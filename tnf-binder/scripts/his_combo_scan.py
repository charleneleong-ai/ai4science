"""Build multi-histidine pH switches and re-predict pKa for the COMBINATION.

Single-mutant pKa values do not compose. Two histidines placed near each other shift one
another's pKa (a protonated imidazole destabilises its neighbour's protonation), so a set of
positions that each look ideal alone can collapse once combined. Multiple switchable His are
still worth pursuing -- they give a steeper, cooperative transition between pH 7.4 and 6.0 than
any single residue -- but the combination has to be measured, not inferred.

For each combination this reports every introduced histidine's pKa in the mutant context, and
counts how many remain inside the switch window.
"""

from __future__ import annotations

import itertools
import tempfile
from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.rosetta import HisMutator
from tnf_binder.switch import HisVariant, his_tag, measure_his, position_list
from tnf_binder.structure import Complex
from tnf_binder.tables import write_rows

app = typer.Typer(add_completion=False)
console = Console()

COLUMNS = ["design", "parent", "n_his", "positions", "pkas", "n_in_window", "n_ideal",
           "sequence"]


def parse_sites(positions: str) -> list[tuple[str, int]]:
    """'A2,A5,A8' -> [('A', 2), ('A', 5), ('A', 8)]."""
    sites = []
    for token in positions.split(","):
        token = token.strip()
        sites.append((token[0], int(token[1:])))
    return sites


def combination_row(variant: HisVariant) -> dict[str, object]:
    """One row per combination, carrying every introduced histidine's own measurement.

    pKa values do not compose -- A35H+A38H drops A35 from 6.07 to 5.94, out of the window --
    so what is reported is what the combination measured, never a sum of singles.
    """
    return {
        "design": variant.name,
        "parent": variant.parent,
        "n_his": len(variant.measurements),
        "positions": position_list(variant.positions),
        "pkas": " ".join(str(m) for m in variant.measurements),
        "n_in_window": variant.in_window,
        "n_ideal": variant.ideal,
        "sequence": Complex.from_file(variant.structure).binder_sequence,
    }


@app.command()
def main(
    parent: Path = typer.Argument(..., help="the parent complex, as .pdb"),
    positions: str = typer.Option(..., "--positions",
                                  help="comma list like A2,A5,A8,A15 (chain letter + resnum)"),
    out: Path = typer.Option(..., "--out", help="CSV of every combination measured"),
    max_his: int = typer.Option(4, "--max-his"),
    keep_dir: Path | None = typer.Option(None, "--keep-dir",
                                         help="keep the mutant PDBs here for refolding"),
) -> None:
    """Measure every His combination up to `max_his`, rather than composing single-mutant pKas."""
    sites = parse_sites(positions)
    combos = [c for k in range(1, max_his + 1) for c in itertools.combinations(sites, k)]
    console.print(f"{len(combos)} combinations from {len(sites)} positions")
    if keep_dir:
        keep_dir.mkdir(parents=True, exist_ok=True)

    mutator = HisMutator()
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        for combo in combos:
            variant = measure_his(parent, combo, mutator=mutator, workdir=workdir,
                                  dest=keep_dir)
            if variant is None:
                console.print(f"  {his_tag(combo)}: mutation failed")
                continue
            rows.append(combination_row(variant))
            console.print(f"  {his_tag(combo):28s} pKa: {rows[-1]['pkas']:48s} "
                          f"in_window={variant.in_window}/{len(combo)} "
                          f"ideal={variant.ideal}")

    if not rows:
        typer.echo("no results", err=True)
        raise typer.Exit(code=1)
    rows.sort(key=lambda r: (-r["n_ideal"], -r["n_in_window"], -r["n_his"]))
    write_rows(out, rows, COLUMNS)
    best = rows[0]
    console.print(f"\nbest: [bold]{best['design']}[/] -- {best['n_ideal']} ideal / "
                  f"{best['n_in_window']} in-window of {best['n_his']} His")
    console.print(f"wrote [cyan]{out}[/]")


if __name__ == "__main__":
    app()
