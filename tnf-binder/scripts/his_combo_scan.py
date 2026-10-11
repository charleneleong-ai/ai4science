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

from tnf_binder.pka import run_propka
from tnf_binder.rosetta import HisMutator
from tnf_binder.structure import Complex
from tnf_binder.tables import write_rows
from tnf_binder.verifiers.pka_switch import pka_window

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


def combination_row(parent: Path, combo: tuple[tuple[str, int], ...], mutant: Path,
                    workdir: Path) -> dict[str, object]:
    tag = "_".join(f"{chain}{num}H" for chain, num in combo)
    pkas = run_propka(mutant, workdir)
    values = {f"{chain}{num}": pkas.get((chain, num)) for chain, num in combo}
    predicted = [v for v in values.values() if v is not None]
    in_window = [v for v in predicted if pka_window(v) in ("ideal", "usable")]
    return {
        "design": f"{parent.stem}__{tag}",
        "parent": parent.stem,
        "n_his": len(combo),
        "positions": " ".join(f"{chain}{num}" for chain, num in combo),
        "pkas": " ".join(f"{k}={v:.2f}" if v is not None else f"{k}=NA"
                         for k, v in values.items()),
        "n_in_window": len(in_window),
        "n_ideal": sum(1 for v in predicted if pka_window(v) == "ideal"),
        "sequence": Complex.from_file(mutant).binder_sequence,
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
            tag = "_".join(f"{chain}{num}H" for chain, num in combo)
            mutant = (keep_dir or workdir) / f"{parent.stem}__{tag}.pdb"
            pose = mutator.load(parent)
            if not mutator.mutate(pose, list(combo)):
                console.print(f"  {tag}: mutation failed")
                continue
            mutator.write(pose, mutant)
            rows.append(combination_row(parent, combo, mutant, workdir))
            console.print(f"  {tag:28s} pKa: {rows[-1]['pkas']:48s} "
                          f"in_window={rows[-1]['n_in_window']}/{len(combo)} "
                          f"ideal={rows[-1]['n_ideal']}")

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
