"""Scan every interface position for histidine substitutions with a switchable pKa.

This replaces mechanism-based reasoning with a measurement. Objective 3 needs binding at pH 7.4
and release at pH 6.0, which requires a histidine whose pKa lies BETWEEN them: neutral at 7.4,
cationic at 6.0. Both mechanisms reasoned about earlier push pKa the WRONG way --

  PROPKA on a real design put three interface histidines at pKa 3.15, 4.70 and 5.71, all far
  below the 6.0 trigger. Interface burial destabilises the protonated form, and so does an
  adjacent cation. A buried His near a lysine is therefore the LEAST switchable position, yet
  it scored highest under the earlier structural proxy.

So: substitute His at every interface position in turn, predict pKa for each, and keep only
those landing in the usable window. PROPKA is CPU-only, so this runs without competing for GPU.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import typer
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE_TO_ONE
from rich.console import Console

from tnf_binder.pka import run_propka
from tnf_binder.rosetta import HisMutator
from tnf_binder.structure import Complex
from tnf_binder.tables import write_rows
from tnf_binder.verifiers.pka_switch import IDEAL_WINDOW, USABLE_WINDOW, pka_window

app = typer.Typer(add_completion=False)
console = Console()

MIN_INTERFACE_DSASA = 15.0
# Gly and Pro set backbone geometry and Cys may be disulphide-bonded, so none is a free slot.
KEEP = frozenset("GPC")
MARK = {"ideal": "[green]<== IDEAL[/]", "usable": "[yellow]<-- usable[/]", "dead": ""}
COLUMNS = ["parent", "chain", "resnum", "wt", "dsasa", "pka", "verdict", "switchable"]

InterfaceSite = tuple[str, int, str, float]


def interface_positions(complex_: Complex, min_dsasa: float) -> list[InterfaceSite]:
    """(chain, residue number, wild-type aa, buried area) for substitutable binder residues."""
    _target_buried, binder_buried = complex_.buried_area()
    found: list[InterfaceSite] = []
    for chain in complex_.binder_chains:
        for residue in (r for r in chain if r.id[0] == " "):
            buried = binder_buried.get((chain.id, residue.id[1]), 0.0)
            aa = THREE_TO_ONE.get(residue.get_resname(), "X")
            if buried >= min_dsasa and aa not in KEEP and aa not in ("H", "X"):
                found.append((chain.id, residue.id[1], aa, round(buried, 1)))
    found.sort(key=lambda site: -site[3])
    return found


def scan_design(source: Path, sites: list[InterfaceSite], mutator: HisMutator,
                workdir: Path) -> list[dict[str, object]]:
    """One row per substitution whose pKa PROPKA could predict."""
    rows: list[dict[str, object]] = []
    for chain, resnum, wild_type, dsasa in sites:
        mutant = workdir / f"{source.stem}_{chain}{resnum}H.pdb"
        pose = mutator.load(source)
        if not mutator.mutate(pose, [(chain, resnum)]):
            continue
        mutator.write(pose, mutant)
        pka = run_propka(mutant, workdir).get((chain, resnum))
        if pka is None:
            console.print(f"   {chain}{resnum} {wild_type}->H  propka gave no pKa")
            continue
        verdict = pka_window(pka)
        rows.append({"parent": source.stem, "chain": chain, "resnum": resnum,
                     "wt": wild_type, "dsasa": dsasa, "pka": round(pka, 2),
                     "verdict": verdict, "switchable": verdict in ("ideal", "usable")})
        console.print(f"   {chain}{resnum} {wild_type}->H  dSASA={dsasa:5.1f}  "
                      f"pKa={pka:5.2f}  {verdict:6s} {MARK[verdict]}")
    return rows


@app.command()
def main(
    structures: list[Path] = typer.Argument(..., help=".pdb complexes, or directories of them"),
    out: Path = typer.Option(..., "--out", help="CSV of every scanned substitution"),
    max_positions: int = typer.Option(14, "--max-positions",
                                      help="most-buried positions to scan per design"),
) -> None:
    """Substitute His at each interface position in turn and predict its pKa."""
    paths: list[Path] = []
    for argument in structures:
        paths.extend(sorted(argument.rglob("*.pdb")) if argument.is_dir() else [argument])
    if not paths:
        typer.echo("no .pdb inputs (propka needs PDB, not mmCIF)", err=True)
        raise typer.Exit(code=2)

    mutator = HisMutator()
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        for source in paths:
            complex_ = Complex.from_file(source)
            if not complex_.is_resolved:
                console.print(f"{source.stem}: could not separate binder from target")
                continue
            sites = interface_positions(complex_, MIN_INTERFACE_DSASA)[:max_positions]
            console.print(f"\n[bold]{source.stem}[/]: {len(sites)} interface positions to scan")
            rows += scan_design(source, sites, mutator, workdir)

    if not rows:
        typer.echo("no pKa results", err=True)
        raise typer.Exit(code=1)
    write_rows(out, rows, COLUMNS)
    ideal = sum(1 for r in rows if r["verdict"] == "ideal")
    usable = sum(1 for r in rows if r["verdict"] == "usable")
    console.print(f"\nscanned {len(rows)} substitutions: [green]{ideal}[/] ideal "
                  f"({IDEAL_WINDOW[0]}-{IDEAL_WINDOW[1]}), [yellow]{usable}[/] usable "
                  f"({USABLE_WINDOW[0]}-{USABLE_WINDOW[1]}) -> [cyan]{out}[/]")


if __name__ == "__main__":
    app()
