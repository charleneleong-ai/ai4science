"""The two stages of engineering a pH switch: find the positions, then test them together.

    his.py scan   <complexes> --out scan.csv        # which positions could switch at all
    his.py combos <parent> --positions A35,A72      # what they measure when combined

Stage one replaces mechanism-based reasoning with a measurement, because the reasoning was
wrong three times over. Objective 3 needs a histidine whose pKa lies between the two assay
points -- neutral at 7.4, cationic at 6.0 -- and every structural proxy argued for pushed pKa
the wrong way: PROPKA put three interface histidines at 3.15, 4.70 and 5.71, all below the 6.0
trigger, because burial destabilises the protonated form and so does an adjacent cation. A
buried His beside a lysine is the LEAST switchable position, and it scored highest under the
proxy. So: substitute at every interface position in turn and predict each pKa.

Stage two exists because those pKas do not compose. A35H+A38H drops A35 from 6.07 to 5.94, out
of the window, so a combination has to be built and measured as a combination. The two stages
share `tnf_binder.switch.measure_his`; what differs is that one discovers positions from
structure and the other is told which to try.

PROPKA is CPU-only, so neither stage competes for the GPU.
"""

from __future__ import annotations

import itertools
import tempfile
from pathlib import Path

import typer
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE_TO_ONE
from rich.console import Console

from tnf_binder.rosetta import HisMutator
from tnf_binder.structure import Complex
from tnf_binder.switch import HisVariant, his_tag, measure_his, position_list
from tnf_binder.tables import write_rows
from tnf_binder.verifiers.pka_switch import IDEAL_WINDOW, USABLE_WINDOW

app = typer.Typer(add_completion=False)
console = Console()

MIN_INTERFACE_DSASA = 15.0
# Gly and Pro set backbone geometry and Cys may be disulphide-bonded, so none is a free slot.
KEEP = frozenset("GPC")
MARK = {"ideal": "[green]<== IDEAL[/]", "usable": "[yellow]<-- usable[/]", "dead": ""}
SCAN_COLUMNS = ["parent", "chain", "resnum", "wt", "dsasa", "pka", "verdict", "switchable"]
COMBO_COLUMNS = ["design", "parent", "n_his", "positions", "pkas", "n_in_window", "n_ideal",
                 "sequence"]

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
        variant = measure_his(source, [(chain, resnum)], mutator=mutator, workdir=workdir)
        if variant is None:
            continue
        measurement = variant.measurements[0]
        if measurement.pka is None:
            console.print(f"   {chain}{resnum} {wild_type}->H  propka gave no pKa")
            continue
        rows.append({"parent": source.stem, "chain": chain, "resnum": resnum,
                     "wt": wild_type, "dsasa": dsasa, "pka": round(measurement.pka, 2),
                     "verdict": measurement.window,
                     "switchable": measurement.in_window})
        console.print(f"   {chain}{resnum} {wild_type}->H  dSASA={dsasa:5.1f}  "
                      f"pKa={measurement.pka:5.2f}  {measurement.window:6s} "
                      f"{MARK[measurement.window]}")
    return rows

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
def scan(
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


@app.command()
def combos(
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
