"""Log one design's full verifier stack to W&B — the adjudication, visualized.

Three views of the same design's metal site: the geometry-trusted structure, an
independent Boltz-2 co-fold, and a physics (xtb) optimization — as comparable 3D
molecules plus each verifier's verdict. The W&B wiring lives in `touchstone.tracking`.

    uv run --extra viz python scripts/log_stack_wandb.py <pack.pdb> <boltz.pdb> <opt.pdb>
"""

from __future__ import annotations

import numpy as np
import typer
import wandb
from rich.console import Console
from rich.table import Table

from touchstone import BinderDesign, GeometryVerifier, PDBReference, coordination_site_from_pdb, tracking

console = Console()
STYLE = {"trust": "green", "weak": "yellow", "defer": "red", "no site": "red"}


def main(pack_pdb: str, boltz_pdb: str, opt_pdb: str, name: str = "design_4",
         metal_element: str = "NI", metal_label: str = "Ni2+") -> None:
    verifier = GeometryVerifier(PDBReference())
    stages = [
        ("LigandMPNN pack", "geometry oracle", pack_pdb),
        ("Boltz-2 co-fold", "independent co-fold", boltz_pdb),
        ("xtb GFN-FF opt", "physics", opt_pdb),
    ]

    run = tracking.init(f"verifier-stack-{name}")
    wb_table = wandb.Table(columns=["stage", "verifier", "CN", "verdict", "bonds (Å)", "detail"])
    rich_table = Table(title=f"verifier stack — {name} ({metal_label})")
    for col in ("stage", "verifier", "CN", "verdict", "bonds (Å)"):
        rich_table.add_column(col)

    mols = {}
    for stage, role, pdb in stages:
        try:
            site = coordination_site_from_pdb(pdb, metal_element, metal_label)
            vd = verifier.verify(BinderDesign(stage, site, "x", float("nan")))
            bonds = str(np.round(np.sort(site.bond_lengths()), 2).tolist())
            wb_table.add_data(stage, role, site.coordination_number, vd.label, bonds, vd.reason)
            rich_table.add_row(stage, role, str(site.coordination_number), vd.label, bonds,
                               style=STYLE.get(vd.label))
        except ValueError as e:
            wb_table.add_data(stage, role, 0, "no site", "[]", str(e))
            rich_table.add_row(stage, role, "0", "no site", "[]", style="red")
        mols[f"3d/{stage.replace(' ', '_')}"] = pdb

    run.log({"verifier_stack": wb_table})
    tracking.log_molecules(run, mols)
    console.print(rich_table)
    console.print(f"[bold]logged →[/] {run.url}")
    run.finish()


if __name__ == "__main__":
    typer.run(main)
