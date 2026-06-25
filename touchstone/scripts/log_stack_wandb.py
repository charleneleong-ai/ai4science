"""Log one design's full verifier stack to W&B — the adjudication, visualized.

Three views of the same design's metal site: the geometry-trusted structure, an
independent Boltz-2 co-fold, and a physics (xtb GFN-FF) optimization — as comparable
3D molecules plus each verifier's verdict. Disagreement resolved by physics.

    uv run --extra viz python scripts/log_stack_wandb.py <pack.pdb> <boltz.pdb> <opt.pdb>

Reads WANDB_API_KEY / WANDB_PROJECT from touchstone/.env (gitignored).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import typer
import wandb
from dotenv import load_dotenv

from touchstone import BinderDesign, GeometryVerifier, PDBReference, coordination_site_from_pdb

load_dotenv(Path(__file__).parent.parent / ".env")


def main(pack_pdb: str, boltz_pdb: str, opt_pdb: str, name: str = "design_4",
         metal_element: str = "NI", metal_label: str = "Ni2+") -> None:
    verifier = GeometryVerifier(PDBReference())
    stages = [
        ("LigandMPNN pack", "geometry oracle", pack_pdb),
        ("Boltz-2 co-fold", "independent co-fold", boltz_pdb),
        ("xtb GFN-FF opt", "physics", opt_pdb),
    ]

    run = wandb.init(project="touchstone", name=f"verifier-stack-{name}")
    table = wandb.Table(columns=["stage", "verifier", "CN", "verdict", "bonds (Å)", "detail"])
    mols = {}
    print(f"verifier stack — {name} ({metal_label})\n")
    for stage, role, pdb in stages:
        try:
            site = coordination_site_from_pdb(pdb, metal_element, metal_label)
            vd = verifier.verify(BinderDesign(stage, site, "x", float("nan")))
            flag = "TRUST" if vd.trust else ("DEFER" if vd.ood else "weak")
            bonds = np.round(np.sort(site.bond_lengths()), 2).tolist()
            table.add_data(stage, role, site.coordination_number, flag, str(bonds), vd.reason)
            print(f"  {stage:18} [{role:20}] CN{site.coordination_number} {flag:6} {bonds}")
        except ValueError as e:
            table.add_data(stage, role, 0, "no site", "[]", str(e))
            print(f"  {stage:18} [{role:20}] no metal site")
        mols[f"3d/{stage.replace(' ', '_')}"] = wandb.Molecule(open(pdb))

    run.log({"verifier_stack": table, **mols})
    print(f"\nlogged → {run.url}")
    run.finish()


if __name__ == "__main__":
    typer.run(main)
