"""Extract a metal-coordination cluster from a structure for QM (xtb / cluster GFN2-xTB).

Pulls the residues near the metal + the metal itself into a small cluster PDB, and
reports the net charge (Ni²⁺ + basic residues − deprotonated carboxylates) — the input
a cluster GFN2-xTB optimization needs to sharpen the geometry verdict on a single site.

    uv run python scripts/extract_cluster.py <design.pdb> <out_cluster.pdb> --radius 5.0

Then, on a box with xtb + openbabel:
    obabel out_cluster.pdb -O cluster.xyz -p 7.4
    xtb cluster.xyz --gfn 2 --opt --chrg <reported charge>   # ulimit -s unlimited
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import typer


def main(pdb: str, out: str, metal_element: str = "NI", radius: float = 5.0) -> None:
    lines = [l for l in Path(pdb).read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])  # noqa: E731
    metal = next(xyz(l) for l in lines if l[76:78].strip().upper() == metal_element.upper())

    keep = {(l[21], l[22:26]) for l in lines if np.linalg.norm(xyz(l) - metal) <= radius}
    cluster = [l for l in lines if (l[21], l[22:26]) in keep]
    Path(out).write_text("\n".join(cluster) + "\nEND\n")

    res = [l[17:20].strip() for l in cluster if l[12:16].strip() == "CA"]
    c = Counter(res)
    charge = 2 + (c["ARG"] + c["LYS"]) - (c["ASP"] + c["GLU"])  # Ni²⁺ + basic − acidic (deprotonated)
    print(f"cluster: {len(res)} residues, {len(cluster)} heavy atoms, est net charge {charge}")
    print(f"residues: {sorted(res)}  →  wrote {out}")


if __name__ == "__main__":
    typer.run(main)
