"""Track a metal's coordination across an xtb MD trajectory (thermal-stability check).

Parses a multi-frame xyz trajectory (xtb.trj), and for each frame counts the donor
atoms (N/O/S) within `cutoff` of the metal — answering whether the coordination
survives 300 K dynamics or the metal dissociates.

    uv run python scripts/md_coordination.py <xtb.trj> --metal Ni --cutoff 2.8
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import typer
from rich.console import Console

console = Console()
_DONORS = {"N", "O", "S"}


def main(trajectory: str, metal: str = "Ni", cutoff: float = 2.8, dump_fs: float = 20.0) -> None:
    lines = Path(trajectory).read_text().strip().splitlines()
    cns, closest = [], []
    i = 0
    while i < len(lines):
        n = int(lines[i].split()[0])
        atoms = [lines[i + 2 + j].split() for j in range(n)]
        i += 2 + n
        els = [a[0] for a in atoms]
        xyz = np.array([[float(a[1]), float(a[2]), float(a[3])] for a in atoms])
        m = xyz[els.index(metal)]
        d = np.linalg.norm(np.array([xyz[k] for k, e in enumerate(els) if e in _DONORS]) - m, axis=1)
        cns.append(int((d <= cutoff).sum()))
        closest.append(float(d.min()))

    cns, closest = np.array(cns), np.array(closest)
    console.print(f"[bold]{metal} coordination over {len(cns)} frames "
                  f"(~{len(cns) * dump_fs / 1000:.1f} ps)[/]")
    console.print(f"  CN: mean [green]{cns.mean():.2f}[/]  (min {cns.min()}, max {cns.max()})")
    console.print(f"  frames with CN≥4: {100 * (cns >= 4).mean():.0f}%   CN≥3: {100 * (cns >= 3).mean():.0f}%")
    console.print(f"  closest {metal}–donor bond: mean {closest.mean():.2f} Å, "
                  f"max {closest.max():.2f} Å  → {'never dissociates' if closest.max() < 3.0 else 'dissociates'}")


if __name__ == "__main__":
    typer.run(main)
