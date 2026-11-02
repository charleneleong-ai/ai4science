"""Can a designed side chain physically reach the one accessible carboxylate oxygen?

OD1 of the G12D aspartate is the only atom that distinguishes this target from the wild-type
complex and is also solvent-accessible, at 6.4 A^2. Accessible is not the same as reachable: an
atom at the bottom of a narrow pit between the HLA helices has surface area yet no approach
vector wide enough for an arginine. Rays are cast from the atom over a uniform sphere and the
ones that leave the complex without striking a heavy atom are counted, which gives the escape
cone an incoming side chain would have to use.

A cone of a few percent means a design would have to thread a needle; a broad cone means the
carboxylate can be engaged with ordinary geometry. The same probe is run on the wild-type
complex's position-6 backbone carbonyl, so the number has something to be compared against.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import typer
from common import HYDROGEN, load, write_rows

app = typer.Typer(add_completion=False)

# Heavy-atom vdW radii, in A; a probe radius is added on top to stand for the approaching atom.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
DEFAULT_VDW = 1.70


def fibonacci_directions(count: int) -> np.ndarray:
    """`count` roughly equally spaced unit vectors, so solid angle maps to a ray fraction."""
    index = np.arange(count) + 0.5
    phi = np.arccos(1 - 2 * index / count)
    theta = np.pi * (1 + 5**0.5) * index
    return np.column_stack(
        (np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi))
    )


@dataclass(frozen=True)
class Reach:
    structure: str
    atom: str
    free_fraction: float
    rays: int
    probe: float
    reach: float

    @property
    def row(self) -> dict[str, str | int | float]:
        return {
            "structure": self.structure,
            "atom": self.atom,
            "escape_cone_fraction": round(self.free_fraction, 7),
            "escape_cone_percent": round(100 * self.free_fraction, 5),
            "free_rays": round(self.free_fraction * self.rays),
            "rays": self.rays,
            "probe_radius": self.probe,
            "reach_distance": self.reach,
        }


def heavy_atoms(path: Path, exclude: tuple[str, int] | None) -> tuple[np.ndarray, np.ndarray]:
    """Coordinates and radii of every heavy atom, optionally dropping one residue's own atoms."""
    st = load(path)
    coords: list[list[float]] = []
    radii: list[float] = []
    for chain in st[0]:
        for position, residue in enumerate(chain, start=1):
            if exclude and chain.name == exclude[0] and position == exclude[1]:
                continue
            for atom in residue:
                if atom.element == HYDROGEN:
                    continue
                coords.append(atom.pos.tolist())
                radii.append(VDW.get(atom.element.name, DEFAULT_VDW))
    return np.array(coords), np.array(radii)


def atom_position(path: Path, chain: str, position: int, name: str) -> np.ndarray | None:
    """None when the residue has no such atom, so a control position can share one atom list."""
    residue = list(load(path)[0][chain])[position - 1]
    if not (atom := residue.find_atom(name, "*")):
        return None
    return np.array(atom.pos.tolist())


def escape_fraction(
    origin: np.ndarray,
    coords: np.ndarray,
    radii: np.ndarray,
    directions: np.ndarray,
    probe: float,
    reach: float,
) -> float:
    """Fraction of directions along which a probe sphere can travel `reach` A unobstructed.

    A ray is blocked when some atom's centre lies within (its radius + probe) of the ray and
    ahead of the origin, which is the standard ray-sphere test written for all atoms at once.

    Atoms too far away to block any ray are dropped first. A blocking atom has its centre
    within `reach` along the ray and within (radius + probe) across it, so its distance from
    the origin cannot exceed the hypotenuse of those two legs -- 8.67 A here, which keeps 115
    of 3,134 atoms and is exact, not an approximation. Without it the float64 intermediates
    are 308 MB per call and `--rays 262144` needs 20 GB, which is the ray count the escape
    cone would have to be re-measured at: the 0.02% quoted at the default 4096 rests on a
    single surviving ray.
    """
    offsets = coords - origin
    if len(offsets):
        limit = np.hypot(reach, radii.max() + probe)
        near = (offsets**2).sum(axis=1) <= limit**2
        offsets, radii = offsets[near], radii[near]
    along = offsets @ directions.T
    perpendicular_sq = (offsets**2).sum(axis=1)[:, None] - along**2
    threshold = (radii + probe)[:, None] ** 2
    blocked = (along > 0) & (along < reach) & (perpendicular_sq < threshold)
    return float((~blocked.any(axis=0)).mean())


@app.command()
def main(
    data: Path = typer.Option(Path("data")),
    out: Path = typer.Option(Path("reachability.csv")),
    mutant: str = typer.Option("7OW6_pmhc.cif"),
    wildtype: str = typer.Option("7OW5_pmhc.cif"),
    peptide: str = typer.Option("C"),
    position: int = typer.Option(6),
    atoms: str = typer.Option("OD1,O", help="atoms of that residue to probe"),
    rays: int = typer.Option(
        262144,
        help="directions sampled; the default resolves a cone of a few thousandths of a percent",
    ),
    probe: float = typer.Option(1.55, help="radius of the approaching atom; 1.55 is nitrogen"),
    reach: float = typer.Option(8.0, help="how far the probe must travel to count as escaped"),
) -> None:
    directions = fibonacci_directions(rays)
    rows: list[Reach] = []
    targets = [
        (filename, name)
        for filename in (mutant, wildtype)
        for name in (a.strip() for a in atoms.split(","))
    ]
    for filename, atom_name in targets:
        path = data / filename
        origin = atom_position(path, peptide, position, atom_name)
        if origin is None:
            typer.echo(f"{filename:16s} {peptide}{position}:{atom_name:4s} absent, skipped")
            continue
        coords, radii = heavy_atoms(path, exclude=(peptide, position))
        fraction = escape_fraction(origin, coords, radii, directions, probe, reach)
        rows.append(
            Reach(
                structure=filename,
                atom=f"{peptide}{position}:{atom_name}",
                free_fraction=fraction,
                rays=rays,
                probe=probe,
                reach=reach,
            )
        )
        typer.echo(
            f"{filename:16s} {peptide}{position}:{atom_name:4s}"
            f" escape cone {100 * fraction:8.4f}%  ({round(fraction * rays)} of {rays} rays"
            f" free, probe {probe} A, reach {reach} A)"
        )

    typer.echo(f"wrote {out} ({write_rows(rows, out)} rows)")


if __name__ == "__main__":
    app()
