"""Does G12D change the peptide's shape, or only its charge?

The receptor in 7OW6 discriminates G12D from wild type, yet at 4 A it contacts position 6 only
through that residue's backbone carbonyl, not its carboxylate. Two mechanisms fit that
observation and they imply different designs:

- charge readout: the design needs a basic group aimed at OD1, the one accessible carboxylate
  oxygen, over 6.4 A^2 of surface;
- shape readout: the mutation reshapes the peptide's central bulge, and the design reads a
  backbone conformation, over far more surface and without depending on a rotamer that no
  receptor-free crystal resolves.

Superposing the two complexes on the HLA heavy chain, which is chemically identical in both,
leaves any peptide difference as the signal. Deviations are reported per position so a local
change at the bulge is distinguishable from a rigid shift of the whole peptide.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np
import typer
from common import load, write_rows

app = typer.Typer(add_completion=False)

BACKBONE = ("N", "CA", "C", "O")


@dataclass(frozen=True)
class Deviation:
    position: int
    mutant_resname: str
    wildtype_resname: str
    backbone_rmsd: float
    ca_shift: float
    carbonyl_shift: float

    @property
    def row(self) -> dict[str, str | int | float]:
        return {
            "position": self.position,
            "mutant": self.mutant_resname,
            "wildtype": self.wildtype_resname,
            "backbone_rmsd": round(self.backbone_rmsd, 3),
            "ca_shift": round(self.ca_shift, 3),
            "carbonyl_o_shift": round(self.carbonyl_shift, 3),
        }


def ca_by_seqid(structure: gemmi.Structure, chain: str) -> dict[int, gemmi.Position]:
    out: dict[int, gemmi.Position] = {}
    for residue in structure[0][chain]:
        if atom := residue.find_atom("CA", "*"):
            out[residue.seqid.num] = atom.pos
    return out


def superpose_on(
    mobile: gemmi.Structure, fixed: gemmi.Structure, chain: str
) -> tuple[gemmi.Transform, float]:
    """Least-squares fit of mobile onto fixed over the CA atoms they share in `chain`."""
    a, b = ca_by_seqid(mobile, chain), ca_by_seqid(fixed, chain)
    shared = sorted(set(a) & set(b))
    if len(shared) < 50:
        raise typer.BadParameter(f"only {len(shared)} shared CA atoms in chain {chain}")
    fit = gemmi.superpose_positions([b[s] for s in shared], [a[s] for s in shared])
    return fit.transform, fit.rmsd


def atoms_by_position(structure: gemmi.Structure, chain: str) -> list[tuple[str, dict[str, np.ndarray]]]:
    out: list[tuple[str, dict[str, np.ndarray]]] = []
    for residue in structure[0][chain]:
        out.append((residue.name, {a.name: np.array(a.pos.tolist()) for a in residue}))
    return out


def shift(a: dict[str, np.ndarray], b: dict[str, np.ndarray], name: str) -> float:
    if name not in a or name not in b:
        return float("nan")
    return float(np.linalg.norm(a[name] - b[name]))


def backbone_rmsd(a: dict[str, np.ndarray], b: dict[str, np.ndarray]) -> float:
    shared = [n for n in BACKBONE if n in a and n in b]
    if not shared:
        return float("nan")
    diffs = np.array([a[n] - b[n] for n in shared])
    return float(np.sqrt((diffs**2).sum() / len(shared)))


@app.command()
def main(
    data: Path = typer.Option(Path("data")),
    out: Path = typer.Option(Path("conformer_deviation.csv")),
    mutant: str = typer.Option("7OW6_pmhc.cif"),
    wildtype: str = typer.Option("7OW5_pmhc.cif"),
    hla: str = typer.Option("A", help="chain to superpose on"),
    peptide: str = typer.Option("C"),
) -> None:
    mut, wt = load(data / mutant), load(data / wildtype)
    transform, fit_rmsd = superpose_on(wt, mut, hla)
    for chain in wt[0]:
        for residue in chain:
            for atom in residue:
                atom.pos = gemmi.Position(transform.apply(atom.pos))
    typer.echo(f"superposed on chain {hla}: {fit_rmsd:.3f} A over the shared CA atoms")

    rows: list[Deviation] = []
    for position, ((mut_name, mut_atoms), (wt_name, wt_atoms)) in enumerate(
        zip(atoms_by_position(mut, peptide), atoms_by_position(wt, peptide)), start=1
    ):
        rows.append(
            Deviation(
                position=position,
                mutant_resname=mut_name,
                wildtype_resname=wt_name,
                backbone_rmsd=backbone_rmsd(mut_atoms, wt_atoms),
                ca_shift=shift(mut_atoms, wt_atoms, "CA"),
                carbonyl_shift=shift(mut_atoms, wt_atoms, "O"),
            )
        )

    typer.echo("\n pos  mut  wt   bb_rmsd   CA     C=O")
    for r in rows:
        typer.echo(
            f"  {r.position:2d}  {r.mutant_resname:3s} {r.wildtype_resname:3s}"
            f"  {r.backbone_rmsd:6.3f}  {r.ca_shift:6.3f}  {r.carbonyl_shift:6.3f}"
        )

    typer.echo(f"\nwrote {out} ({write_rows(rows, out)} rows)")


if __name__ == "__main__":
    app()
