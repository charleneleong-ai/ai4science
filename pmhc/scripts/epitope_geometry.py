"""How the one proven reader of this epitope discriminates it, measured rather than assumed.

Three questions decide whether a design can reach ranking tier 1, and all three are answerable
from the depositions:

1. Is the mutated side chain's conformation experimentally known without a receptor bound? The
   challenge page states 7OW4, the receptor-free G12D crystal, lacks coordinates for peptide
   positions 5 and 6. If so, the rotamer any design must read is modelled, not measured.
2. How much peptide surface is on offer at all, against the HLA surface a binder must also cover?
   This sets the floor on how much of any interface is spent on the shared, non-discriminating
   part of the target.
3. The TCR in 7OW6 does tell G12D from wild type (Poole et al., Nat Commun 2022), so the
   discrimination is achievable. Which of its residues touch position 6, from which direction,
   and over how much area? That is the geometric prescription, and it is measured here as a
   reference for where a de novo binder's own reading element has to sit. The receptor itself is
   not a template: the challenge requires de novo, zero-shot designs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from pathlib import Path

import gemmi
import numpy as np
import typer
from common import chain_area, load, residues_of, sasa_model, write_rows

app = typer.Typer(add_completion=False)

CARBOXYLATE = ("OD1", "OD2", "CG")
CONTACT_CUTOFF = 4.0


@dataclass
class Contact:
    tcr_chain: str
    tcr_resname: str
    tcr_seqid: int
    tcr_atom: str
    peptide_position: int
    peptide_resname: str
    peptide_atom: str
    distance: float

    @property
    def row(self) -> dict[str, str | int | float]:
        return {
            "tcr_chain": self.tcr_chain,
            "tcr_residue": f"{self.tcr_resname}{self.tcr_seqid}",
            "tcr_atom": self.tcr_atom,
            "peptide_position": self.peptide_position,
            "peptide_residue": f"{self.peptide_resname}{self.peptide_position}",
            "peptide_atom": self.peptide_atom,
            "distance": round(self.distance, 2),
        }


@dataclass
class Chain:
    """One polymer chain's residues, as (seqid, name, {atom: coordinate})."""

    name: str
    residues: list[tuple[int, str, dict[str, np.ndarray]]] = field(default_factory=list)

    @classmethod
    def read(cls, structure: gemmi.Structure, name: str) -> Chain:
        chain = cls(name=name)
        for residue in structure[0][name]:
            atoms = {a.name: np.array(a.pos.tolist()) for a in residue}
            chain.residues.append((residue.seqid.num, residue.name, atoms))
        return chain

    def coordinates(self) -> np.ndarray:
        return np.array([xyz for _, _, atoms in self.residues for xyz in atoms.values()])


def resolved_positions(path: Path, chain_name: str) -> list[tuple[int, str, int]]:
    """(seqid, resname, atom count) per residue, so unmodelled side chains are visible."""
    chain = Chain.read(load(path), chain_name)
    return [(seqid, name, len(atoms)) for seqid, name, atoms in chain.residues]


def contacts(path: Path, peptide: str, partners: tuple[str, ...], cutoff: float) -> list[Contact]:
    st = load(path)
    pep = Chain.read(st, peptide)
    found: list[Contact] = []
    for partner_name in partners:
        partner = Chain.read(st, partner_name)
        for (pos, (_, pep_name, pep_atoms)), (seqid, res_name, res_atoms) in product(
            enumerate(pep.residues, start=1), partner.residues
        ):
            for (pep_atom, pep_xyz), (atom, xyz) in product(pep_atoms.items(), res_atoms.items()):
                if (d := float(np.linalg.norm(pep_xyz - xyz))) <= cutoff:
                    found.append(
                        Contact(
                            tcr_chain=partner_name,
                            tcr_resname=res_name,
                            tcr_seqid=seqid,
                            tcr_atom=atom,
                            peptide_position=pos,
                            peptide_resname=pep_name,
                            peptide_atom=pep_atom,
                            distance=d,
                        )
                    )
    return found


def atom_sasa(model: object, chain_name: str, position: int) -> dict[str, float]:
    return {a.get_id(): round(a.sasa, 2) for a in residues_of(model, chain_name)[position - 1]}


@app.command()
def main(
    data: Path = typer.Option(Path("data")),
    out: Path = typer.Option(Path("tcr_readout.csv")),
    ternary: str = typer.Option("7OW6", help="G12D complex with the TCR present"),
    receptor_free: str = typer.Option("7OW4", help="G12D crystal without a TCR"),
    peptide: str = typer.Option("C"),
    tcr: str = typer.Option("D,E"),
    cutoff: float = typer.Option(CONTACT_CUTOFF),
) -> None:
    typer.echo(f"=== {receptor_free}: is the mutated side chain resolved without a receptor? ===")
    for seqid, name, natoms in resolved_positions(data / f"{receptor_free}.cif", peptide):
        typer.echo(f"  seqid {seqid:3d} {name:3s} atoms={natoms}")

    typer.echo(f"\n=== {ternary}: surface on offer, receptor removed ===")
    # One parse and one Shrake-Rupley pass for all four numbers below; this previously ran
    # four full passes over a ~4,800-atom complex.
    model = sasa_model(data / f"{ternary}_pmhc.cif")
    pep_area = chain_area(model, peptide)
    hla_area = chain_area(model, "A") + chain_area(model, "B")
    typer.echo(f"  peptide accessible   {pep_area:8.1f} A^2")
    typer.echo(f"  HLA + b2m accessible {hla_area:8.1f} A^2")
    typer.echo(f"  peptide share of the complex surface {100 * pep_area / (pep_area + hla_area):.1f}%")
    typer.echo(f"  position 6 per-atom SASA: {atom_sasa(model, peptide, 6)}")

    typer.echo(f"\n=== {ternary}: how the proven receptor reads position 6 ===")
    partners = tuple(c.strip() for c in tcr.split(","))
    found = contacts(data / f"{ternary}.cif", peptide, partners, cutoff)
    per_position: dict[int, int] = {}
    for c in found:
        per_position[c.peptide_position] = per_position.get(c.peptide_position, 0) + 1
    for pos in sorted(per_position):
        typer.echo(f"  peptide position {pos:2d}: {per_position[pos]:3d} receptor contacts")
    six = [c for c in found if c.peptide_position == 6]
    typer.echo(f"  contacts to position 6: {len(six)}")
    for c in sorted(six, key=lambda c: c.distance):
        flag = " <- carboxylate" if c.peptide_atom in CARBOXYLATE else ""
        typer.echo(
            f"    {c.tcr_chain} {c.tcr_resname}{c.tcr_seqid:<4d} {c.tcr_atom:4s}"
            f" -- {c.peptide_atom:4s} {c.distance:5.2f} A{flag}"
        )

    typer.echo(f"\nwrote {out} ({write_rows(found, out)} contacts)")


if __name__ == "__main__":
    app()
