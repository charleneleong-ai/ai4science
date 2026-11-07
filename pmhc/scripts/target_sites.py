"""Measure what a binder can actually see of the KRAS peptide in the HLA-A*11:01 groove.

Challenge 5 ranks mutation specificity above affinity: a design has to tell VVVGADGVGK from
VVVGAGGVGK, which differ only at peptide position 6 (KRAS Gly12 to Asp). That discrimination is
only available to a binder if the position-6 side chain is solvent-exposed in the complex, so the
exposure is measured here rather than assumed. A buried carboxylate would make the first ranking
tier unreachable by any design and is worth knowing before a GPU is booked.

Run against 7OW6 (G12D) and 7OW5 (wild type), both with the T-cell receptor deleted.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi
import typer
from common import load, one_letter, residues_of, sasa_model, sidechain_area, write_rows

app = typer.Typer(add_completion=False)

# Shrake-Rupley reference areas for a residue free in solution, computed on an extended
# Gly-X-Gly tripeptide; relative exposure is the complex area over this.
TRIPEPTIDE_AREA = {
    "ALA": 106.0, "ARG": 248.0, "ASN": 157.0, "ASP": 163.0, "CYS": 135.0,
    "GLN": 198.0, "GLU": 194.0, "GLY": 84.0, "HIS": 184.0, "ILE": 169.0,
    "LEU": 164.0, "LYS": 205.0, "MET": 188.0, "PHE": 197.0, "PRO": 136.0,
    "SER": 130.0, "THR": 142.0, "TRP": 227.0, "TYR": 222.0, "VAL": 142.0,
}

@dataclass(frozen=True)
class Exposure:
    structure: str
    chain: str
    position: int
    resname: str
    total: float
    sidechain: float
    relative: float
    polar_sidechain: float

    @property
    def row(self) -> dict[str, str | int | float]:
        return {
            "structure": self.structure,
            "chain": self.chain,
            "position": self.position,
            "resname": self.resname,
            "sasa_total": round(self.total, 1),
            "sasa_sidechain": round(self.sidechain, 1),
            "relative_exposure": round(self.relative, 3),
            "sasa_polar_sidechain": round(self.polar_sidechain, 1),
        }


def carve(path: Path, keep: tuple[str, ...], out: Path) -> Path:
    """Drop every chain but the pMHC itself and write the result.

    Both depositions are ternary complexes: chains D and E are the T-cell receptor alpha and
    beta. SASA measured with the receptor in place describes the peptide buried under a TCR, not
    the surface a binder is offered, and understates position-6 exposure several-fold. The
    challenge names chains A-C with the receptor removed, and this is that structure.
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_ligands_and_waters()
    present = {ch.name for ch in st[0]}
    if missing := set(keep) - present:
        raise typer.BadParameter(f"{path} lacks chain(s) {sorted(missing)}")
    for name in present - set(keep):
        st[0].remove_chain(name)
    st.setup_entities()
    st.make_mmcif_document().write_file(str(out))
    return out


def peptide_chain(path: Path) -> str:
    """The shortest polymer chain is the 10-mer epitope; the others are HLA heavy and b2m."""
    st = load(path)
    polymers = [(len(ch), ch.name) for ch in st[0] if len(ch) > 1]
    if not polymers:
        raise typer.BadParameter(f"{path} has no polymer chains")
    return min(polymers)[1]


def measure(path: Path, chain_id: str, label: str) -> list[Exposure]:
    """Per-residue SASA of one chain inside the assembled complex."""
    model = sasa_model(path)
    out: list[Exposure] = []
    for index, residue in enumerate(residues_of(model, chain_id)):
        name = residue.get_resname()
        total = sum(a.sasa for a in residue)
        side = sidechain_area(residue)
        reference = TRIPEPTIDE_AREA.get(name)
        out.append(
            Exposure(
                structure=label,
                chain=chain_id,
                position=index + 1,
                resname=name,
                total=total,
                sidechain=side,
                relative=total / reference if reference else float("nan"),
                polar_sidechain=sidechain_area(residue, polar_only=True),
            )
        )
    return out


@app.command()
def main(
    data: Path = typer.Option(Path("data"), help="directory holding the downloaded mmCIFs"),
    out: Path = typer.Option(Path("peptide_exposure.csv")),
    mutant: str = typer.Option("7OW6", help="G12D complex"),
    wildtype: str = typer.Option("7OW5", help="wild-type complex"),
    keep: str = typer.Option("A,B,C", help="chains to retain: HLA heavy, b2m, peptide"),
) -> None:
    chains = tuple(c.strip() for c in keep.split(","))
    rows: list[Exposure] = []
    for pdb_id in (mutant, wildtype):
        path = carve(data / f"{pdb_id}.cif", chains, data / f"{pdb_id}_pmhc.cif")
        chain_id = peptide_chain(path)
        measured = measure(path, chain_id, pdb_id)
        rows.extend(measured)
        seq = "".join(one_letter(r.resname) for r in measured)
        typer.echo(f"{pdb_id} peptide chain {chain_id}: {seq}")
        for r in measured:
            typer.echo(
                f"  {r.position:2d} {r.resname:3s} "
                f"total={r.total:6.1f} side={r.sidechain:6.1f} "
                f"rel={r.relative:5.2f} polar_side={r.polar_sidechain:5.1f}"
            )
    typer.echo(f"wrote {out} ({write_rows(rows, out)} rows)")


if __name__ == "__main__":
    app()
