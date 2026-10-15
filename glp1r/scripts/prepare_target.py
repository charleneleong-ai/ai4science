"""Cut the design inputs out of 6X18.

Binder design tools want a single-chain target plus a hotspot list, and the two
GLP1R sites need different cuts: the extracellular domain is a soluble domain
that can be handed to any binder pipeline, while the transmembrane bundle only
makes sense in the context of the whole receptor. This writes both, and the
hotspots measured by target_sites.py alongside them, so no downstream script
re-derives a residue list.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Annotated

import gemmi
import typer
from target_sites import ECD, PEPTIDE_CHAIN, RECEPTOR_CHAIN

# The assay construct (Sino Biological 13944-H18H-DD) is stated as 24-432.
ASSAY_CONSTRUCT = range(24, 433)

app = typer.Typer(add_completion=False)


def load_uniprot_sequence(path: Path) -> str:
    return json.loads(path.read_text())["sequence"]["value"]


def write_chain(model: gemmi.Model, chain_name: str, residues: range | None, out: Path) -> int:
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model("1"))
    chain = gemmi.Chain(chain_name)
    kept = 0
    for residue in model[chain_name].get_polymer():
        if residues is not None and residue.seqid.num not in residues:
            continue
        chain.add_residue(residue)
        kept += 1
    structure[0].add_chain(chain)
    structure.setup_entities()
    structure.write_pdb(str(out))
    return kept


def hotspots_by_site(contacts: Path, sequence: str) -> dict[str, list[int]]:
    """Split the measured contacts into the two sites.

    The structure and the UniProt entry have to agree on numbering or every
    hotspot downstream is silently off; the CSV carries the structure's own
    residue code, so the check is free.
    """
    groove: list[int] = []
    activation: list[int] = []
    with contacts.open() as handle:
        for row in csv.DictReader(handle):
            number = int(row["residue"])
            if sequence[number - 1] != row["code"]:
                raise typer.BadParameter(
                    f"residue {number} is {row['code']} in the structure but"
                    f" {sequence[number - 1]} in UniProt: numbering disagrees"
                )
            if row["touches_groove"]:
                groove.append(number)
            if row["touches_activation"]:
                activation.append(number)
    return {
        "groove_all": groove,
        "groove_ecd_only": [n for n in groove if n in ECD],
        "activation": activation,
    }


@app.command()
def main(
    structure: Annotated[Path, typer.Option(help="6X18 mmCIF")],
    contacts: Annotated[Path, typer.Option(help="CSV from target_sites.py")],
    uniprot: Annotated[Path, typer.Option(help="P43220 JSON from the UniProt REST API")],
    out_dir: Annotated[Path, typer.Option(help="where to write design inputs")],
) -> None:
    model = gemmi.read_structure(str(structure))[0]
    out_dir.mkdir(parents=True, exist_ok=True)

    receptor_n = write_chain(model, RECEPTOR_CHAIN, None, out_dir / "receptor.pdb")
    ecd_n = write_chain(model, RECEPTOR_CHAIN, ECD, out_dir / "ecd.pdb")
    peptide_n = write_chain(model, PEPTIDE_CHAIN, None, out_dir / "glp1_reference.pdb")

    sequence = load_uniprot_sequence(uniprot)
    fasta = out_dir / "target_sequences.fasta"
    records = {
        f"GLP1R_assay_construct_{ASSAY_CONSTRUCT.start}_{ASSAY_CONSTRUCT.stop - 1}": sequence[
            ASSAY_CONSTRUCT.start - 1 : ASSAY_CONSTRUCT.stop - 1
        ],
        f"GLP1R_ecd_{ECD.start}_{ECD.stop - 1}": sequence[ECD.start - 1 : ECD.stop - 1],
    }
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))

    sites = hotspots_by_site(contacts, sequence)
    (out_dir / "hotspots.json").write_text(json.dumps(sites, indent=2) + "\n")

    typer.echo(f"receptor.pdb          {receptor_n} residues (chain R as modelled)")
    typer.echo(f"ecd.pdb               {ecd_n} residues ({ECD.start}-{ECD.stop - 1})")
    typer.echo(f"glp1_reference.pdb    {peptide_n} residues (GLP-1, for geometry reference only)")
    for name, seq in records.items():
        typer.echo(f"{name}  {len(seq)} aa")
    for name, numbers in sites.items():
        typer.echo(f"hotspots {name:16s} n={len(numbers):2d} {numbers}")
    typer.echo(f"wrote {out_dir}")


if __name__ == "__main__":
    app()
