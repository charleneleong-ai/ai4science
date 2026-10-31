"""Map the two GLP1R ligand sites from the active-state structure.

Class B1 receptors bind their peptide agonists in two parts (the "two-domain"
model): the peptide's C-terminal half lies in a groove on the extracellular
domain, and its N-terminal half drops into the transmembrane bundle, which is
what actually triggers signalling. Challenge 3 scores binding first and
activation second, so the design needs both sites located precisely rather
than taken from a review.

This reads 6X18 (GLP-1 / GLP1R / Gs / Nb35) and reports, per receptor residue,
which half of the peptide it contacts.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Annotated

import gemmi
import typer

RECEPTOR_CHAIN = "R"
PEPTIDE_CHAIN = "P"

# UniProt P43220 topology. The extracellular domain is a soluble, independently
# folding domain; everything from TM1 on is membrane-embedded and only reachable
# through the loop vestibule.
SIGNAL_PEPTIDE = range(1, 24)
ECD = range(24, 140)
TRANSMEMBRANE = {
    1: range(140, 165),
    2: range(176, 202),
    3: range(228, 252),
    4: range(266, 291),
    5: range(306, 329),
    6: range(349, 371),
    7: range(384, 405),
}
EXTRACELLULAR_LOOPS = {1: range(202, 228), 2: range(291, 306), 3: range(371, 384)}

# GLP-1 numbering follows the proglucagon convention used in the PDB entry, so
# the peptide runs 7-36. Residues 7-13 are the activation segment that enters
# the bundle; 16-36 is the part the extracellular domain grips. 14/15 are the
# hinge and are deliberately left unassigned.
ACTIVATION_SEGMENT = range(7, 14)
GROOVE_SEGMENT = range(16, 37)

CONTACT_CUTOFF = 4.5

FIELDNAMES = ("residue", "name", "code", "region", "touches_activation", "touches_groove")

app = typer.Typer(add_completion=False)


def heavy_atoms(residue: gemmi.Residue) -> list[gemmi.Atom]:
    return [a for a in residue if a.element.atomic_number != 1]


def min_distance(a: gemmi.Residue, b: gemmi.Residue) -> float:
    return min(
        (x.pos.dist(y.pos) for x in heavy_atoms(a) for y in heavy_atoms(b)),
        default=float("inf"),
    )


def polymer_residues(model: gemmi.Model, chain_name: str) -> list[gemmi.Residue]:
    chain = model.find_chain(chain_name)
    if chain is None:
        raise typer.BadParameter(f"chain {chain_name} not in model")
    return list(chain.get_polymer())


def segment_label(peptide_number: int) -> str | None:
    if peptide_number in ACTIVATION_SEGMENT:
        return "activation"
    if peptide_number in GROOVE_SEGMENT:
        return "groove"
    return None


def receptor_region(number: int) -> str:
    if number in ECD:
        return "ECD"
    for index, span in TRANSMEMBRANE.items():
        if number in span:
            return f"TM{index}"
    for index, span in EXTRACELLULAR_LOOPS.items():
        if number in span:
            return f"ECL{index}"
    return "other"


def centroid(residues: list[gemmi.Residue]) -> gemmi.Position:
    atoms = [a.pos for r in residues for a in heavy_atoms(r)]
    return gemmi.Position(
        sum(p.x for p in atoms) / len(atoms),
        sum(p.y for p in atoms) / len(atoms),
        sum(p.z for p in atoms) / len(atoms),
    )


@app.command()
def main(
    structure: Annotated[Path, typer.Option(help="6X18 mmCIF")],
    out: Annotated[Path, typer.Option(help="per-residue contact CSV")],
) -> None:
    model = gemmi.read_structure(str(structure))[0]
    receptor = polymer_residues(model, RECEPTOR_CHAIN)
    peptide = polymer_residues(model, PEPTIDE_CHAIN)

    rows: list[dict[str, object]] = []
    for res in receptor:
        contacts: dict[str, float] = {}
        for pep in peptide:
            label = segment_label(pep.seqid.num)
            if label is None:
                continue
            distance = min_distance(res, pep)
            if distance <= CONTACT_CUTOFF:
                contacts[label] = min(contacts.get(label, float("inf")), distance)
        if not contacts:
            continue
        rows.append(
            {
                "residue": res.seqid.num,
                "name": res.name,
                "code": gemmi.find_tabulated_residue(res.name).one_letter_code.upper(),
                "region": receptor_region(res.seqid.num),
                "touches_activation": round(contacts["activation"], 2) if "activation" in contacts else "",
                "touches_groove": round(contacts["groove"], 2) if "groove" in contacts else "",
            }
        )

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    activation = [r for r in rows if r["touches_activation"] != ""]
    groove = [r for r in rows if r["touches_groove"] != ""]

    typer.echo(f"{len(rows)} receptor residues within {CONTACT_CUTOFF} A of the peptide")
    for label, subset in (
        ("activation segment (GLP-1 7-13)", activation),
        ("groove segment (GLP-1 16-36)", groove),
    ):
        tally = Counter(str(row["region"]) for row in subset)
        regions = " ".join(f"{region}={count}" for region, count in sorted(tally.items()))
        typer.echo(f"  {label}: {len(subset)} residues  {regions}")
        typer.echo("    " + " ".join(f"{r['code']}{r['residue']}" for r in subset))

    shared = {r["residue"] for r in activation} & {r["residue"] for r in groove}
    typer.echo(f"  receptor residues shared by both halves: {sorted(shared) or 'none'}")

    # How far a design has to reach. The two sites are scored separately
    # (affinity on the groove, activation in the bundle) but a single molecule
    # has to span them, so the span is the binding constraint on the scaffold.
    by_number = {r.seqid.num: r for r in receptor}
    groove_site = [by_number[int(r["residue"])] for r in groove if int(r["residue"]) in ECD]
    bundle_site = [by_number[int(r["residue"])] for r in activation if str(r["region"]).startswith("TM")]
    span = centroid(groove_site).dist(centroid(bundle_site))
    typer.echo(f"  ECD groove centroid to TM pocket centroid: {span:.1f} A")

    peptide_by_number = {r.seqid.num: r for r in peptide}
    n_terminus = peptide_by_number[min(peptide_by_number)]
    spread, farthest = max(
        (min_distance(n_terminus, by_number[int(r["residue"])]), int(r["residue"])) for r in activation
    )
    typer.echo(
        f"  every activation contact lies within {spread:.1f} A of the peptide N-terminus"
        f" ({n_terminus.name}{n_terminus.seqid.num}), farthest at residue {farthest}"
    )

    # Insertion depth. The loop residues that touch the peptide form the rim of
    # the vestibule the N-terminus passes through, so the N-terminus' distance
    # below that rim is what a design's own N-terminal arm would have to cover.
    vestibule = [by_number[int(r["residue"])] for r in rows if str(r["region"]).startswith("ECL")]
    depth = centroid([n_terminus]).dist(centroid(vestibule))
    typer.echo(
        f"  N-terminus sits {depth:.1f} A below the centroid of the {len(vestibule)}"
        " loop residues forming the vestibule rim"
    )
    typer.echo(f"wrote {out}")


if __name__ == "__main__":
    app()
