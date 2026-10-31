"""Score a ligand's engagement of the GLP1R activation machinery.

No tool predicts agonism, and the published record says confidence scores do not
separate binders from non-binders, let alone agonists from antagonists. What did
discriminate, in the one campaign that produced de novo GPCR agonists, was
whether a design independently reproduced the natural ligand's activating
geometry. So this measures that geometry explicitly, feature by feature, with
cutoffs declared up front rather than fitted.

The pharmacophore is ranked by the mutational cost of losing each feature; see
docs/strategy.md for the sources. Everything is in UniProt P43220 numbering.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import gemmi
import typer

# Vestibule rim: the extracellular-loop residues that contact GLP-1 in 6X18.
# Depth is measured relative to their centroid, so "inside the bundle" is
# defined by the receptor itself rather than by a membrane we are not modelling.
VESTIBULE_RIM = (205, 210, 211, 214, 215, 221, 298, 299, 300, 372, 380)

# The floor of the orthosteric cavity. The rim-to-floor vector is the axis.
POCKET_FLOOR = (234, 237, 241, 306, 310, 313, 387)

ECD = range(24, 140)
INTRACELLULAR = (range(165, 176), range(252, 266), range(329, 349), range(405, 464))


@dataclass(frozen=True)
class Feature:
    """One pharmacophore element, with the receptor atoms it must reach."""

    name: str
    rank: int
    residues: tuple[int, ...]
    receptor_atoms: tuple[str, ...]
    ligand_atoms: tuple[str, ...]
    cutoff: float


# Cutoffs are 1.0 A looser than the distance the native complex achieves, so a
# design is credited for reaching the same partner without having to match a
# 2.1 A crystal structure atom for atom. They are declared here, before any
# design exists, and are not tuned afterwards.
FEATURES = (
    Feature("anion_to_R190", 1, (190,), ("NE", "NH1", "NH2"), ("OD1", "OD2", "OE1", "OE2"), 4.2),
    Feature("aromatic_at_Y148", 2, (148, 141, 144, 388), ("OH", "CZ", "CE1", "CE2"), ("CZ", "CE1", "CE2", "CD1", "CD2", "CG", "NE1", "CH2"), 4.5),
    Feature("small_at_W306", 3, (306,), ("CZ2", "CH2", "CZ3", "NE1"), (), 4.5),
    Feature("hbond_to_Q234", 4, (234,), ("OE1", "NE2"), ("N", "O", "OD1", "OE1", "ND1", "NE2", "OG", "OG1", "OH", "NZ"), 4.0),
    Feature("hydrophobe_at_E387", 5, (387, 388), ("OE1", "OE2", "CD1", "CD2"), (), 4.5),
    Feature("ecl3_contact", 6, (372, 373, 380), (), (), 4.5),
)

app = typer.Typer(add_completion=False)


def heavy(residue: gemmi.Residue) -> list[gemmi.Atom]:
    return [a for a in residue if a.element.atomic_number != 1]


def centroid(positions: list[gemmi.Position]) -> gemmi.Position:
    n = len(positions)
    return gemmi.Position(
        sum(p.x for p in positions) / n, sum(p.y for p in positions) / n, sum(p.z for p in positions) / n
    )


def unit(a: gemmi.Position, b: gemmi.Position) -> tuple[float, float, float]:
    dx, dy, dz = b.x - a.x, b.y - a.y, b.z - a.z
    norm = (dx * dx + dy * dy + dz * dz) ** 0.5
    return dx / norm, dy / norm, dz / norm


def selected_atoms(
    residues: dict[int, gemmi.Residue], numbers: tuple[int, ...], names: tuple[str, ...]
) -> list[gemmi.Atom]:
    out: list[gemmi.Atom] = []
    for number in numbers:
        residue = residues.get(number)
        if residue is None:
            continue
        out.extend(a for a in heavy(residue) if not names or a.name in names)
    return out


def min_pair_distance(left: list[gemmi.Atom], right: list[gemmi.Atom]) -> float:
    return min((a.pos.dist(b.pos) for a in left for b in right), default=float("inf"))


def receptor_residues(model: gemmi.Model, chain_name: str, offset: int) -> dict[int, gemmi.Residue]:
    """Index the receptor by UniProt number.

    A predicted complex numbers its chains from 1, so `offset` is the UniProt
    number of the chain's first residue. For an experimental GLP1R entry the
    author numbering already is UniProt numbering and offset is 0.
    """
    chain = model.find_chain(chain_name)
    if chain is None:
        raise typer.BadParameter(f"receptor chain {chain_name} not found")
    polymer = list(chain.get_polymer())
    if offset == 0:
        return {r.seqid.num: r for r in polymer}
    return {offset + i: r for i, r in enumerate(polymer)}


def depth_profile(
    ligand: list[gemmi.Residue], rim: gemmi.Position, axis: tuple[float, float, float]
) -> tuple[float, float]:
    """Deepest and shallowest ligand atom, projected onto the rim-to-floor axis."""
    projections = [
        (a.pos.x - rim.x) * axis[0] + (a.pos.y - rim.y) * axis[1] + (a.pos.z - rim.z) * axis[2]
        for residue in ligand
        for a in heavy(residue)
    ]
    return max(projections), min(projections)


def buried_split(
    receptor: dict[int, gemmi.Residue], ligand: list[gemmi.Residue], cutoff: float = 4.5
) -> dict[str, int]:
    ligand_atoms = [a for r in ligand for a in heavy(r)]
    tally = {"ecd": 0, "tmd": 0, "intracellular": 0}
    for number, residue in receptor.items():
        if min_pair_distance(heavy(residue), ligand_atoms) > cutoff:
            continue
        if number in ECD:
            tally["ecd"] += 1
        elif any(number in span for span in INTRACELLULAR):
            tally["intracellular"] += 1
        else:
            tally["tmd"] += 1
    return tally


def score_complex(
    model: gemmi.Model, receptor_chain: str, ligand_chain: str, offset: int
) -> dict[str, object]:
    receptor = receptor_residues(model, receptor_chain, offset)
    chain = model.find_chain(ligand_chain)
    if chain is None:
        raise typer.BadParameter(f"ligand chain {ligand_chain} not found")
    ligand = list(chain.get_polymer())
    ligand_atoms = [a for r in ligand for a in heavy(r)]

    rim_atoms = selected_atoms(receptor, VESTIBULE_RIM, ())
    floor_atoms = selected_atoms(receptor, POCKET_FLOOR, ())
    if not rim_atoms or not floor_atoms:
        raise typer.BadParameter("receptor is missing the rim or floor residues; check --offset")
    rim = centroid([a.pos for a in rim_atoms])
    axis = unit(rim, centroid([a.pos for a in floor_atoms]))

    deepest, shallowest = depth_profile(ligand, rim, axis)
    split = buried_split(receptor, ligand)

    row: dict[str, object] = {
        "ligand_len": len(ligand),
        "depth_max": round(deepest, 2),
        "depth_min": round(shallowest, 2),
        "contacts_ecd": split["ecd"],
        "contacts_tmd": split["tmd"],
        "contacts_intracellular": split["intracellular"],
    }

    satisfied = 0
    for feature in FEATURES:
        target = selected_atoms(receptor, feature.residues, feature.receptor_atoms)
        probe = [a for a in ligand_atoms if not feature.ligand_atoms or a.name in feature.ligand_atoms]
        distance = min_pair_distance(target, probe) if target and probe else float("inf")
        row[feature.name] = round(distance, 2) if distance < 99 else ""
        if distance <= feature.cutoff:
            satisfied += 1
    row["features_met"] = satisfied
    row["tier"] = assign_tier(deepest, split, satisfied)
    return row


# Set from the control ladder in docs/calibration.md, not from any design: the
# native peptide truncated to its groove-binding half reaches 7.3 A below the rim
# centroid with 0/6 features, the experimentally determined weak partial agonist
# GLP-1(9-36) reaches 14.2 A with 1/6, and the three full agonists span 18.9-23.6 A
# with 4-6/6. The rim centroid averages over loops that sit below the top of the
# ECD, so a groove-only ligand registers a non-zero depth; these numbers absorb
# that offset rather than assuming a zero point.
DEEP_DEPTH = 16.0
DEEP_FEATURES = 3
SHELF_DEPTH = 10.0


def assign_tier(deepest: float, split: dict[str, int], satisfied: int) -> str:
    """Which of the three portfolio tiers this geometry corresponds to."""
    if split["intracellular"] > split["tmd"] and split["intracellular"] > split["ecd"]:
        return "wrong_face"
    if deepest >= DEEP_DEPTH and satisfied >= DEEP_FEATURES:
        return "deep"
    if deepest >= SHELF_DEPTH:
        return "shelf"
    return "groove"


@app.command()
def main(
    structures: Annotated[list[Path], typer.Argument(help="complexes to score")],
    receptor_chain: Annotated[str, typer.Option(help="receptor chain id")] = "R",
    ligand_chain: Annotated[str, typer.Option(help="ligand/design chain id")] = "P",
    offset: Annotated[int, typer.Option(help="UniProt number of the receptor chain's first residue; 0 if already UniProt-numbered")] = 0,
    out: Annotated[Path | None, typer.Option(help="write a CSV here as well")] = None,
) -> None:
    rows: list[dict[str, object]] = []
    for path in structures:
        model = gemmi.read_structure(str(path))[0]
        try:
            row = score_complex(model, receptor_chain, ligand_chain, offset)
        except typer.BadParameter as exc:
            typer.echo(f"{path.stem:22s} SKIPPED: {exc}")
            continue
        row = {"structure": path.stem} | row
        rows.append(row)
        feature_text = " ".join(
            f"{f.name.split('_')[0]}={row[f.name] or '-':>5}" for f in FEATURES
        )
        typer.echo(
            f"{path.stem:22s} len={row['ligand_len']:3d} depth={row['depth_max']:6.2f}"
            f" ecd/tmd/ic={row['contacts_ecd']:2d}/{row['contacts_tmd']:2d}/{row['contacts_intracellular']:2d}"
            f" met={row['features_met']}/6 {row['tier']:10s} {feature_text}"
        )

    if out and rows:
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        typer.echo(f"wrote {out}")


if __name__ == "__main__":
    app()
