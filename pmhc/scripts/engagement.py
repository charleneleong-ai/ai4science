"""Does the design touch the peptide, or just the groove it sits in?

Challenge 5's second ranking tier is peptide-dependence: a binder that grips the HLA helices and
ignores the peptide binds the peptide-free receptor just as well, and the peptide is only
**1.2%** of the complex's accessible surface, so that is the easy thing for a generator to do.
A three-way co-fold is the real test, but it costs two extra predictions per design across two
predictors. This reads the generator's own pose instead, which is free, and narrows the field
before any of that is spent.

The bar is not invented. The T-cell receptor in 7OW6 discriminates these two peptides, so its
own contact split is what a peptide-reading binder looks like on this target; `--control`
measures it with the same code. Thresholds are expressed relative to that measurement rather
than chosen, because a share that sounds low in the abstract may be all this epitope affords.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import typer
from common import binder_chain, chain_sequence, heavy_coords, load, target_sequences, write_rows

app = typer.Typer(add_completion=False)

# Heavy-atom contact distance. 4.5 A is what tfr1 and glp1r call a contact, and this follows
# them so the three campaigns' numbers can be compared; epitope_geometry.py quotes the receptor
# readout at an explicit 4.0 and 5.5 A instead, because that analysis is about which atom class
# a known receptor touches rather than about how much interface there is.
CONTACT = 4.5

# Which peptide positions each generation variant nominated, from prepare_target.py.
VARIANT_HOTSPOTS = {
    "bulge_core": (5, 6, 7),
    "bulge_wide": (4, 5, 6, 7, 8),
    "mut_proximal": (5, 6),
    "hub": (7, 8),
}

MUTATION_POSITION = 6


@dataclass(frozen=True)
class Engagement:
    design: str
    variant: str
    binder_length: int
    peptide_contacts: int
    hla_contacts: int
    hotspot_contacts: int
    mutation_contacts: int
    peptide_residues: int
    hotspots_touched: int
    hotspots_nominated: int

    @property
    def peptide_share(self) -> float:
        """Peptide contacts as a fraction of all target contacts.

        The discriminating quantity: a design reading the groove scores near zero here however
        many contacts it makes in total, and a raw peptide-contact count cannot tell the two
        apart because a larger interface raises both.
        """
        total = self.peptide_contacts + self.hla_contacts
        return self.peptide_contacts / total if total else 0.0

    @property
    def row(self) -> dict[str, object]:
        return {
            "design": self.design,
            "variant": self.variant,
            "binder_length": self.binder_length,
            "peptide_contacts": self.peptide_contacts,
            "hla_contacts": self.hla_contacts,
            "peptide_share": round(self.peptide_share, 4),
            "hotspot_contacts": self.hotspot_contacts,
            "mutation_contacts": self.mutation_contacts,
            "peptide_residues_touched": self.peptide_residues,
            "hotspots_touched": self.hotspots_touched,
            "hotspots_nominated": self.hotspots_nominated,
        }


def contact_pairs(
    binder: np.ndarray, target: np.ndarray, cutoff: float = CONTACT
) -> np.ndarray:
    """Boolean matrix of binder-atom x target-atom pairs within `cutoff`."""
    if not len(binder) or not len(target):
        return np.zeros((len(binder), len(target)), dtype=bool)
    gaps = binder[:, None, :] - target[None, :, :]
    return (gaps**2).sum(axis=2) <= cutoff**2


def variant_of(design: str) -> str:
    for name in VARIANT_HOTSPOTS:
        if design.startswith(name):
            return name
    raise typer.BadParameter(f"{design}: no known variant prefix")


def peptide_chain_name(structure, epitopes: tuple[str, ...]) -> str:
    for chain in structure[0]:
        if chain_sequence(chain) in epitopes:
            return chain.name
    raise typer.BadParameter("no chain matches either epitope sequence")


def stack(model, names: list[str]) -> np.ndarray:
    parts = [heavy_coords(model[name])[0] for name in names]
    return np.concatenate(parts) if parts else np.empty((0, 3))


def measure_pose(
    model,
    label: str,
    variant: str,
    binder_names: list[str],
    peptide_name: str,
    hotspots: tuple[int, ...],
) -> Engagement:
    """Contact split for one pose, given which chains play which role.

    The roles are the caller's to decide: a design is the single chain that is not the target's,
    while the 7OW6 control's binder is two chains, the receptor alpha and beta.
    """
    binder = stack(model, binder_names)
    peptide, peptide_positions = heavy_coords(model[peptide_name])
    hla = stack(
        model,
        [c.name for c in model if c.name not in set(binder_names) | {peptide_name}],
    )

    touching_peptide = contact_pairs(binder, peptide)
    positions = np.array(peptide_positions)
    touched = {int(p) for p in positions[touching_peptide.any(axis=0)]}

    return Engagement(
        design=label,
        variant=variant,
        binder_length=sum(len(model[name]) for name in binder_names),
        peptide_contacts=int(touching_peptide.sum()),
        hla_contacts=int(contact_pairs(binder, hla).sum()),
        hotspot_contacts=int(touching_peptide[:, np.isin(positions, hotspots)].sum()),
        mutation_contacts=int(touching_peptide[:, positions == MUTATION_POSITION].sum()),
        peptide_residues=len(touched),
        hotspots_touched=len(touched & set(hotspots)),
        hotspots_nominated=len(hotspots),
    )


def measure(
    pose: Path, targets: frozenset[str], epitopes: tuple[str, ...], hotspots: tuple[int, ...]
) -> Engagement:
    structure = load(pose)
    return measure_pose(
        structure[0],
        pose.stem,
        variant_of(pose.stem),
        [binder_chain(pose, structure, targets)],
        peptide_chain_name(structure, epitopes),
        hotspots,
    )


def summarise(label: str, rows: list[Engagement]) -> None:
    shares = sorted(r.peptide_share for r in rows)
    if not shares:
        return
    typer.echo(
        f"{label}: n={len(shares)} peptide share"
        f" min={shares[0]:.3f} median={shares[len(shares) // 2]:.3f} max={shares[-1]:.3f}"
    )


@app.command()
def main(
    poses: Path = typer.Argument(..., help="directory of generator poses (.cif)"),
    out: Path = typer.Option(Path("engagement.csv")),
    references: Path = typer.Option(Path("targets/novelty_references.json")),
    control: Path = typer.Option(
        None, help="a deposition holding a known binder, e.g. the carved-free 7OW6"
    ),
    control_binder: str = typer.Option("D,E", help="its binder chains; 7OW6's TCR is D and E"),
    control_hotspots: str = typer.Option("5,6,7,8", help="positions to score the control on"),
) -> None:
    targets = target_sequences(references)
    epitopes = tuple(
        s for s in targets if len(s) == 10 and s.startswith("VVVGA")
    )
    if len(epitopes) != 2:
        raise typer.BadParameter(f"expected both 10-mer epitopes in {references}, got {epitopes}")

    if control is not None:
        structure = load(control)
        measured = measure_pose(
            structure[0],
            control.stem,
            "control",
            [c.strip() for c in control_binder.split(",")],
            peptide_chain_name(structure, epitopes),
            tuple(int(p) for p in control_hotspots.split(",")),
        )
        for key, value in measured.row.items():
            typer.echo(f"  {key:26s} {value}")
        typer.echo(
            f"\nthe known binder's peptide share is {measured.peptide_share:.3f}"
            f" over {measured.peptide_contacts + measured.hla_contacts} contacts"
        )
        typer.echo(f"wrote {out} ({write_rows([measured], out)} rows)")
        return

    found = sorted(p for p in poses.rglob("*.cif") if p.is_file())
    if not found:
        raise typer.BadParameter(f"{poses} holds no .cif poses")

    rows = [
        measure(pose, targets, epitopes, VARIANT_HOTSPOTS[variant_of(pose.stem)])
        for pose in found
    ]

    for row in sorted(rows, key=lambda r: -r.peptide_share)[:15]:
        typer.echo(
            f"  {row.design:20s} share={row.peptide_share:.3f}"
            f" peptide={row.peptide_contacts:4d} hla={row.hla_contacts:5d}"
            f" hotspots={row.hotspots_touched}/{row.hotspots_nominated}"
            f" P6={row.mutation_contacts:3d}"
        )
    summarise("all designs", rows)
    for variant in VARIANT_HOTSPOTS:
        summarise(f"  {variant}", [r for r in rows if r.variant == variant])
    typer.echo(f"wrote {out} ({write_rows(rows, out)} rows)")


if __name__ == "__main__":
    app()
