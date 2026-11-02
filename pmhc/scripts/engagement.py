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

A caveat the sibling campaign paid for: tnf-binder abandoned raw contact counts because the
count scales with interpenetration, so clashing poly-alanine designs out-scored real ones. The
ratio reported here is immune to the scale half of that, since a bigger interface raises both
sides, but not to interpenetration itself -- a pose that clashes into the groove inflates the
denominator. The co-fold screen, not this gate, is what separates a plausible complex from a
geometric coincidence.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from statistics import median

import gemmi
import numpy as np
import typer
from common import (
    VARIANTS,
    ChainRoles,
    chain_roles,
    chain_sequences,
    epitopes_from,
    heavy_coords,
    load,
    pose_files,
    target_references,
    write_rows,
    wrote,
)

app = typer.Typer(add_completion=False)

# Heavy-atom contact distance. 4.5 A is what tfr1 and glp1r call a contact, and this follows
# them so the three campaigns' numbers can be compared; epitope_geometry.py quotes the receptor
# readout at an explicit 4.0 and 5.5 A instead, because that analysis is about which atom class
# a known receptor touches rather than about how much interface there is. The value travels in
# the output as `contact_cutoff`: the calibrated 0.200 bar is only meaningful against shares
# measured at the same cutoff, and a re-measure at 4.0 A would otherwise be silently comparable.
CONTACT = 4.5

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
    contact_cutoff: float = CONTACT

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
            "contact_cutoff": self.contact_cutoff,
        }


def contact_pairs(
    binder: np.ndarray, target: np.ndarray, cutoff: float = CONTACT
) -> np.ndarray:
    """Boolean matrix of binder-atom x target-atom pairs within `cutoff`.

    Squared distances via `(a-b)^2 = a^2 + b^2 - 2ab`, which reaches the (n, m) answer through
    one BLAS call instead of materialising the (n, m, 3) difference. On this cohort's real
    shapes -- 672 binder atoms against 3,079 receptor atoms -- that is 19.0 ms down to 1.5 ms
    and 116 MB down to 33 MB per pose, 4.5 s and 82 MB over 256 poses, with every row of the
    output identical.

    float64 throughout, deliberately. The identity loses precision near the cutoff, and in
    float32 that flips 52k-180k of 69.6M pairs placed within 1e-7 A of the boundary, which is
    enough to move a gate value; in float64 it flips none at coordinates in any deposition's
    range.
    """
    if not len(binder) or not len(target):
        return np.zeros((len(binder), len(target)), dtype=bool)
    binder = np.asarray(binder, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    squared = (binder * binder).sum(axis=1)[:, None] + (target * target).sum(axis=1)[None, :]
    squared -= 2.0 * (binder @ target.T)
    return squared <= cutoff**2


def variant_of(design: str) -> str:
    """The generation variant a design came from, by name prefix.

    BoltzGen names designs after the spec stem and prepare_target names each spec after its
    variant, so the prefix is reliable for this cohort -- but the coupling is implicit. The
    variant is also authoritative in the pose path (`gen_out/<variant>/...`), which is the
    sounder source and is recorded as a follow-up.
    """
    for name in VARIANTS:
        if design.startswith(name):
            return name
    raise typer.BadParameter(f"{design}: no known variant prefix")


def stack(model: gemmi.Model, names: tuple[str, ...]) -> np.ndarray:
    parts = [heavy_coords(model[name])[0] for name in names]
    return np.concatenate(parts) if parts else np.empty((0, 3))


def measure_pose(
    model: gemmi.Model,
    label: str,
    variant: str,
    roles: ChainRoles,
    hotspots: tuple[int, ...],
    cutoff: float = CONTACT,
) -> Engagement:
    """Contact split for one pose, given the partition of its chains."""
    binder = stack(model, roles.design)
    peptide, peptide_positions = heavy_coords(model[roles.peptide])
    hla = stack(model, roles.receptor)

    touching_peptide = contact_pairs(binder, peptide, cutoff)
    positions = np.array(peptide_positions)
    touched = {int(p) for p in positions[touching_peptide.any(axis=0)]}

    return Engagement(
        design=label,
        variant=variant,
        binder_length=sum(len(model[name]) for name in roles.design),
        peptide_contacts=int(touching_peptide.sum()),
        hla_contacts=int(contact_pairs(binder, hla, cutoff).sum()),
        hotspot_contacts=int(touching_peptide[:, np.isin(positions, hotspots)].sum()),
        mutation_contacts=int(touching_peptide[:, positions == MUTATION_POSITION].sum()),
        peptide_residues=len(touched),
        hotspots_touched=len(touched & set(hotspots)),
        hotspots_nominated=len(hotspots),
        contact_cutoff=cutoff,
    )


def measure(
    pose: Path, targets: frozenset[str], epitopes: tuple[str, ...], cutoff: float = CONTACT
) -> Engagement:
    structure = load(pose)
    variant = variant_of(pose.stem)
    return measure_pose(
        structure[0],
        pose.stem,
        variant,
        chain_roles(structure, epitopes, targets=targets, pose=pose),
        VARIANTS[variant],
        cutoff,
    )


def control_chains(structure: gemmi.Structure, controls: Path) -> tuple[str, ...]:
    """The control's binder chains, resolved by sequence against prepare_target's controls.json.

    Resolved rather than named on the command line so that the control goes through the same
    sequence-based identification as a design, which is the point of measuring it at all.
    """
    wanted = set(json.loads(controls.read_text())["tcr"].values())
    found = tuple(name for name, seq in chain_sequences(structure).items() if seq in wanted)
    if len(found) != len(wanted):
        raise typer.BadParameter(
            f"{controls} lists {len(wanted)} control chains; matched {found} in the structure"
        )
    return found


def summarise(label: str, rows: list[Engagement]) -> None:
    shares = sorted(r.peptide_share for r in rows)
    if not shares:
        return
    typer.echo(
        f"{label}: n={len(shares)} peptide share"
        f" min={shares[0]:.3f} median={median(shares):.3f} max={shares[-1]:.3f}"
    )


@app.command()
def main(
    poses: Path = typer.Argument(..., help="directory of generator poses (.cif)"),
    out: Path = typer.Option(Path("engagement_cohort.csv")),
    references: Path = typer.Option(Path("targets/novelty_references.json")),
    control: Path | None = typer.Option(
        None, help="a deposition holding a known binder, e.g. the carved-free 7OW6"
    ),
    control_out: Path = typer.Option(Path("engagement_control.csv")),
    controls: Path = typer.Option(Path("targets/controls.json")),
    cutoff: float = typer.Option(CONTACT, help="heavy-atom contact distance, Angstrom"),
) -> None:
    epitopes = epitopes_from(target_references(references))
    targets = frozenset(target_references(references).values())

    if control is not None:
        structure = load(control)
        # The control nominates no hotspots: the threshold read off it is the peptide share,
        # which is hotspot-independent.
        measured = measure_pose(
            structure[0],
            control.stem,
            "control",
            chain_roles(structure, epitopes, design=control_chains(structure, controls)),
            (),
            cutoff,
        )
        for key, value in measured.row.items():
            typer.echo(f"  {key:26s} {value}")
        typer.echo(
            f"\nthe known binder's peptide share is {measured.peptide_share:.3f}"
            f" over {measured.peptide_contacts + measured.hla_contacts} contacts"
            f" at {cutoff} A"
        )
        typer.echo(wrote(control_out, write_rows([measured], control_out)))
        return

    rows = [measure(pose, targets, epitopes, cutoff) for pose in pose_files(poses)]

    for row in sorted(rows, key=lambda r: -r.peptide_share)[:15]:
        typer.echo(
            f"  {row.design:20s} share={row.peptide_share:.3f}"
            f" peptide={row.peptide_contacts:4d} hla={row.hla_contacts:5d}"
            f" hotspots={row.hotspots_touched}/{row.hotspots_nominated}"
            f" P6={row.mutation_contacts:3d}"
        )
    summarise("all designs", rows)
    for variant in VARIANTS:
        summarise(f"  {variant}", [r for r in rows if r.variant == variant])
    typer.echo(wrote(out, write_rows(rows, out)))


if __name__ == "__main__":
    app()
