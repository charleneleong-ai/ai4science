"""Engineer histidine combinations into the de novo designs aimed at a conserved epitope.

Why not the 3DS18 chassis for these slots: 3DS18's footprint on TfR1 is 44% identical in mouse
(8 of 18 contact residues), and the substitutions land on the contacts that carry the interface --
Lys371, which 3DS18's most buried position packs against, is a glutamine there, and Arg208 an
asparagine. A 3DS18-derived binder therefore cannot satisfy the cross-reactivity objective the
brief ranks second, whatever its switch does. The designs here target the protease-like patch the
target survey measured at 90% identity between the two species.

Why every candidate has to be engineered rather than selected: of the 300 designs the generator
produced across six specs, none carries three histidines and only four carry two, so no design
has a multi-site switch to begin with.

The substitution route, the imidazole gate, the wild-type assertion and the replication
requirement are all as in `his_variants`, which this imports rather than restates. The one
difference is where the replicates come from: a de novo design replicates over the structures its
own generator wrote, and a structure is accepted here only if one of its chains carries exactly
the recorded sequence *and* its side chains were actually placed -- the generator writes an
inverse-folded structure whose side chains carry the designed sequence with every atom at the
origin, which passes a completeness check on atom names and is scored as a side chain collapsed
onto one point.

Which positions: the most buried mutable interface positions, by the burial criterion the TNF
coupling cohort validated -- rho +0.74 (propka) and +0.51 (pKAI) of the bound-versus-free pKa
difference against distance to the target, over 114 readings, against no consistent effect for
proximity to a receptor carboxylate.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Any

import numpy as np
import typer
from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
from csvio import read_designs
from his_variants import (
    Campaign,
    ChainPair,
    Site,
    Source,
    Variant,
    combinations,
    one_letter,
    report_variants,
    write_outputs,
)
from target_sites import MAX_ASA
from tnf_binder.rosetta import HisMutator

# Glycine and proline carry the backbone, cysteine would add an unpaired thiol, and a histidine
# already there needs no substituting. Aspartate and glutamate are also left alone: glutamate
# fraction is the only expression predictor that replicated on the public corpus (AUROC 0.724,
# better than any structural score tested), these designs already sit at a median 0.130 against
# the 0.283 of the reference chassis, and spending acidic residues on histidines would push the
# one number that predicts whether the protein can be made at all further the wrong way.
UNTOUCHABLE = frozenset("GPCHDE")

STRUCTURE_SUFFIXES = (".cif", ".cif.gz")

app = typer.Typer(add_completion=False)


@dataclass(frozen=True)
class Burial:
    """How much of one binder residue the interface buries, and how close it sits."""

    number: int
    one_letter: str
    dsasa: float
    rsa_free: float
    distance: float

    @property
    def mutable(self) -> bool:
        return self.one_letter not in UNTOUCHABLE

    @property
    def site(self) -> Site:
        return Site(self.one_letter, self.number)


def chain_sequence(chain: Any) -> str:
    return "".join(one_letter(r) for r in chain if is_aa(r, standard=True))


def placeholder_atoms(model: Any) -> int:
    """Atoms parked at the origin, which a completeness check on atom names cannot see.

    The Biopython-side twin of `target_sites.placeholder_coordinates`, kept here because this
    script runs in the PyRosetta environment where gemmi is absent. Both exist because the
    failure is silent in opposite directions: the verifier must refuse to measure such a
    structure, and the generator must refuse to build on one -- PyRosetta crashes on the
    zero-length bond vector, so without this the whole run dies on one bad input file.
    """
    return sum(
        1 for chain in model for residue in chain for atom in residue
        if all(abs(value) < 1e-6 for value in atom.coord)
    )


def structure_paths(root: Path, name: str) -> list[Path]:
    """Every structure file belonging to this design, matched on the file stem.

    The direction is the opposite of `ph_switch.index_structures`, which tests the design name
    inside the file stem: the pooled table prefixes each design with the spec that produced it,
    so a design recorded as `conserved_s_conserved_s_45` is written to disk as
    `conserved_s_45.cif` and it is the *stem* that is contained here. Delimiting by underscores
    stops `conserved_s_4` claiming `conserved_s_45`, and the longest matching stem wins so the
    result does not depend on directory order. Compressed structures count, as they do there.
    """
    padded = f"_{name}_"
    by_stem: dict[str, list[Path]] = {}
    for path in sorted(root.rglob("*.cif*")):
        if not path.name.endswith(STRUCTURE_SUFFIXES):
            continue
        stem = path.name.removesuffix(".gz").removesuffix(".cif")
        if f"_{stem}_" in padded:
            by_stem.setdefault(stem, []).append(path)
    return by_stem[max(by_stem, key=len)] if by_stem else []


def sources_for(root: Path, name: str, sequence: str, work: Path) -> list[Source]:
    """Every structure of this design that is a measurable complex carrying its sequence."""
    sources = []
    parser, writer = MMCIFParser(QUIET=True), PDBIO()
    for path in structure_paths(root, name):
        model = parser.get_structure(name, str(path))[0]
        sequences = {chain.id: chain_sequence(chain) for chain in model}
        populated = [chain_id for chain_id, text in sequences.items() if text]
        matches = [chain_id for chain_id in populated if sequences[chain_id] == sequence]
        if len(populated) < 2 or len(matches) != 1 or placeholder_atoms(model):
            continue
        label = path.parent.name
        destination = work / f"src_{name}_{label}.pdb"
        writer.set_structure(model.get_parent())
        writer.save(str(destination), ChainPair(populated))
        sources.append(Source(label, destination, matches[0]))
    return sources


def burials(pdb: Path, binder: str) -> list[Burial]:
    """Per-residue burial on binding, from the complex and then the binder alone."""
    sasa = ShrakeRupley()
    model = PDBParser(QUIET=True).get_structure("complex", str(pdb))[0]
    sasa.compute(model, level="R")

    residues = [r for r in model[binder] if is_aa(r, standard=True)]
    bound = {r.id[1]: float(r.sasa) for r in residues}

    # One vectorised pass per residue rather than ~440k Atom.__sub__ calls.
    target = np.array([atom.coord for chain in model if chain.id != binder
                       for residue in chain for atom in residue])
    distances = {
        r.id[1]: float(np.linalg.norm(
            np.array([atom.coord for atom in r])[:, None, :] - target[None, :, :], axis=-1).min())
        for r in residues
    }

    # The bound numbers are already copied out, so the free state can reuse this parse.
    for chain_id in [chain.id for chain in model if chain.id != binder]:
        model.detach_child(chain_id)
    sasa.compute(model, level="R")

    return [
        Burial(
            number=r.id[1],
            one_letter=one_letter(r),
            dsasa=float(r.sasa) - bound[r.id[1]],
            rsa_free=float(r.sasa) / MAX_ASA.get(one_letter(r), 200.0),
            distance=distances[r.id[1]],
        )
        for r in model[binder] if is_aa(r, standard=True)
    ]


@dataclass(frozen=True)
class Floors:
    """What a position must clear to be substituted, recorded beside the positions chosen."""

    count: int
    min_dsasa: float
    min_rsa: float

    def as_dict(self) -> dict[str, object]:
        return {"positions_per_design": self.count, "min_dsasa": self.min_dsasa,
                "min_rsa_free": self.min_rsa}


def pick_sites(measured: Sequence[Burial], floors: Floors) -> list[Burial]:
    eligible = [
        b for b in measured
        if b.mutable and b.dsasa >= floors.min_dsasa and b.rsa_free >= floors.min_rsa
    ]
    return sorted(eligible, key=lambda b: -b.dsasa)[:floors.count]


@dataclass
class Chassis:
    """One de novo design and everything decided about it before any substitution is made."""

    name: str
    sources: list[Source]
    chosen: list[Burial]

    @property
    def note(self) -> str:
        if not self.sources:
            return "no structure carries this sequence with real side-chain coordinates"
        if len(self.chosen) < 2:
            return f"only {len(self.chosen)} position(s) clear the burial floor"
        return ""

    @property
    def usable(self) -> bool:
        return not self.note

    def summary(self) -> dict[str, object]:
        return {
            "design": self.name,
            "structures": [source.label for source in self.sources],
            "positions": [
                {"site": b.site.label, "dsasa": round(b.dsasa, 1),
                 "rsa_free": round(b.rsa_free, 2), "distance": round(b.distance, 1)}
                for b in self.chosen
            ],
            "usable": self.usable,
            "note": self.note,
        }


def prepare(name: str, sequence: str, root: Path, work: Path, floors: Floors) -> Chassis:
    sources = sources_for(root, name, sequence, work)
    chosen = pick_sites(burials(sources[0].path, sources[0].binder), floors) if sources else []
    return Chassis(name, sources, chosen)


def build_variants(entry: Chassis, mutator: HisMutator, out_dir: Path,
                   repacks: int) -> list[Variant]:
    campaign = Campaign(mutator, entry.sources, out_dir, stem=entry.name, prefix="dn",
                        repacks=repacks)
    return [campaign.build(c) for c in combinations([b.site for b in entry.chosen])]


def report_positions(chassis: Sequence[Chassis]) -> None:
    typer.echo("\nper design, the positions chosen on burial:")
    for entry in chassis:
        if not entry.usable:
            typer.echo(f"  {entry.name:26s} SKIPPED - {entry.note}")
            continue
        sites = "  ".join(f"{b.site.label}({b.dsasa:.0f}A2,{b.distance:.1f}A)"
                          for b in entry.chosen)
        typer.echo(f"  {entry.name:26s} {len(entry.sources)} structures  {sites}")


@app.command()
def main(
    designs: Annotated[Path, typer.Option(help="pooled CSV from collect_designs.py")],
    structures: Annotated[Path, typer.Option(help="BoltzGen output root")],
    out_dir: Annotated[Path, typer.Option(help="where the variant complexes are written")],
    out_designs: Annotated[Path, typer.Option(help="name,sequence CSV for ph_switch.py")],
    only: Annotated[str, typer.Option(help="comma-separated design names to engineer")] = "",
    n_sites: Annotated[int, typer.Option(help="positions to substitute per design")] = 5,
    min_dsasa: Annotated[float, typer.Option(help="minimum area buried on binding")] = 20.0,
    min_rsa: Annotated[float, typer.Option(help="minimum exposure in the free binder")] = 0.20,
    repacks: Annotated[int, typer.Option(
        help="repacks per structure; >1 measures rotamer sensitivity and must not be ranked on, "
             "because repacks of one pose are correlated evidence the verifier cannot "
             "distinguish from independent structures")] = 1,
    repack_radius: Annotated[float, typer.Option(
        help="angstroms around each substitution that are repacked")] = 6.0,
) -> None:
    wanted = {name.strip() for name in only.split(",") if name.strip()}
    sequences = read_designs(designs, wanted)
    if missing := wanted - set(sequences):
        raise typer.BadParameter(f"not in {designs}: {', '.join(sorted(missing))}")
    if not sequences:
        raise typer.BadParameter(f"no designs selected from {designs}")

    floors = Floors(n_sites, min_dsasa, min_rsa)
    out_dir.mkdir(parents=True, exist_ok=True)
    prepared = [prepare(name, sequence, structures, out_dir, floors)
                for name, sequence in sequences.items()]

    mutator = HisMutator(repack_radius)
    variants: list[Variant] = []
    for entry in prepared:
        if not entry.usable:
            continue
        built = build_variants(entry, mutator, out_dir, repacks)
        variants.extend(built)
        typer.echo(f"  {entry.name:26s} {len(built)} variants built")

    report_positions(prepared)
    sidecar = out_designs.with_suffix(".chassis.json")
    write_outputs(out_designs, variants, sidecar,
                  [{"floors": floors.as_dict(), "repacks": repacks},
                   *(entry.summary() for entry in prepared)])
    report_variants(variants, width=34)


if __name__ == "__main__":
    app()
