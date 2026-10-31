"""Engineer histidine combinations into the 3DS18 chassis and write them as scorable complexes.

Why mutate rather than co-fold: 3DS18's pose on TfR1 is crystallographic (6WRV). Re-predicting
the complex would discard the one property that makes this chassis worth using over the de novo
designs, whose best interface confidence only reaches the bottom of the range the same pipeline
produced on an easier target. The backbone is held fixed here and only the substituted side chain
plus its neighbourhood is repacked, so every variant is measured in the validated pose.

Why PyRosetta and not a renamed residue: propka cannot score an incomplete residue. Renaming a
residue to HIS and keeping backbone plus CB leaves no imidazole, propka silently omits it, and
every pKa comes back missing -- so a variant whose switch was never measured would be
indistinguishable from one that has no switch. Every structure written here is gated on its
imidazole being complete, and a variant that loses one is reported rather than scored.

Why three structures per variant: 6WRV's asymmetric unit holds three independent copies of the
complex (binder C on target B, D on A, F on E). The verifier returns a verdict only where at
least two structures agree, and crystallographically independent copies are a real replication
axis rather than three draws from one model's ensemble. Chain F also contacts A -- a packing
neighbour, excluded so that each variant is scored against exactly one receptor.

Which positions: the five most buried mutable interface positions. The TNF coupling cohort
measured burial as the lever that drives the bound-versus-free pKa difference negative -- rho
+0.74 (propka) and +0.51 (pKAI) against distance to the target over 114 readings, buried sites
sitting at a median delta of -0.98 against -0.04 for the rest -- while proximity to a receptor
carboxylate showed no consistent effect and the two methods disagreed even on its sign. Arg84 and
the interface acids rank high on burial and are still left alone: they carry the interface's
charge pairs, and the brief ranks affinity third rather than nowhere.
"""

from __future__ import annotations

import itertools
import json
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Any

import typer
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa
from csvio import write_rows
from tnf_binder.rosetta import HisMutator

# The ring propka needs to see; backbone plus CB is not a histidine to it.
IMIDAZOLE = frozenset({"CG", "ND1", "CD2", "CE1", "NE2"})

# Each binder chain in 6WRV's asymmetric unit, with the receptor copy it actually binds.
COPIES: tuple[tuple[str, str], ...] = (("C", "B"), ("D", "A"), ("F", "E"))

DEFAULT_POSITIONS = "Y87,Q26,L33,V30,N37"
DESIGN_COLUMNS = ["name", "sequence", "substitutions", "structures"]

app = typer.Typer(add_completion=False)


def one_letter(residue: Any) -> str:
    """One-letter code for a Bio.PDB residue, uppercase, or `X` if it is not standard."""
    return protein_letters_3to1.get(residue.get_resname().capitalize(), "X")


@dataclass(frozen=True)
class Site:
    """One wild-type interface position, as written on the command line."""

    wild_type: str
    number: int

    @property
    def label(self) -> str:
        return f"{self.wild_type}{self.number}H"


@dataclass(frozen=True)
class Source:
    """One structure a variant can be built on, and which of its chains is the binder.

    The replication axis differs by chassis -- 3DS18 replicates over the independent copies in
    one crystal, a de novo design over the structures its own generator wrote -- so where a
    `Source` comes from is the producer's business and not encoded in the type.
    """

    label: str
    path: Path
    binder: str


@dataclass
class Variant:
    """What one histidine combination produced, including the structures that refused."""

    name: str
    sites: tuple[Site, ...]
    sequence: str = ""
    structures: list[str] = field(default_factory=list)
    failures: list[str] = field(default_factory=list)

    @property
    def labels(self) -> tuple[str, ...]:
        return tuple(site.label for site in self.sites)

    @property
    def verdict(self) -> str:
        if len(self.structures) >= 2:
            return "ok"
        return "unreplicated" if self.structures else "failed"

    def row(self) -> dict[str, object]:
        return {
            "name": self.name,
            "sequence": self.sequence,
            "substitutions": " ".join(self.labels),
            "structures": len(self.structures),
        }

    def summary(self) -> dict[str, object]:
        return self.row() | {
            "substitutions": list(self.labels),
            "structures": self.structures,
            "failures": self.failures,
            "verdict": self.verdict,
        }


class ChainPair(Select):
    """The named chains, standard residues, one altloc -- what Rosetta reads without guessing."""

    def __init__(self, chains: Sequence[str]) -> None:
        self.chains = frozenset(chains)

    def accept_chain(self, chain: Any) -> bool:
        return chain.id in self.chains

    def accept_residue(self, residue: Any) -> bool:
        return residue.id[0] == " " and is_aa(residue, standard=True)

    def accept_atom(self, atom: Any) -> bool:
        return atom.get_altloc() in (" ", "A")


def parse_sites(spec: str) -> list[Site]:
    sites = []
    for token in (part.strip() for part in spec.split(",") if part.strip()):
        if len(token) < 2 or not token[1:].isdigit():
            raise typer.BadParameter(f"{token!r} is not a position like Y87")
        sites.append(Site(token[0].upper(), int(token[1:])))
    return sites


def extract_copies(cif: Path, work: Path, sites: Sequence[Site]) -> list[Source]:
    """Check the named positions against the crystal, then split it into one PDB per copy.

    The check and the split are one function so the gate cannot be bypassed by a caller that
    only wants the structures, and so the 3.5 MB mmCIF is parsed once rather than twice.
    """
    structure = MMCIFParser(QUIET=True).get_structure("6wrv", str(cif))
    model = structure[0]
    wrong = [
        f"chain {binder} {site.number} is {one_letter(model[binder][site.number])}, "
        f"not {site.wild_type}"
        for binder, _ in COPIES
        for site in sites
        if one_letter(model[binder][site.number]) != site.wild_type
    ]
    if wrong:
        raise typer.BadParameter(
            "the structure disagrees with the positions asked for: " + "; ".join(wrong)
            + ". Mutating a position that is not what it was chosen for would move a different "
              "side chain than the burial ranking measured."
        )

    writer = PDBIO()
    writer.set_structure(structure)
    sources = []
    for binder, target in COPIES:
        label = f"{binder}{target}"
        destination = work / f"copy_{label}.pdb"
        writer.save(str(destination), ChainPair((binder, target)))
        sources.append(Source(label, destination, binder))
    return sources


def finalise(pdb: Path, destination: Path, chain: str, numbers: Sequence[int]) -> list[str]:
    """Convert to CIF, but only once every substituted residue carries a full imidazole."""
    structure = PDBParser(QUIET=True).get_structure(destination.stem, str(pdb))
    model = structure[0]
    incomplete = [
        f"{chain}{number}"
        for number in numbers
        if model[chain][number].get_resname() != "HIS"
        or not IMIDAZOLE <= {atom.get_name() for atom in model[chain][number]}
    ]
    if not incomplete:
        writer = MMCIFIO()
        writer.set_structure(structure)
        writer.save(str(destination))
    return incomplete


def combinations(sites: Sequence[Site]) -> list[tuple[Site, ...]]:
    """Every subset including the wild type, so the chassis is measured on the same footing."""
    return [
        tuple(combination)
        for size in range(len(sites) + 1)
        for combination in itertools.combinations(sites, size)
    ]


def design_name(stem: str, combination: Sequence[Site]) -> str:
    ordered = sorted(combination, key=lambda site: site.number)
    return stem + ("".join(f"-{site.label}" for site in ordered) if ordered else "-wt")


@dataclass
class Campaign:
    """One PyRosetta session against one set of source complexes."""

    mutator: HisMutator
    sources: Sequence[Source]
    out_dir: Path
    stem: str = "3ds18"
    prefix: str = "var"
    repacks: int = 1
    """Repacks of each source. The packer is stochastic, so repeated runs sample the rotamers a
    substitution could adopt -- which measures how much a verdict depends on the rotamer chosen,
    *not* how much it depends on the structure. Those are correlated draws from one pose: they
    share its backbone, its neighbourhood and its burial, so they agree with each other far more
    readily than independent structures do. The verifier weighs all structures equally, so
    feeding it repacks would convert `unreplicated` into a verdict on correlated evidence,
    preferentially for the thinnest chassis. Left at 1 for that reason; raise it only to measure
    rotamer sensitivity, and keep that measurement out of the pass that ranks designs."""

    def build(self, combination: Sequence[Site]) -> Variant:
        variant = Variant(design_name(self.stem, combination), tuple(combination))
        # Repacking an unsubstituted pose would write the same coordinates twice and the verifier
        # would read the duplicate as agreement.
        attempts = self.repacks if combination else 1
        for source in self.sources:
            for attempt in range(attempts):
                label = source.label if attempts == 1 else f"{source.label}p{attempt + 1}"
                self.add(variant, source, label)
        return variant

    def add(self, variant: Variant, source: Source, label: str) -> None:
        pose = self.mutator.load(source.path)
        if wrong := misplaced(pose, source.binder, variant.sites):
            variant.failures.append(f"{label}: {'; '.join(wrong)}")
            return

        numbers = [site.number for site in variant.sites]
        sites = [(source.binder, number) for number in numbers]
        if sites and not self.mutator.mutate(pose, sites):
            variant.failures.append(f"{label}: a substituted position is absent from the pose")
            return

        scratch = self.out_dir / f"pose_{variant.name}_{label}.pdb"
        destination = self.out_dir / f"{self.prefix}_{variant.name}_{label}.cif"
        self.mutator.write(pose, scratch)
        try:
            if incomplete := finalise(scratch, destination, source.binder, numbers):
                variant.failures.append(f"{label}: imidazole incomplete at {','.join(incomplete)}")
                return
        finally:
            scratch.unlink(missing_ok=True)

        variant.structures.append(destination.name)
        variant.sequence = variant.sequence or self.mutator.chain_sequence(pose, source.binder)


def misplaced(pose: Any, chain: str, sites: Sequence[Site]) -> list[str]:
    """Positions that do not hold the wild-type residue they were chosen for, in this pose.

    Checked per structure rather than once per chassis, because structures of one design are
    matched on *sequence* and need not agree on numbering. A de novo design's positions are
    picked from one of its structures and then substituted in all of them, so a template offset
    in another structure would silently mutate a different side chain -- and the imidazole gate
    could not tell, because the residue it finds there really is a complete histidine.
    """
    info = pose.pdb_info()
    wrong = []
    for site in sites:
        index = info.pdb2pose(chain, site.number)
        if index == 0:
            wrong.append(f"{site.label} is absent from the pose")
        elif (found := pose.residue(index).name1()) != site.wild_type:
            wrong.append(f"{site.label} found {found}, not {site.wild_type}")
    return wrong


def report_variants(variants: Sequence[Variant], width: int = 28) -> None:
    scorable = [variant for variant in variants if variant.verdict == "ok"]
    typer.echo(f"\n{len(variants)} variants, {sum(len(v.structures) for v in variants)} "
               f"structures written, {len(scorable)} with the two structures the verifier needs "
               f"before it will return a verdict")
    for variant in variants:
        if variant.verdict != "ok":
            typer.echo(f"  {variant.name:{width}s} {variant.verdict}: "
                       f"{'; '.join(variant.failures)}")


def write_outputs(out_designs: Path, variants: Sequence[Variant],
                  sidecar: Path, payload: Sequence[dict[str, object]]) -> None:
    """The designs table plus its provenance sidecar; only variants with a structure are scored."""
    write_rows(out_designs, [v.row() for v in variants if v.structures], DESIGN_COLUMNS)
    sidecar.write_text(json.dumps(payload, indent=2) + "\n")
    typer.echo(f"wrote {out_designs} and {sidecar}")


@app.command()
def main(
    cif: Annotated[Path, typer.Option(help="6WRV mmCIF, the 3DS18-TfR1 complex")],
    out_dir: Annotated[Path, typer.Option(help="where the variant complexes are written")],
    out_designs: Annotated[Path, typer.Option(help="name,sequence CSV for ph_switch.py")],
    positions: Annotated[str, typer.Option(
        help="comma-separated wild-type positions to substitute")] = DEFAULT_POSITIONS,
    repack_radius: Annotated[float, typer.Option(
        help="angstroms around each substitution that are repacked")] = 6.0,
) -> None:
    sites = parse_sites(positions)
    if not sites:
        raise typer.BadParameter("no positions given")

    out_dir.mkdir(parents=True, exist_ok=True)
    sources = extract_copies(cif, out_dir, sites)
    typer.echo(f"extracted {len(sources)} independent complexes: "
               f"{', '.join(f'{b} on {t}' for b, t in COPIES)}")

    campaign = Campaign(HisMutator(repack_radius), sources, out_dir)
    variants = []
    for combination in combinations(sites):
        variant = campaign.build(combination)
        variants.append(variant)
        typer.echo(f"  {variant.name:28s} {len(variant.structures)}/{len(sources)} structures")

    typer.echo(f"\n{len(variants)} variants over {len(sites)} positions "
               f"({', '.join(site.label for site in sites)})")
    write_outputs(out_designs, variants, out_designs.with_suffix(".variants.json"),
                  [variant.summary() for variant in variants])
    report_variants(variants)


if __name__ == "__main__":
    app()
