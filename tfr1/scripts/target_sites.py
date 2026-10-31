"""Which epitope on TfR1, measured rather than taken from the brief.

The challenge recommends the apical domain and states the problem with it in the
same breath: mouse differs there at six of nine residues. Affinity is ranked
*third*; mouse cross-reactivity is ranked second and pH-selectivity first. So the
recommended epitope optimises the objective that matters least and fights the one
ranked above it.

This script measures the trade-off instead of assuming it:

- the two epitopes anyone has actually hit on this receptor, from structures:
  parvovirus B19 VP1u (11RN) and the designed binder 3DS18 (6WRV);
- human/mouse identity at every contact, from a real alignment of P02786 against
  Q62351, so "six of nine" is checked rather than cited;
- a scan of the whole ectodomain surface for patches that are **conserved**,
  **exposed on the dimer**, and **clear of the transferrin site** (619-670) and of
  the N-glycans. The assay protein is glycosylated and the models are not, so a
  patch under a glycan looks bindable here and is not;
- the titratable context of each patch. A pH switch needs a protonatable residue
  *in* the interface: the TNF-alpha campaign placed six histidines 9-21 A from the
  target, found them uncoupled from binding on every sample, and lost its central
  claim to that. Carboxylates facing the patch are what a designed histidine can
  pair with, so they are counted here, before any design exists.

Dimer SASA is used throughout, not monomer: the receptor is supplied as a
homodimer, so the dimer interface is not available to bind and must not be scored
as surface.

**Every threshold that decides the epitope is a CLI option**, the resolved values
are written beside the results, and a sensitivity table reports how the candidate
count moves with the two that bind. An earlier revision hardcoded them, and two of
the three epitopes this campaign designs against turned out to be excluded by the
script's own gate while being carried forward anyway -- a stale-conclusion failure
with no error path. `prepare_target.py` now refuses to cut a target whose centre is
not usable in the scan it reads.
"""

from __future__ import annotations

import csv
import json
import statistics as st
from collections import Counter
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import gemmi
import typer
from Bio import Align
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

# Tien et al. 2013, theoretical maximum solvent-accessible area per residue (A^2).
MAX_ASA: dict[str, float] = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "Q": 225.0,
    "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0, "L": 201.0, "K": 236.0,
    "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0,
    "Y": 263.0, "V": 174.0,
}

HYDROGEN = gemmi.Element("H")
CONTACT_CUTOFF = 4.5

# Transferrin footprint, Lawrence et al. 1999 / Cheng et al. 2004 (1SUV). UniProt
# P02786 carries no binding-region feature, so this is the one conclusion-bearing
# constant here that cannot be parsed from the data.
TRANSFERRIN_SITE = frozenset(range(619, 671))

# N-linked sites parsed from the UniProt features; the O-linked T104 is excluded
# deliberately, being both outside the construct and not glycan-shielding here.
N_LINKED = "N-linked"

# TfR1 domain boundaries (Lawrence et al. 1999): the apical domain is an insertion
# in the protease-like domain, which is why its sequence range is discontinuous.
DOMAINS: dict[str, frozenset[int]] = {
    "protease-like": frozenset(range(121, 189)) | frozenset(range(384, 607)),
    "apical": frozenset(range(189, 384)),
    "helical": frozenset(range(607, 761)),
}

# A binder histidine needs a partner to titrate against; these are what it can pair with.
ACIDIC = frozenset("DE")
TITRATABLE = frozenset("DEH")

# The brief names six divergent apical residues. 210 and 211 differ only in where
# the gap falls inside the 210-214 indel, which is alignment-ambiguous and leaves
# both non-conserved either way; anything else differing is a real disagreement.
BRIEF_DIVERGENT: dict[int, str] = {204: "Q", 208: "N", 210: "D", 211: "P", 348: "K", 371: "Q"}
INDEL_AMBIGUOUS = frozenset({210, 211})

EPITOPE_SOURCES: tuple[tuple[str, str, frozenset[str]], ...] = (
    ("VP1u (11RN)", "apical", frozenset({"X", "Y"})),
    ("3DS18 designed (6WRV)", "designed", frozenset({"C", "D", "F"})),
)

app = typer.Typer(add_completion=False)


@dataclass(frozen=True)
class Thresholds:
    patch_radius: float
    min_patch_size: int
    surface_min_rsa: float
    transferrin_clearance: float
    glycan_clearance: float

    def as_dict(self) -> dict[str, float | int]:
        return {
            "patch_radius": self.patch_radius,
            "min_patch_size": self.min_patch_size,
            "surface_min_rsa": self.surface_min_rsa,
            "transferrin_clearance": self.transferrin_clearance,
            "glycan_clearance": self.glycan_clearance,
        }


@dataclass(frozen=True)
class Residue:
    number: int
    amino_acid: str
    rsa: float
    centroid: gemmi.Position


@dataclass(frozen=True)
class Patch:
    centre: int
    members: tuple[int, ...]
    acidic_members: tuple[int, ...]
    conserved: int
    counted: int
    mean_rsa: float
    transferrin_distance: float
    glycan_distance: float

    @property
    def divergent(self) -> int:
        return self.counted - self.conserved

    @property
    def conservation(self) -> float:
        return self.conserved / self.counted if self.counted else 0.0

    @property
    def domain(self) -> str:
        return domain_of(self.centre)

    def rejections(self, limits: Thresholds) -> tuple[str, ...]:
        reasons = []
        if len(self.members) < limits.min_patch_size:
            reasons.append(f"size {len(self.members)}")
        if self.transferrin_distance < limits.transferrin_clearance:
            reasons.append(f"transferrin {self.transferrin_distance:.1f} A")
        if self.glycan_distance < limits.glycan_clearance:
            reasons.append(f"glycan {self.glycan_distance:.1f} A")
        return tuple(reasons)

    def row(self) -> dict[str, object]:
        return {
            "centre": self.centre,
            "domain": self.domain,
            "size": len(self.members),
            "conserved": self.conserved,
            "divergent": self.divergent,
            "conservation": f"{self.conservation:.3f}",
            "acidic_count": len(self.acidic_members),
            "mean_rsa": f"{self.mean_rsa:.3f}",
            "transferrin_dist": f"{self.transferrin_distance:.1f}",
            "glycan_dist": f"{self.glycan_distance:.1f}",
            "acidic_members": " ".join(str(n) for n in self.acidic_members),
            "residues": " ".join(str(n) for n in self.members),
        }


def domain_of(number: int) -> str:
    return next((name for name, span in DOMAINS.items() if number in span), "outside")


def one_letter(residue: gemmi.Residue) -> str:
    info = gemmi.find_tabulated_residue(residue.name)
    return info.one_letter_code.upper() if info and info.is_amino_acid() else ""


def centroid(residue: gemmi.Residue) -> gemmi.Position | None:
    """`statistics.mean` is ~100x slower here -- it uses exact Fraction arithmetic,
    which a three-component coordinate average does not need."""
    atoms = [a.pos for a in residue if a.element != HYDROGEN]
    if not atoms:
        return None
    count = len(atoms)
    return gemmi.Position(
        sum(p.x for p in atoms) / count,
        sum(p.y for p in atoms) / count,
        sum(p.z for p in atoms) / count,
    )


def read_uniprot(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def read_sequence(path: Path) -> str:
    return str(read_uniprot(path)["sequence"]["value"])  # type: ignore[index]


def n_linked_sites(path: Path) -> tuple[int, ...]:
    """N-glycosylation positions, parsed rather than copied out of the entry."""
    features = read_uniprot(path).get("features", [])
    return tuple(
        sorted(
            f["location"]["start"]["value"]
            for f in features  # type: ignore[union-attr]
            if f["type"] == "Glycosylation" and N_LINKED in f.get("description", "")
        )
    )


def map_human_to_mouse(human: str, mouse: str) -> dict[int, tuple[int | None, str]]:
    """Human residue number -> (mouse number, mouse residue). Gaps map to None."""
    aligner = Align.PairwiseAligner(scoring="blastp", mode="global")
    alignment = max(aligner.align(human, mouse), key=lambda a: a.score)
    mapping: dict[int, tuple[int | None, str]] = {}
    for (h_start, h_end), (m_start, m_end) in zip(*alignment.aligned):
        for offset in range(h_end - h_start):
            mapping[h_start + offset + 1] = (m_start + offset + 1, mouse[m_start + offset])
    for position in range(1, len(human) + 1):
        mapping.setdefault(position, (None, "-"))
    return mapping


def sasa_by_residue(path: Path, keep_chains: set[str]) -> dict[int, float]:
    """Per-residue SASA computed on the assembly, so buried interfaces stay buried."""
    structure = MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    model = next(iter(structure))
    for chain in list(model):
        if chain.id not in keep_chains:
            model.detach_child(chain.id)
    ShrakeRupley().compute(model, level="R")
    first = next(iter(model))
    return {r.id[1]: float(r.sasa) for r in first if r.id[0] == " "}


def load_chain(path: Path, name: str) -> gemmi.Chain:
    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    return structure[0][name]


def receptor_residues(chain: gemmi.Chain, sasa: dict[int, float]) -> dict[int, Residue]:
    out: dict[int, Residue] = {}
    for residue in chain:
        code = one_letter(residue)
        centre = centroid(residue) if code in MAX_ASA else None
        if centre is None:
            continue
        number = residue.seqid.num
        out[number] = Residue(number, code, sasa.get(number, 0.0) / MAX_ASA[code], centre)
    return out


def mismatches_against_uniprot(chain: gemmi.Chain, sequence: str) -> dict[int, tuple[str, str]]:
    """Residue number -> (structure code, UniProt code) wherever the two disagree."""
    out: dict[int, tuple[str, str]] = {}
    for residue in chain:
        code = one_letter(residue)
        number = residue.seqid.num
        if not code or not 1 <= number <= len(sequence):
            continue
        if code != sequence[number - 1]:
            out[number] = (code, sequence[number - 1])
    return out


def contact_numbers(
    path: Path, receptor_chain: str, partner_chains: Collection[str]
) -> dict[int, int]:
    """Receptor residue number -> count of partner heavy atoms within the cutoff.

    Indexed with an EMPTY UnitCell, not `structure.cell`: with the real cell the
    search also returns symmetry images, which is a crystal contact rather than an
    assembly contact. On 6WRV (P 4 2 2, 256x256x128 A) that silently adds three
    residues to the epitope -- and 11RN agrees either way, so checking one
    structure would not catch it.
    """
    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    model = structure[0]
    search = gemmi.NeighborSearch(model, gemmi.UnitCell(), CONTACT_CUTOFF)
    for chain in model:
        if chain.name == receptor_chain:
            search.add_chain(chain, include_h=False)
    counts: Counter[int] = Counter()
    for chain in model:
        if chain.name not in partner_chains:
            continue
        for residue in chain:
            for atom in residue:
                if atom.element == HYDROGEN:
                    continue
                for mark in search.find_atoms(atom.pos, "\0", radius=CONTACT_CUTOFF):
                    found = mark.to_cra(model)
                    if one_letter(found.residue):
                        counts[found.residue.seqid.num] += 1
    return dict(counts)


def describe_epitope(
    label: str,
    contacts: dict[int, int],
    human: str,
    mapping: dict[int, tuple[int | None, str]],
    structure_codes: dict[int, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for number in sorted(contacts):
        if number > len(human):
            continue
        mouse_number, mouse_aa = mapping[number]
        rows.append(
            {
                "epitope": label,
                "residue": number,
                "human": human[number - 1],
                "structure": structure_codes.get(number, "?"),
                "mouse": mouse_aa,
                "mouse_number": mouse_number if mouse_number else "",
                "conserved": "yes" if human[number - 1] == mouse_aa else "no",
                "contact_atoms": contacts[number],
                "in_transferrin_site": "yes" if number in TRANSFERRIN_SITE else "no",
            }
        )
    return rows


def build_patch(
    centre: Residue,
    residues: dict[int, Residue],
    human: str,
    mapping: dict[int, tuple[int | None, str]],
    glycans: Collection[int],
    limits: Thresholds,
) -> Patch:
    members = [
        r for r in residues.values()
        if r.rsa >= limits.surface_min_rsa
        and centre.centroid.dist(r.centroid) <= limits.patch_radius
    ]
    in_uniprot = [m for m in members if m.number <= len(human)]
    conserved = sum(1 for m in in_uniprot if human[m.number - 1] == mapping[m.number][1])
    return Patch(
        centre=centre.number,
        members=tuple(sorted(m.number for m in members)),
        acidic_members=tuple(sorted(m.number for m in members if m.amino_acid in ACIDIC)),
        conserved=conserved,
        counted=len(in_uniprot),
        mean_rsa=st.mean(m.rsa for m in members) if members else 0.0,
        transferrin_distance=min_distance(centre.centroid, residues, TRANSFERRIN_SITE),
        glycan_distance=min_distance(centre.centroid, residues, glycans),
    )


def min_distance(
    centre: gemmi.Position, residues: dict[int, Residue], numbers: Collection[int]
) -> float:
    present = [residues[n].centroid for n in numbers if n in residues]
    return min((centre.dist(p) for p in present), default=float("inf"))


def report_brief_alignment(human: str, mapping: dict[int, tuple[int | None, str]]) -> None:
    differing = {n for n, aa in BRIEF_DIVERGENT.items() if mapping[n][1] != aa}
    unexpected = differing - INDEL_AMBIGUOUS
    if unexpected:
        raise typer.BadParameter(
            f"alignment disagrees with the brief at {sorted(unexpected)}, outside the known "
            f"210-214 indel ambiguity -- the alignment, not the brief, should be suspected"
        )
    note = (
        f"; {sorted(differing)} differ only in where the gap falls inside the 210-214 "
        "indel, and are non-conserved either way" if differing else ""
    )
    typer.echo(f"brief's six divergent apical residues: "
               f"{len(BRIEF_DIVERGENT) - len(differing)}/6 reproduce exactly{note}")
    identity = sum(1 for n in range(121, 761) if human[n - 1] == mapping[n][1]) / 640
    typer.echo(f"ectodomain identity, human vs mouse, 121-760: {identity:.1%}\n")


def report_epitope(label: str, rows: list[dict[str, object]]) -> None:
    conserved = sum(1 for r in rows if r["conserved"] == "yes")
    titratable = sum(1 for r in rows if r["human"] in TITRATABLE)
    typer.echo(
        f"{label:28s} {len(rows):3d} contacts  {conserved:3d} conserved "
        f"({conserved / len(rows):4.0%})  {titratable:2d} titratable (D/E/H)"
    )


def report_sensitivity(scanned: list[Patch], limits: Thresholds) -> None:
    typer.echo("\nsensitivity of the candidate count to the two binding gates:")
    for size in (8, 10, 12, 14):
        counts = []
        for clearance in (0.0, 6.0, 8.0):
            trial = Thresholds(limits.patch_radius, size, limits.surface_min_rsa,
                               clearance, limits.glycan_clearance)
            counts.append(sum(1 for p in scanned if not p.rejections(trial)))
        typer.echo(f"  size >= {size:2d}:  transferrin clearance 0 A: {counts[0]:3d}   "
                   f"6 A: {counts[1]:3d}   8 A: {counts[2]:3d}")


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@app.command()
def main(
    apical_complex: Annotated[Path, typer.Option("--apical", help="VP1u complex, 11RN")],
    designed: Annotated[Path, typer.Option(help="designed-binder complex, 6WRV")],
    dimer: Annotated[Path, typer.Option(help="unbound ectodomain homodimer, 11QC")],
    human_json: Annotated[Path, typer.Option(help="UniProt P02786 JSON")],
    mouse_json: Annotated[Path, typer.Option(help="UniProt Q62351 JSON")],
    out_contacts: Annotated[Path, typer.Option(help="per-contact CSV")],
    out_patches: Annotated[Path, typer.Option(help="surface-patch scan CSV")],
    patch_radius: Annotated[float, typer.Option(help="patch radius, A")] = 12.0,
    min_patch_size: Annotated[int, typer.Option(help="minimum exposed residues")] = 10,
    surface_min_rsa: Annotated[float, typer.Option(help="relative exposure cut")] = 0.25,
    transferrin_clearance: Annotated[float, typer.Option(help="minimum distance, A")] = 6.0,
    glycan_clearance: Annotated[float, typer.Option(help="minimum distance, A")] = 12.0,
) -> None:
    limits = Thresholds(patch_radius, min_patch_size, surface_min_rsa,
                        transferrin_clearance, glycan_clearance)
    human, mouse = read_sequence(human_json), read_sequence(mouse_json)
    mapping = map_human_to_mouse(human, mouse)
    glycans = n_linked_sites(human_json)
    typer.echo(f"N-glycosylation parsed from P02786: {list(glycans)}")
    report_brief_alignment(human, mapping)

    rows: list[dict[str, object]] = []
    for label, which, partners in EPITOPE_SOURCES:
        path = apical_complex if which == "apical" else designed
        chain = load_chain(path, "A")
        wrong = mismatches_against_uniprot(chain, human)
        codes = {r.seqid.num: one_letter(r) for r in chain if one_letter(r)}
        epitope = describe_epitope(label, contact_numbers(path, "A", partners),
                                   human, mapping, codes)
        bad = {n: v for n, v in wrong.items() if n in {int(r["residue"]) for r in epitope}}  # type: ignore[arg-type]
        if bad:
            raise typer.BadParameter(f"{path.name} disagrees with UniProt inside the epitope: {bad}")
        if wrong:
            typer.echo(f"  note: {path.name} differs from UniProt at {sorted(wrong)} "
                       f"(outside the epitope; the Ser142 allele is expected in 6WRV)")
        report_epitope(label, epitope)
        rows += epitope

    vp1u = [r for r in rows if r["epitope"] == EPITOPE_SOURCES[0][0]]
    conserved_core = sum(1 for r in vp1u if r["conserved"] == "yes")
    typer.echo(f"\nthe brief's apical patch, measured: {len(vp1u)} VP1u contacts "
               f"{sorted(int(r['residue']) for r in vp1u)}, "  # type: ignore[arg-type]
               f"{conserved_core} conserved ({conserved_core / len(vp1u):.0%})")
    write_rows(out_contacts, rows)

    sasa = sasa_by_residue(dimer, {"A", "B"})
    residues = receptor_residues(load_chain(dimer, "A"), sasa)
    surface = [r for r in residues.values() if r.rsa >= limits.surface_min_rsa]
    typer.echo(f"\n{len(surface)} of {len(residues)} residues exposed on the dimer "
               f"(RSA >= {limits.surface_min_rsa})")

    scanned = [build_patch(r, residues, human, mapping, glycans, limits) for r in surface]
    usable = sorted((p for p in scanned if not p.rejections(limits)),
                    key=lambda p: (-p.conservation, -len(p.acidic_members)))
    write_rows(out_patches, [p.row() for p in usable])
    out_patches.with_suffix(".thresholds.json").write_text(
        json.dumps(limits.as_dict(), indent=2) + "\n"
    )

    counted = Counter(reason.split()[0] for p in scanned for reason in p.rejections(limits))
    typer.echo(f"\n{len(usable)} usable patches of {len(scanned)} scanned; rejections by gate: "
               + ", ".join(f"{gate} {n}" for gate, n in counted.most_common()))
    typer.echo(f"\n{'centre':>6s} {'domain':>14s} {'size':>4s} {'cons':>5s} {'D/E':>4s} "
               f"{'rsa':>5s} {'Tf':>5s} {'glyc':>5s}")
    for patch in usable[:14]:
        typer.echo(f"{patch.centre:6d} {patch.domain:>14s} {len(patch.members):4d} "
                   f"{patch.conservation:5.0%} {len(patch.acidic_members):4d} "
                   f"{patch.mean_rsa:5.2f} {patch.transferrin_distance:5.1f} "
                   f"{patch.glycan_distance:5.1f}")
    report_sensitivity(scanned, limits)
    typer.echo(f"\nwrote {out_contacts}\nwrote {out_patches}")


if __name__ == "__main__":
    app()
