"""Reading designed complexes: which chain is the binder, and what it buries.

Chain identity is decided by TNF-alpha substring containment, not fractional identity.
BindCraft concatenates the target protomers into one chain, so identity/len(chain) falls to
~0.5 and an identity threshold misclassifies the target as the binder -- which silently
attributed target histidines to the design and zeroed every interface measurement.
"""

from __future__ import annotations

import copy
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from Bio import Align, PDB
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE_TO_ONE
from Bio.PDB.Residue import Residue
from Bio.PDB.SASA import ShrakeRupley

from tnf_binder.target import HUMAN, MOUSE

ResidueKey = tuple[str, int]  # (chain id, author residue number)
CLASH_DISTANCE = 2.2
_ALIGNER = Align.PairwiseAligner(scoring="blastp", open_gap_score=-11, extend_gap_score=-1)


def read_model(path: Path) -> Model:
    parser = (PDB.MMCIFParser(QUIET=True) if path.suffix in (".cif", ".mmcif")
              else PDB.PDBParser(QUIET=True))
    return parser.get_structure(path.stem, path)[0]


def chain_sequence(chain: Chain) -> str:
    return "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in chain if r.id[0] == " ")


def contains_tnf(sequence: str, window: int = 25) -> bool:
    """True if a long exact stretch of either ortholog's TNF-alpha occurs in the chain.

    Both orthologs, because mouse TNF contains no exact 25-mer of the human sequence, so a
    human-only test classifies the mouse target as the binder. The GPU-box copy guarded that
    with `identity_to_human > 0.8`, which could never fire: mouse scores 0.795 under its own
    normalisation (by the chain length) and 0.790 under `TnfTarget.identity`.
    """
    return any(sequence[i:i + window] in ortholog
               for ortholog in (HUMAN, MOUSE)
               for i in range(0, max(1, len(sequence) - window), 5))


def mature_numbering(chains: Sequence[Chain]) -> dict[Residue, int]:
    """Target residue -> mature TNF number, by alignment.

    Alignment rather than an offset: BoltzGen renumbers target chains to 1..N, and 1TNF carries
    Leu where UniProt P01375 has Asp143, so neither stored numbering nor substring matching
    survives every input.
    """
    numbering: dict[Residue, int] = {}
    for chain in chains:
        residues = [r for r in chain if r.id[0] == " "]
        sequence = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in residues)
        human_row, chain_row = (str(x) for x in _ALIGNER.align(HUMAN, sequence)[0])
        hi = si = 0
        for h_char, c_char in zip(human_row, chain_row):
            hi += h_char != "-"
            if c_char != "-":
                if h_char != "-":
                    numbering[residues[si]] = hi
                si += 1
    return numbering


def residue_sasa(model: Model, chain_ids: set[str]) -> dict[ResidueKey, float]:
    """Per-residue SASA with only `chain_ids` present, i.e. that subset taken alone."""
    subset = copy.deepcopy(model)
    for chain in list(subset):
        if chain.id not in chain_ids:
            subset.detach_child(chain.id)
    ShrakeRupley().compute(subset, level="R")
    return {(c.id, r.id[1]): r.sasa for c in subset for r in c if r.id[0] == " "}


class BinderChains(PDB.Select):
    """PDBIO selector keeping only the named chains."""

    def __init__(self, keep: set[str]) -> None:
        self.keep = keep

    def accept_chain(self, chain: Chain) -> bool:
        return chain.id in self.keep


@dataclass(frozen=True)
class Complex:
    """A designed binder bound to TNF-alpha, with the two roles separated."""

    model: Model
    target_chains: tuple[Chain, ...]
    binder_chains: tuple[Chain, ...]

    @classmethod
    def from_file(cls, path: Path) -> "Complex":
        return cls.from_model(read_model(path))

    @classmethod
    def from_model(cls, model: Model) -> "Complex":
        """Split an already-loaded model into target and binder chains.

        Separate from `from_file` because the GPU-box scripts load a model once and pass it
        around; without this they each kept their own copy of the chain-splitting rule.
        """
        target, binder = [], []
        for chain in model:
            sequence = chain_sequence(chain)
            if len(sequence) < 5:
                continue
            (target if contains_tnf(sequence) else binder).append(chain)
        return cls(model, tuple(target), tuple(binder))

    @property
    def binder_sequence(self) -> str:
        return "".join(chain_sequence(c) for c in self.binder_chains)

    @property
    def is_resolved(self) -> bool:
        return bool(self.target_chains and self.binder_chains)

    def mature_numbering(self) -> dict[Residue, int]:
        return mature_numbering(self.target_chains)

    def sasa(self, chain_ids: set[str]) -> dict[ResidueKey, float]:
        return residue_sasa(self.model, chain_ids)

    def buried_area(self) -> tuple[dict[ResidueKey, float], dict[ResidueKey, float]]:
        """(target-side, binder-side) buried area per residue, in A^2.

        Buried area rather than contact counts: counts scale with interpenetration, so
        clashing poly-alanine designs out-scored real ones on every ratio built from them.
        """
        target_ids = {c.id for c in self.target_chains}
        binder_ids = {c.id for c in self.binder_chains}
        bound = self.sasa(target_ids | binder_ids)
        target_free = self.sasa(target_ids)
        binder_free = self.sasa(binder_ids)
        target_buried = {k: v - bound.get(k, v) for k, v in target_free.items()
                         if v - bound.get(k, v) > 0.1}
        binder_buried = {k: v - bound.get(k, v) for k, v in binder_free.items()}
        return target_buried, binder_buried

    def buried_by_position(self) -> dict[int, float]:
        """Target-side buried area summed per mature TNF position.

        Residues the alignment could not place are dropped rather than guessed at.
        """
        target_buried, _ = self.buried_area()
        numbering = self.mature_numbering()
        position_of: dict[ResidueKey, int] = {
            (c.id, r.id[1]): numbering[r]
            for c in self.target_chains for r in c if r in numbering}
        areas: defaultdict[int, float] = defaultdict(float)
        for key, area in target_buried.items():
            if (position := position_of.get(key)) is not None:
                areas[position] += area
        return dict(areas)

    def write_binder(self, dest: Path) -> Path:
        """Write the binder chains alone, which is the only form a fold search can use.

        `touchstone.search.FoldSearch` must be queried with binder-only structures: a whole
        complex matches its own target, and reports a fold hit that says nothing about the
        design. Nothing else in the package produced one, so a fold-novelty lookup had no
        input it could legally take.
        """
        keep = {c.id for c in self.binder_chains}
        io = PDB.PDBIO()
        io.set_structure(self.model)
        io.save(str(dest), BinderChains(keep))
        return dest

    def clashes(self) -> int:
        """Heavy-atom pairs closer than 2.2 A, i.e. physically interpenetrating."""
        target_atoms = [a for c in self.target_chains for r in c if r.id[0] == " "
                        for a in r if a.element != "H"]
        search = PDB.NeighborSearch(target_atoms)
        return sum(len(search.search(atom.coord, CLASH_DISTANCE))
                   for chain in self.binder_chains
                   for residue in chain if residue.id[0] == " "
                   for atom in residue if atom.element != "H")
