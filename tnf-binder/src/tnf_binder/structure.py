"""Reading designed complexes: which chain is the binder, and what it buries.

Chain identity is decided by TNF-alpha substring containment, not fractional identity.
BindCraft concatenates the target protomers into one chain, so identity/len(chain) falls to
~0.5 and an identity threshold misclassifies the target as the binder -- which silently
attributed target histidines to the design and zeroed every interface measurement.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path

from Bio import Align, PDB
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE_TO_ONE
from Bio.PDB.SASA import ShrakeRupley

from tnf_binder.target import HUMAN

CLASH_DISTANCE = 2.2
_ALIGNER = Align.PairwiseAligner(scoring="blastp", open_gap_score=-11, extend_gap_score=-1)


def read_model(path: Path):
    parser = (PDB.MMCIFParser(QUIET=True) if path.suffix in (".cif", ".mmcif")
              else PDB.PDBParser(QUIET=True))
    return parser.get_structure(path.stem, path)[0]


def chain_sequence(chain) -> str:
    return "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in chain if r.id[0] == " ")


def contains_tnf(sequence: str, window: int = 25) -> bool:
    """True if a long exact stretch of human TNF-alpha occurs anywhere in the chain."""
    return any(sequence[i:i + window] in HUMAN
               for i in range(0, max(1, len(sequence) - window), 5))


@dataclass(frozen=True)
class Complex:
    """A designed binder bound to TNF-alpha, with the two roles separated."""

    model: object
    target_chains: tuple
    binder_chains: tuple

    @classmethod
    def from_file(cls, path: Path) -> "Complex":
        model = read_model(path)
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

    def mature_numbering(self) -> dict:
        """Target residue -> mature TNF number, by alignment.

        Alignment rather than an offset: BoltzGen renumbers target chains to 1..N, and 1TNF
        carries Leu where UniProt P01375 has Asp143, so neither stored numbering nor substring
        matching survives every input.
        """
        numbering = {}
        for chain in self.target_chains:
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

    def _sasa(self, chain_ids: set[str]) -> dict[tuple[str, int], float]:
        subset = copy.deepcopy(self.model)
        for chain in list(subset):
            if chain.id not in chain_ids:
                subset.detach_child(chain.id)
        ShrakeRupley().compute(subset, level="R")
        return {(c.id, r.id[1]): r.sasa for c in subset for r in c if r.id[0] == " "}

    def buried_area(self) -> tuple[dict, dict]:
        """(target-side, binder-side) buried area per residue, in A^2.

        Buried area rather than contact counts: counts scale with interpenetration, so
        clashing poly-alanine designs out-scored real ones on every ratio built from them.
        """
        target_ids = {c.id for c in self.target_chains}
        binder_ids = {c.id for c in self.binder_chains}
        bound = self._sasa(target_ids | binder_ids)
        target_free = self._sasa(target_ids)
        binder_free = self._sasa(binder_ids)
        target_buried = {k: v - bound.get(k, v) for k, v in target_free.items()
                         if v - bound.get(k, v) > 0.1}
        binder_buried = {k: v - bound.get(k, v) for k, v in binder_free.items()}
        return target_buried, binder_buried

    def clashes(self) -> int:
        """Heavy-atom pairs closer than 2.2 A, i.e. physically interpenetrating."""
        target_atoms = [a for c in self.target_chains for r in c if r.id[0] == " "
                        for a in r if a.element != "H"]
        search = PDB.NeighborSearch(target_atoms)
        return sum(len(search.search(atom.coord, CLASH_DISTANCE))
                   for chain in self.binder_chains
                   for residue in chain if residue.id[0] == " "
                   for atom in residue if atom.element != "H")
