"""Generator input: the numbering guard, and the co-folding spec text.

BoltzGen's `parse_range` is 1-based within the chain AS STORED IN THE FILE, not by residue
number. 1TNF starts at residue 6, so writing mature numbers straight into a spec shifts the
epitope by five and aims the campaign at a different surface -- which happened once, silently.
`StoredNumbering` reads the mapping off the structure and asserts the residues it lands on,
instead of doing the arithmetic and hoping.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE_TO_ONE

from tnf_binder.structure import ResidueKey, read_model


def collapse_ranges(numbers: Iterable[int]) -> str:
    """Sorted ints in BoltzGen's 'a..b,c' range syntax.

    Deduplicated: a repeated position would be emitted as "5,5..6", which names position 5
    twice and is malformed input rather than a harmless redundancy.
    """
    ordered, out, i = sorted(set(numbers)), [], 0
    while i < len(ordered):
        j = i
        while j + 1 < len(ordered) and ordered[j + 1] == ordered[j] + 1:
            j += 1
        out.append(f"{ordered[i]}..{ordered[j]}" if j > i else f"{ordered[i]}")
        i = j + 1
    return ",".join(out)


@dataclass(frozen=True)
class StoredNumbering:
    """Chain positions and residue identities as one structure file stores them."""

    positions: dict[ResidueKey, int]
    residues: dict[ResidueKey, str]

    @classmethod
    def from_file(cls, path: Path) -> "StoredNumbering":
        positions: dict[ResidueKey, int] = {}
        residues: dict[ResidueKey, str] = {}
        for chain in read_model(path):
            modelled = (r for r in chain if r.id[0] == " ")
            for index, residue in enumerate(modelled, start=1):
                key = (chain.id, residue.id[1])
                positions[key] = index
                residues[key] = THREE_TO_ONE.get(residue.get_resname(), "X")
        return cls(positions, residues)

    def sequence(self, chain: str, mature_numbers: Sequence[int]) -> str:
        return "".join(self.residues[(chain, n)] for n in mature_numbers)

    def ranges(self, chain: str, mature_numbers: Sequence[int],
               expect: str | None = None) -> str:
        """Chain-position ranges, refusing to emit them if the residues are not what `expect` says.

        The expected letters are written out rather than read from `target.HUMAN`, because the
        target files are 1TNF and 1TNF carries Leu where UniProt P01375 has Asp143 -- deriving
        the expectation from the canonical sequence would fail on a correct file.

        `expect` is omitted only for anti-hotspots, which steer the generator away from a
        surface rather than onto one; a shift there costs a worse spec, not a wrong epitope.
        """
        found = self.sequence(chain, mature_numbers)
        if expect is not None and found != expect:
            raise ValueError(
                f"numbering guard failed on chain {chain}: mature {list(mature_numbers)} "
                f"spells {found!r}, expected {expect!r}")
        return collapse_ranges(self.positions[(chain, n)] for n in mature_numbers)


def chain_block(chain_id: str, sequence: str) -> str:
    return f"  - protein:\n      id: {chain_id}\n      sequence: {sequence}\n      msa: empty"


def cofold_spec(target: str, *, binder: str | None = None, protomers: int = 3) -> str:
    """Boltz spec text: the binder as chain X, the target as protomers A..C.

    `msa: empty` throughout -- every sequence is given explicitly, so there is no MSA server
    round-trip to wait on.
    """
    chains = [] if binder is None else [chain_block("X", binder)]
    chains += [chain_block(chain_id, target) for chain_id in "ABC"[:protomers]]
    return "version: 1\nsequences:\n" + "\n".join(chains) + "\n"
