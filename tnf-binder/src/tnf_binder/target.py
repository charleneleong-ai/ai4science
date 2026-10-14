"""TNF-alpha target: sequences, the receptor epitope, and the numbering that keeps biting.

Numbering is the dangerous part of this domain and so it is owned here rather than re-derived
per script. Three conventions coexist:

  mature      1..157 of the soluble domain (UniProt P01375-1 residues 77..233)
  1TNF file   6..157 -- residues 1-5 are disordered, so stored number == mature number
  BoltzGen    output chains renumbered 1..152, so stored number == mature - 5

An unguarded conversion between the last two silently shifted an epitope by five residues and
targeted the wrong surface. Every accessor here is asserted against residue identity.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

from Bio import Align

HUMAN = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTH"
         "TISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")
MOUSE = ("LRSSSQNSSDKPVAHVVANHQVEEQLEWLSQRANALLANGMDLKDNQLVVPADGLYLVYSQVLFKGQGCPDYVLLTH"
         "TVSRFAISYQEKVNLLSAVKSPCPKDTPEGAELKPWYEPIYLGGVFQLEKGDQLSAEVNLPKYLDFAESGQVYFGVIAL")

#: TNF residues within 4.5 A of TNFR2 in 3ALQ, mature numbering.
TNFR2_EPITOPE = (20, 21, 23, 31, 32, 33, 63, 67, 71, 72, 73, 75, 77, 85, 86, 87, 88, 89, 90,
                 91, 92, 97, 111, 113, 115, 128, 137, 143, 144, 145, 146, 147, 149)
#: Conserved, concave sub-patch used as the design epitope.
POCKET = (113, 115, 143, 144, 146, 147)
#: The original conserved lobe, carrying K90.
LOBE = (86, 87, 88, 90, 91, 92)
#: Generation hotspots of the first campaign: the lobe plus the second-protomer contacts.
#: Kept because the geometry tables of that campaign report focus as a share of these, and a
#: number cannot be compared against those tables under any other definition.
HOTSPOTS = (33, 86, 87, 88, 90, 91, 92, 147)

_ALIGNER = Align.PairwiseAligner(scoring="blastp", open_gap_score=-11, extend_gap_score=-1)


@dataclass(frozen=True)
class Ortholog:
    """One species' soluble TNF-alpha domain."""

    name: str
    sequence: str

    def residue(self, mature_number: int) -> str:
        return self.sequence[mature_number - 1]


HUMAN_TNF = Ortholog("human", HUMAN)
MOUSE_TNF = Ortholog("mouse", MOUSE)


class TnfTarget:
    """Epitope definition and the human/mouse divergence that constrains it."""

    def __init__(self, human: Ortholog = HUMAN_TNF, mouse: Ortholog = MOUSE_TNF) -> None:
        self.human = human
        self.mouse = mouse

    @cached_property
    def _alignment(self) -> tuple[str, str]:
        best = _ALIGNER.align(self.human.sequence, self.mouse.sequence)[0]
        return tuple(str(row) for row in best)  # type: ignore[return-value]

    @cached_property
    def human_to_mouse(self) -> dict[int, int]:
        """Mature human number -> mature mouse number, for aligned positions only."""
        human_row, mouse_row = self._alignment
        mapping, hi, mi = {}, 0, 0
        for h_char, m_char in zip(human_row, mouse_row):
            hi += h_char != "-"
            mi += m_char != "-"
            if h_char != "-" and m_char != "-":
                mapping[hi] = mi
        return mapping

    @cached_property
    def divergent(self) -> frozenset[int]:
        """Human positions whose mouse counterpart differs, or has no counterpart."""
        human_row, mouse_row = self._alignment
        out, hi = set(), 0
        for h_char, m_char in zip(human_row, mouse_row):
            if h_char == "-":
                continue
            hi += 1
            if h_char != m_char:
                out.add(hi)
        return frozenset(out)

    @cached_property
    def conserved_epitope(self) -> tuple[int, ...]:
        return tuple(p for p in TNFR2_EPITOPE if p not in self.divergent)

    @property
    def identity(self) -> float:
        """Fraction identical over the aligned soluble domain."""
        human_row, mouse_row = self._alignment
        same = sum(1 for h, m in zip(human_row, mouse_row) if h == m and h != "-")
        return same / len(self.human.sequence)

    def mouse_equivalent(self, positions: tuple[int, ...]) -> tuple[int, ...]:
        """Mouse numbering for hotspots, asserting the residue identity is preserved.

        Dual-ortholog design needs hotspots in each target's own numbering, and a silent
        off-by-one here would aim the mouse campaign at a different surface.
        """
        out = []
        for p in positions:
            m = self.human_to_mouse.get(p)
            if m is None:
                raise KeyError(f"human {p} has no aligned mouse position")
            if self.human.residue(p) != self.mouse.residue(m):
                raise ValueError(
                    f"human {p}{self.human.residue(p)} maps to mouse {m}"
                    f"{self.mouse.residue(m)}: not conserved, so not a valid shared hotspot")
            out.append(m)
        return tuple(out)

    @staticmethod
    def chain_position(mature_number: int, first_stored: int) -> int:
        """1-based position within a chain as stored in a file.

        BoltzGen's parse_range indexes chain POSITION, not residue number. With 1TNF starting
        at residue 6, writing mature numbers directly shifts the epitope by five.
        """
        return mature_number - first_stored + 1
