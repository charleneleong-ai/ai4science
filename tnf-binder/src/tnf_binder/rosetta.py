"""Histidine substitution with a real, repacked side chain, via PyRosetta.

Needed because propka cannot score an incomplete residue: renaming a residue to HIS and keeping
only backbone + CB leaves no imidazole, so propka silently omits it and every pKa comes back
missing. PyRosetta builds the side chain and repacks the local neighbourhood, so the predicted
pKa reflects a physically plausible rotamer in its actual environment.

Runs in the BindCraft conda env, which is where pyrosetta lives. PYTHONNOUSERSITE=1 is required
there: ~/.local numpy shadows that env otherwise.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

INIT_FLAGS = "-mute all -ignore_unrecognized_res -ignore_zero_occupancy false"


class HisMutator:
    """One PyRosetta session: initialise once, then mutate as many poses as asked.

    Stateful on purpose. `pyrosetta.init` is process-global and must not be called per mutation,
    and the score function it needs is expensive to build, so both are held here instead of
    being re-created by every caller.
    """

    def __init__(self, repack_radius: float = 6.0) -> None:
        self.repack_radius = repack_radius
        self._rosetta: Any = None
        self._score_function: Any = None

    @property
    def rosetta(self) -> Any:
        # Imported lazily: pyrosetta is licensed, heavy, and installed only on the design box,
        # while the rest of the package must stay importable without it.
        if self._rosetta is None:
            import pyrosetta

            pyrosetta.init(INIT_FLAGS, silent=True)
            self._rosetta = pyrosetta
            self._score_function = pyrosetta.get_score_function()
        return self._rosetta

    def load(self, pdb: Path) -> Any:
        return self.rosetta.pose_from_pdb(str(pdb))

    def mutate(self, pose: Any, sites: Sequence[tuple[str, int]]) -> bool:
        """Mutate every (chain, residue number) to HIS in place. False if any site is missing.

        All-or-nothing by design: a partially-mutated pose would be scored as if it carried the
        whole combination, and multi-histidine pKa values do not compose.
        """
        rosetta = self.rosetta
        info = pose.pdb_info()
        for chain, number in sites:
            index = info.pdb2pose(chain, number)
            if index == 0:
                return False
            rosetta.toolbox.mutants.mutate_residue(pose, index, "H", self.repack_radius,
                                                   self._score_function)
            if pose.residue(index).name3() != "HIS":
                return False
        return True

    @staticmethod
    def write(pose: Any, dest: Path) -> None:
        pose.dump_pdb(str(dest))

    @staticmethod
    def chain_sequence(pose: Any, chain: str) -> str:
        info = pose.pdb_info()
        return "".join(pose.residue(i).name1()
                       for i in range(1, pose.total_residue() + 1)
                       if info.chain(i) == chain)

    @staticmethod
    def binder_chain(pose: Any) -> str:
        """The shortest chain is the binder; every TNF-alpha protomer is 152 aa."""
        info = pose.pdb_info()
        lengths: dict[str, int] = {}
        for i in range(1, pose.total_residue() + 1):
            chain = info.chain(i)
            lengths[chain] = lengths.get(chain, 0) + 1
        return min(lengths, key=lambda c: lengths[c])
