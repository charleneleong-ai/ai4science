import numpy as np
import pytest

pytest.importorskip("ase")

from ase import Atoms  # noqa: E402
from ase.calculators.calculator import Calculator, all_changes  # noqa: E402

from touchstone import BinderDesign, TrsVerifier  # noqa: E402
from touchstone.core import CoordinationSite, element_symbol  # noqa: E402

_OCT = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
_ELEMS = ("N", "N", "O", "O", "N", "O")


class OffsetSpring(Calculator):
    """Pulls each donor toward its start position + `dx` Å along x — so relax_apo (metal removed)
    yields a controllable donor drift, exercising the trust/weak/defer thresholds deterministically."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, dx: float = 0.0, k: float = 5.0, **kw):
        super().__init__(**kw)
        self.dx, self.k, self._target = dx, k, None

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        pos, syms = atoms.get_positions(), atoms.get_chemical_symbols()
        if self._target is None:  # freeze targets at the first eval (apo start + offset)
            self._target = pos.copy()
            for i, s in enumerate(syms):
                if s in ("N", "O", "S"):
                    self._target[i] = pos[i] + (self.dx, 0.0, 0.0)
        d = pos - self._target
        self.results = {"energy": 0.5 * self.k * float(np.sum(d * d)), "forces": -self.k * d}


class _Exploding(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, *a, **k):
        raise RuntimeError("diverged")


def _design(tmp_path, metal: str = "Ni2+") -> BinderDesign:
    el = element_symbol(metal)
    atoms = Atoms(symbols=[el, *_ELEMS], positions=np.vstack([[0, 0, 0], _OCT * 2.1]))
    src = str(tmp_path / "d.pdb")
    atoms.write(src)
    site = CoordinationSite(metal, np.zeros(3), atoms.get_positions()[1:], _ELEMS)
    return BinderDesign("SEQ", site, generator="test", generator_confidence=0.5, source=src)


class TestTrs:
    @pytest.mark.parametrize(
        "dx, label", [(0.0, "trust"), (2.0, "weak"), (4.0, "defer")],
        ids=["preorganized", "reorganizes", "collapses"],
    )
    def test_verdict_follows_donor_reorganization(self, tmp_path, dx, label):
        # donors held at their holo positions (dx=0) ⇒ preorganized (trust); pulled far (dx=4) ⇒
        # the shell collapses off the manifold (defer); in between ⇒ weak
        v = TrsVerifier(calculator=OffsetSpring(dx=dx)).verify(_design(tmp_path))
        assert v.label == label

    def test_relaxation_failure_defers(self, tmp_path):
        v = TrsVerifier(calculator=_Exploding()).verify(_design(tmp_path))
        assert v.ood and not v.trust and "failed" in v.reason
