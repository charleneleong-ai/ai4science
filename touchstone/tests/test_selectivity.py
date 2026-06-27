import numpy as np
import pytest

pytest.importorskip("ase")

from ase import Atoms  # noqa: E402
from ase.calculators.calculator import Calculator, all_changes  # noqa: E402

from touchstone import BinderDesign, MLIPSelectivityVerifier  # noqa: E402
from touchstone.core import CoordinationSite, element_symbol  # noqa: E402

_OCT = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
_ELEMS = ("N", "N", "O", "O", "N", "O")


class MetalBiasSpring(Calculator):
    """Springs donors to the metal at r0, plus a per-metal per-bond depth — so the
    binding energy ΔE differs by metal (Ni binds strongest here). A deterministic
    stand-in for an MLIP with real metal discrimination."""

    implemented_properties = ["energy", "forces"]
    _DEPTH = {"Ni": -1.0, "Co": -0.7, "Cu": -0.5}

    def __init__(self, r0: float = 2.1, k: float = 5.0, cutoff: float = 2.8, **kw):
        super().__init__(**kw)
        self.r0, self.k, self.cutoff = r0, k, cutoff

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        syms, pos = atoms.get_chemical_symbols(), atoms.get_positions()
        e, f = 0.0, np.zeros_like(pos)
        for mi, ms in enumerate(syms):
            if ms not in self._DEPTH:
                continue
            for i, s in enumerate(syms):
                d = pos[i] - pos[mi]
                r = float(np.linalg.norm(d))
                if i == mi or s not in ("N", "O", "S") or r < 1e-9:
                    continue
                e += 0.5 * self.k * (r - self.r0) ** 2
                g = self.k * (r - self.r0) * (d / r)
                f[i] -= g
                f[mi] += g
                if r <= self.cutoff:
                    e += self._DEPTH[ms]  # per-bond depth → metal-dependent ΔE
        self.results = {"energy": float(e), "forces": f}


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


class TestMLIPSelectivity:
    def test_profile_ranks_by_binding_energy(self, tmp_path):
        prof = MLIPSelectivityVerifier(calculator=MetalBiasSpring()).profile(_design(tmp_path, "Ni2+"))
        assert prof.preferred == "Ni2+" and prof.margin > 0  # Ni binds strongest

    @pytest.mark.parametrize(
        "metal, expect_trust", [("Ni2+", True), ("Cu2+", False)],
        ids=["target-binds-strongest", "competitor-binds-stronger"],
    )
    def test_trust_follows_binding_preference(self, tmp_path, metal, expect_trust):
        # Ni binds strongest in the calc: a Ni-target design is selective (trust),
        # a Cu-target one is out-competed (weak — disagreement, not defer)
        v = MLIPSelectivityVerifier(calculator=MetalBiasSpring()).verify(_design(tmp_path, metal))
        assert v.trust == expect_trust and (v.score > 0.5) == expect_trust and not v.ood

    def test_failure_defers(self, tmp_path):
        v = MLIPSelectivityVerifier(calculator=_Exploding()).verify(_design(tmp_path, "Ni2+"))
        assert v.ood and not v.trust and "failed" in v.reason

    def test_target_outside_panel_does_not_crash(self, tmp_path):
        # design targets Zn2+, not in the default Ni/Cu/Co panel — must produce a verdict,
        # not a KeyError from the margin lookup
        v = MLIPSelectivityVerifier(calculator=MetalBiasSpring()).verify(_design(tmp_path, "Zn2+"))
        assert v.label in ("trust", "weak", "defer") and not v.trust
