import numpy as np
import pytest

pytest.importorskip("ase")

from ase import Atoms  # noqa: E402
from ase.calculators.calculator import Calculator, all_changes  # noqa: E402

from touchstone import BinderDesign, MLIPSelectivityVerifier  # noqa: E402
from touchstone.core import CoordinationSite, element_symbol  # noqa: E402
from touchstone.physics.mlip import multiplicity  # noqa: E402
from touchstone.physics.selectivity import ranks_irving_williams, swap_metal  # noqa: E402

_OCT = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
_ELEMS = ("N", "N", "O", "O", "N", "O")


class MetalBiasSpring(Calculator):
    """Springs donors to the metal at r0, plus a per-metal per-bond depth — so ΔE differs by
    metal. Depths follow the **Irving–Williams series** (Cu²⁺ the peak), so this stand-in *passes*
    the tier's validity gate: a deterministic model of an MLIP that has real metal discrimination."""

    implemented_properties = ["energy", "forces"]
    DEPTH = {"Mn": -0.30, "Fe": -0.40, "Co": -0.50, "Ni": -0.70, "Cu": -1.00, "Zn": -0.20}

    def __init__(self, r0: float = 2.1, k: float = 5.0, cutoff: float = 2.8, **kw):
        super().__init__(**kw)
        self.r0, self.k, self.cutoff = r0, k, cutoff

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        syms, pos = atoms.get_chemical_symbols(), atoms.get_positions()
        e, f = 0.0, np.zeros_like(pos)
        for mi, ms in enumerate(syms):
            if ms not in self.DEPTH:
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
                    e += self.DEPTH[ms]  # per-bond depth → metal-dependent ΔE
        self.results = {"energy": float(e), "forces": f}


class SpinBlindSpring(MetalBiasSpring):
    """The MACE-MP failure mode: no charge/spin state, so no ligand-field stabilisation — it
    *inverts* Irving–Williams and ranks Mn²⁺ strongest. The gate exists to catch exactly this."""

    DEPTH = {"Mn": -1.00, "Fe": -0.80, "Co": -0.60, "Ni": -0.40, "Cu": -0.50, "Zn": -0.20}


class CuPeakScrambledSpring(MetalBiasSpring):
    """Cu²⁺ on top, but the rising limb is nonsense (Mn²⁺ second-strongest, above Fe/Co/Ni).

    This is the backbone a *peak-only* gate waves through: `argmin == Cu2+` is one bit, and one bit
    is cheap to satisfy by accident. Irving–Williams is a **series**, so the gate checks the series."""

    DEPTH = {"Mn": -0.90, "Fe": -0.30, "Co": -0.50, "Ni": -0.60, "Cu": -1.00, "Zn": -0.20}


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


class TestIrvingWilliamsGate:
    """A metal ranking from a spin-blind potential isn't a weak signal — it's a meaningless one."""

    @pytest.mark.parametrize(
        "calc_cls, valid",
        [(MetalBiasSpring, True), (SpinBlindSpring, False), (CuPeakScrambledSpring, False)],
        ids=["reproduces-the-series", "inverts-the-series-like-MACE", "Cu-peak-but-scrambled-limb"],
    )
    def test_gate_requires_the_whole_series_not_just_the_peak(self, calc_cls, valid):
        # the scrambled case is the point: it puts Cu2+ at the peak, so a one-bit `argmin == Cu`
        # gate would pass it while its Mn > Fe > Co ordering is still physically meaningless
        assert ranks_irving_williams(calc_cls()) is valid

    def test_failing_backbone_makes_the_tier_refuse_to_judge(self, tmp_path):
        # rather than emit a ranking it cannot justify, the tier defers and says why
        v = MLIPSelectivityVerifier(calculator=SpinBlindSpring()).verify(_design(tmp_path, "Ni2+"))
        assert v.ood and not v.trust and "Irving" in v.reason


class TestSpinAndCharge:
    """Spin state *is* the ligand-field physics that decides metal preference, and it depends on the
    oxidation state, not the element. Getting it wrong doesn't fail loudly — it moves the energy."""

    @pytest.mark.parametrize(
        "ion, mult",
        [("Ni2+", 3), ("Cu2+", 2), ("Cu1+", 1), ("Fe2+", 5), ("Fe3+", 6), ("Mn2+", 6), ("Zn2+", 1)],
    )
    def test_multiplicity_follows_the_oxidation_state(self, ion, mult):
        # Fe2+ is d6 (4 unpaired) and Fe3+ is d5 (5): an element-keyed table hands Fe3+ Fe2+'s spin
        assert multiplicity(ion) == mult

    @pytest.mark.parametrize("ion", ["Co3+", "Ru2+", "La3+"])
    def test_untabulated_ions_have_no_spin_rather_than_a_default(self, ion):
        # Co3+ is d6 — high-spin in weak fields, low-spin in most complexes. We don't know, so we
        # don't say. A singlet default would be a silent wrong answer.
        assert multiplicity(ion) is None

    def test_tier_defers_on_an_ion_it_cannot_assign_a_spin(self, tmp_path):
        v = MLIPSelectivityVerifier(calculator=MetalBiasSpring(), metals=("Ni2+", "Co3+")).verify(
            _design(tmp_path, "Ni2+")
        )
        assert v.ood and not v.trust and "Co3+" in v.reason and "spin" in v.reason

    def test_swap_carries_charge_and_spin_with_the_ion(self):
        # a panel mixing oxidation states must move the cluster's total charge too, or the apo leg
        # silently absorbs the difference and every ΔE in the panel is wrong
        atoms = Atoms(symbols=["Fe", "O"], positions=[[0, 0, 0], [2.1, 0, 0]])
        atoms.info.update(charge=0, spin=multiplicity("Fe2+"), ion="Fe2+")
        swapped = swap_metal(atoms, "Fe", "Fe3+")
        assert swapped.info["ion"] == "Fe3+"
        assert swapped.info["spin"] == 6  # d5, not Fe2+'s d6
        assert swapped.info["charge"] == 1  # 0 + (3 - 2)


class TestMLIPSelectivity:
    def test_profile_ranks_by_binding_energy(self, tmp_path):
        prof = MLIPSelectivityVerifier(calculator=MetalBiasSpring()).profile(_design(tmp_path, "Cu2+"))
        assert prof.preferred == "Cu2+" and prof.margin > 0  # Cu binds strongest (Irving–Williams)

    @pytest.mark.parametrize(
        "metal, expect_trust", [("Cu2+", True), ("Ni2+", False)],
        ids=["target-is-the-preferred-metal", "target-out-competed-by-Cu"],
    )
    def test_trust_follows_binding_preference(self, tmp_path, metal, expect_trust):
        # Cu binds strongest in this calc: a Cu-target design is selective (trust); a Ni-target
        # one is out-competed (weak — a disagreement, not a defer)
        v = MLIPSelectivityVerifier(calculator=MetalBiasSpring()).verify(_design(tmp_path, metal))
        assert v.trust == expect_trust and (v.score > 0.5) == expect_trust and not v.ood

    def test_failure_defers(self, tmp_path):
        v = MLIPSelectivityVerifier(calculator=_Exploding()).verify(_design(tmp_path, "Ni2+"))
        assert v.ood and not v.trust

    def test_target_outside_panel_does_not_crash(self, tmp_path):
        # design targets Zn2+, not in the default Ni/Cu/Co panel — must produce a verdict,
        # not a KeyError from the margin lookup
        v = MLIPSelectivityVerifier(calculator=MetalBiasSpring()).verify(_design(tmp_path, "Zn2+"))
        assert v.label in ("trust", "weak", "defer") and not v.trust
