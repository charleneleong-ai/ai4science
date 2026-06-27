import numpy as np
import pytest

pytest.importorskip("ase")  # the MLIP tier rides on the optional [mace]/[uma] extra

from ase import Atoms  # noqa: E402
from ase.calculators.calculator import Calculator, all_changes  # noqa: E402

from touchstone import (  # noqa: E402
    BinderDesign,
    MLIPDynamicsVerifier,
    MLIPVerifier,
    make_backbone,
    relax_site,
)
from touchstone.core import CoordinationSite  # noqa: E402

_OCT = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


class SpringToMetal(Calculator):
    """Every non-metal atom is sprung to `r0` from the metal — a deterministic
    stand-in for an MLIP, so relaxations have a known ground truth."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, r0: float = 2.1, k: float = 5.0, metal: str = "Ni", **kw):
        super().__init__(**kw)
        self.r0, self.k, self.metal = r0, k, metal

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        syms, pos = atoms.get_chemical_symbols(), atoms.get_positions()
        e, f = 0.0, np.zeros_like(pos)
        if self.metal in syms:
            mi = syms.index(self.metal)
            for i in range(len(pos)):
                d = pos[i] - pos[mi]
                r = float(np.linalg.norm(d))
                if i == mi or r < 1e-9:
                    continue
                e += 0.5 * self.k * (r - self.r0) ** 2
                g = self.k * (r - self.r0) * (d / r)
                f[i] -= g
                f[mi] += g
        self.results = {"energy": float(e), "forces": f}


class _Exploding(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, *a, **k):
        raise RuntimeError("diverged")


def _cluster(r: float, elems=("N", "N", "O", "O", "N", "O")) -> Atoms:
    pos = np.vstack([[0, 0, 0], _OCT * r])
    return Atoms(symbols=["Ni", *elems], positions=pos)


def _design(tmp_path, r: float = 2.5, source: bool = True) -> BinderDesign:
    atoms = _cluster(r)
    src = None
    if source:
        src = str(tmp_path / "design.pdb")
        atoms.write(src)
    site = CoordinationSite(
        "Ni2+", np.zeros(3), atoms.get_positions()[1:], tuple(atoms.get_chemical_symbols()[1:])
    )
    return BinderDesign("SEQ", site, generator="test", generator_confidence=0.5, source=src)


class TestRelaxSite:
    def test_relaxes_donors_to_the_potential_minimum(self):
        r = relax_site(_cluster(2.5), SpringToMetal(r0=2.1), metal="Ni", interaction=False)
        assert r.cn_before == r.cn_after == 6
        assert abs(r.mean_bond - 2.1) < 0.05  # donors pulled to the spring length
        assert r.site_drift > 0.3  # they had to move ~0.4 Å to get there

    def test_reports_donors_that_leave_the_shell(self):
        # spring length past the 2.8 Å cutoff ⇒ the whole shell relaxes out
        r = relax_site(_cluster(2.6), SpringToMetal(r0=3.5), metal="Ni", interaction=False)
        assert r.cn_before == 6 and r.cn_after == 0 and r.donors_lost == 6

    def test_interaction_energy_is_finite(self):
        r = relax_site(_cluster(2.3), SpringToMetal(r0=2.1), metal="Ni", interaction=True)
        assert r.interaction_energy is not None and np.isfinite(r.interaction_energy)

    def test_missing_metal_raises(self):
        with pytest.raises(ValueError, match="no Ni"):
            relax_site(Atoms("N2", positions=[[0, 0, 0], [0, 0, 2]]), SpringToMetal(), metal="Ni")


class TestRestraint:
    def test_freezes_the_scaffold_so_donors_cannot_disperse(self):
        # metal + 3 in-shell N donors + 2 out-of-shell "backbone" O atoms (the
        # scaffold), placed symmetrically so they pull the free metal equally. The
        # spring drags every non-metal toward the metal; frozen, the scaffold holds.
        pos = np.vstack([[0, 0, 0], _OCT[:3] * 2.5, [[4.0, 0, 0], [-4.0, 0, 0]]])
        restrained = Atoms(symbols=["Ni", "N", "N", "N", "O", "O"], positions=pos)
        free = restrained.copy()
        r = relax_site(restrained, SpringToMetal(r0=2.1), metal="Ni", interaction=False, restrain=True)
        relax_site(free, SpringToMetal(r0=2.1), metal="Ni", interaction=False, restrain=False)

        assert np.allclose(restrained.get_positions()[4], [4.0, 0, 0], atol=1e-3)  # scaffold held
        assert np.linalg.norm(free.get_positions()[4]) < 3.0  # unrestrained: dragged in
        # the in-shell donors still relax to the spring length against the fixed scaffold
        assert r.cn_after == 3 and abs(r.mean_bond - 2.1) < 0.1


class TestMLIPVerifier:
    def test_trusts_a_site_that_holds(self, tmp_path):
        v = MLIPVerifier(calculator=SpringToMetal(r0=2.1)).verify(_design(tmp_path, r=2.5))
        assert v.trust and not v.ood and v.score > 0.5

    def test_defers_when_the_shell_collapses(self, tmp_path):
        v = MLIPVerifier(calculator=SpringToMetal(r0=3.5)).verify(_design(tmp_path, r=2.6))
        assert not v.trust and v.ood and v.score == 0.0  # cn_after 0 ⇒ defer

    @pytest.mark.parametrize(
        "calc, source",
        [(_Exploding(), True), (SpringToMetal(), False)],
        ids=["relaxation-raises", "missing-source"],
    )
    def test_defers_when_it_cannot_run(self, tmp_path, calc, source):
        # both reach the single except branch in verify(): a diverging relaxation
        # and a missing source structure
        v = MLIPVerifier(calculator=calc).verify(_design(tmp_path, source=source))
        assert v.ood and not v.trust and "failed" in v.reason

    def test_cluster_is_finite_not_periodic(self, tmp_path):
        # BoltzGen PDBs carry a CRYST1 1 Å cell ⇒ ase.io.read marks them periodic; on
        # a 1 Å cell an MLIP builds a vast periodic neighbour list (millions of image
        # edges ⇒ OOM). The cluster must strip the cell to a finite fragment.
        src = tmp_path / "cryst.pdb"
        src.write_text(
            "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1\n"
            "ATOM      1 NI   NI  A   1       0.000   0.000   0.000  1.00  0.00          NI\n"
            "ATOM      2  N   HIS A   2       2.000   0.000   0.000  1.00  0.00           N\n"
            "END\n"
        )
        site = CoordinationSite("Ni2+", np.zeros(3), np.array([[2.0, 0, 0]]), ("N",))
        d = BinderDesign("S", site, generator="t", generator_confidence=0.5, source=str(src))
        cluster = MLIPVerifier(calculator=SpringToMetal(), radius=5.0)._cluster(d)
        assert not cluster.pbc.any()  # periodic cell stripped ⇒ finite cluster

    def test_radius_crops_to_the_local_cluster(self, tmp_path):
        # a far-away stray O beyond `radius` must not enter the relaxed cluster
        atoms = _cluster(2.5)
        atoms += Atoms("O", positions=[[20, 0, 0]])
        src = str(tmp_path / "stray.pdb")
        atoms.write(src)
        site = CoordinationSite("Ni2+", np.zeros(3), atoms.get_positions()[1:7], ("N",) * 6)
        d = BinderDesign("S", site, generator="t", generator_confidence=0.5, source=src)
        r = MLIPVerifier(calculator=SpringToMetal(r0=2.1), radius=5.0).relax(d)
        assert r.cn_before == 6  # the stray O at 20 Å was cropped out


class TestMLIPDynamics:
    def test_bound_site_survives_md(self, tmp_path):
        np.random.seed(0)  # Langevin draws thermal forces from np.random
        v = MLIPDynamicsVerifier(
            calculator=SpringToMetal(r0=2.1, k=10.0), temperature=100.0, steps=200
        ).verify(_design(tmp_path, r=2.1))
        assert v.trust and not v.ood and v.score > 0.8  # shell holds through the run

    def test_unbound_site_defers(self, tmp_path):
        np.random.seed(0)
        v = MLIPDynamicsVerifier(
            calculator=SpringToMetal(r0=3.5, k=10.0), temperature=100.0, steps=200
        ).verify(_design(tmp_path, r=2.6))
        assert v.ood and not v.trust  # donors pulled past the cutoff ⇒ low retention

    def test_defers_when_md_fails(self, tmp_path):
        v = MLIPDynamicsVerifier(calculator=_Exploding()).verify(_design(tmp_path))
        assert v.ood and not v.trust and "failed" in v.reason

    def test_empty_initial_site_defers_not_trusts(self, tmp_path):
        # donors all start beyond the cutoff ⇒ CN=0; retention >= 0 is trivially true,
        # so this must defer, not be promoted to trust
        atoms = _cluster(r=4.0)  # 6 donors at 4.0 Å, all outside the 2.8 Å cutoff
        src = str(tmp_path / "empty.pdb")
        atoms.write(src)
        site = CoordinationSite("Ni2+", np.zeros(3), atoms.get_positions()[1:],
                                tuple(atoms.get_chemical_symbols()[1:]))
        d = BinderDesign("SEQ", site, generator="test", generator_confidence=0.5, source=src)
        v = MLIPDynamicsVerifier(calculator=SpringToMetal(r0=2.1, k=10.0), temperature=100.0, steps=40).verify(d)
        assert v.ood and not v.trust


class TestMLIPScoreBounds:
    def test_score_clamped_when_donors_migrate_in(self, tmp_path):
        # 1 donor in-shell + 5 just outside that the potential pulls in ⇒ cn_after > cn_before;
        # the cn_after/cn_before ratio must not push the score above 1
        pos = np.vstack([[0, 0, 0], _OCT[0] * 2.1, _OCT[1:] * 4.0])  # metal + 1 near + 5 far
        atoms = Atoms(symbols=["Ni", "N", "N", "O", "O", "N", "O"], positions=pos)
        src = str(tmp_path / "migrate.pdb")
        atoms.write(src)
        site = CoordinationSite("Ni2+", np.zeros(3), atoms.get_positions()[1:],
                                tuple(atoms.get_chemical_symbols()[1:]))
        d = BinderDesign("SEQ", site, generator="test", generator_confidence=0.5, source=src)
        # restrain=False so the 5 out-of-shell donors are free to migrate in (the case
        # the clamp guards); with the default restraint they would be frozen out-of-shell.
        v = MLIPVerifier(calculator=SpringToMetal(r0=2.1, k=10.0), restrain=False).verify(d)
        assert 0.0 <= v.score <= 1.0


def test_make_backbone_rejects_unknown():
    with pytest.raises(ValueError, match="unknown MLIP backbone"):
        make_backbone("not_a_model")
