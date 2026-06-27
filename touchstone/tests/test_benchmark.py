"""End-to-end discrimination: the stack must separate a known-good metal site from a
weak design. The per-component tests prove each verifier; this proves the assembled
stack + reward actually rank good above bad (the thing all-weak fixtures couldn't show).
"""

from pathlib import Path

from touchstone.reward import rank_structures, reward_from_result
from touchstone.service import verify_structure

FIX = Path(__file__).parent / "fixtures"
IDEAL = FIX / "ideal_nickel_site.pdb"  # octahedral Ni²⁺ at the PDB mean bond (2.09 Å)
WEAK = FIX / "ligmpnn_nickel_packed.pdb"  # a strained CN-5 design


class TestDiscrimination:
    def test_ideal_site_trusts(self):
        r = verify_structure(IDEAL, "Ni2+")
        assert r["consensus"] == "trust"
        assert all(v["label"] == "trust" for v in r["verifiers"].values() if "label" in v)

    def test_reward_orders_good_above_weak(self):
        good = reward_from_result(verify_structure(IDEAL, "Ni2+"))
        weak = reward_from_result(verify_structure(WEAK, "Ni2+"))
        assert good > weak  # the assembled stack discriminates

    def test_rank_puts_the_good_site_first(self):
        ranked = rank_structures([WEAK, IDEAL], "Ni2+")
        assert Path(ranked[0]["structure"]).name == IDEAL.name
