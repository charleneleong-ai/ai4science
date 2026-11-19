import pytest

from quantum_compilation.dqi_ir import CircuitColumns, Recomputation

# q0, q1 control a product computed twice into q2 and q4; q5 is reset first, so it is clean.
BASE = ["R q5", "CCX q0 q1 q2", "CX q2 q3", "CCX q1 q0 q4"]


def columns(lines: list[str]) -> CircuitColumns:
    return CircuitColumns.from_text("\n".join(lines))


class TestRecomputations:
    def test_unchanged_controls_pair_the_two_ccx(self) -> None:
        assert columns(BASE).recomputations() == [Recomputation(1, 3, (0, 1))]

    @pytest.mark.parametrize("write", ["X q0", "CX q3 q1", "SWAP q0 q6", "HMR q1 b0", "CSWAP q6 q0 q7"])
    def test_a_write_to_a_control_breaks_the_pair(self, write: str) -> None:
        assert columns([*BASE[:3], write, BASE[3]]).recomputations() == []

    @pytest.mark.parametrize(
        "ccx", ["CCX q0 q1 q4 if b3", "PUSH_CONDITION if b3\nCCX q0 q1 q4\nPOP_CONDITION"]
    )
    def test_conditioned_ccx_is_not_paired(self, ccx: str) -> None:
        assert columns([*BASE[:3], *ccx.split("\n")]).recomputations() == []


class TestHoldAssignment:
    def test_clean_idle_wire_holds_the_product(self) -> None:
        circuit = columns(BASE)
        assert circuit.hold_assignment(circuit.recomputations()) == {Recomputation(1, 3, (0, 1)): 5}

    @pytest.mark.parametrize(
        "lines",
        [["R q5", "X q5", *BASE[1:]], [*BASE[:3], "CX q0 q5", BASE[3]]],
        ids=["dirty-at-start", "touched-during"],
    )
    def test_dirty_or_busy_wire_holds_nothing(self, lines: list[str]) -> None:
        circuit = columns(lines)
        assert circuit.hold_assignment(circuit.recomputations()) == {}

    @pytest.mark.parametrize(
        ("name", "held"), [('"S0"', {}), ('"SCRATCH"', {Recomputation(3, 5, (0, 1)): 5})]
    )
    def test_only_syndrome_registers_are_never_spare(
        self, name: str, held: dict[Recomputation, int]
    ) -> None:
        circuit = columns(["APPEND_TO_REGISTER q5 r0", f"REGISTER r0 {name}", *BASE])
        assert circuit.hold_assignment(circuit.recomputations()) == held

    def test_back_to_back_holds_reuse_one_wire(self) -> None:
        lines = [*BASE, "CCX q6 q7 q8", "CX q8 q9", "CCX q6 q7 q10"]
        circuit = columns(lines)
        held = circuit.hold_assignment(circuit.recomputations())
        assert sorted(held.values()) == [5, 5]
