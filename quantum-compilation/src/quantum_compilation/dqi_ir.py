"""A dqi circuit as numpy columns, and the queries that price holding a product on a spare wire.

The dqi circuit is 57M instructions, too many for one `submission.Op` per instruction to stay
alive. `CircuitColumns` streams the text form through `Op`'s tokenising once and keeps one row
per instruction -- opcode, up to three qubits, inline condition bit, condition-block depth --
cached beside the text as `.npz`.

A product is *recomputed* when a CCX repeats a control pair neither of whose qubits has been
written since that pair's previous CCX. The second becomes a free `CX` if the first is computed
into a clean wire that stays idle until then, and that wire is uncomputed for free by
`HMR(w) -> m` and `CZ(a, b) if m`, the measured uncomputation the public line already uses.
`hold_assignment` gives each recomputation a spare wire where one exists, without two holds
sharing a wire at once, and so bounds how many Toffolis the rewrite can reach.

Writes are counted conservatively (any op naming a qubit as its target, conditioned or not), and
only unconditioned CCX outside condition blocks are paired. A wire is clean when it is not a
syndrome qubit and was never touched, or its last touch was `HMR` or `R`.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import kickmix as km
import numpy as np

from quantum_compilation.submission import Op

OPS = (
    "APPEND_TO_REGISTER", "REGISTER", "X", "Z", "CX", "CZ", "SWAP", "CCX", "CCZ", "CSWAP", "HMR", "R",
    "NEG", "BIT_INVERT", "BIT_STORE0", "BIT_STORE1", "PUSH_CONDITION", "POP_CONDITION",
)
OPCODE = {name: code for code, name in enumerate(OPS)}
CLEANING = (OPCODE["HMR"], OPCODE["R"])
WRITES_TARGET = tuple(OPCODE[n] for n in ("X", "CX", "CCX"))  # the last qubit
WRITES_ALL = tuple(OPCODE[n] for n in ("SWAP", "HMR", "R"))  # every qubit named
WRITERS = (*WRITES_TARGET, *WRITES_ALL, OPCODE["CSWAP"])  # CSWAP: all but its control
SYNDROME_REGISTER = re.compile(r"S\d+")  # the verifier's S0..S{n-1}
NONE = -1


@dataclass(frozen=True)
class Recomputation:
    first: int  # instruction index of the earlier CCX on the pair
    second: int  # the CCX that recomputes it
    controls: tuple[int, int]

    @property
    def gap(self) -> int:
        return self.second - self.first


@dataclass
class CircuitColumns:
    op: np.ndarray  # uint8 opcode
    qubits: np.ndarray  # int16 (n, 3), NONE where absent
    cond: np.ndarray  # int32: inline `if` (or PUSH_CONDITION) bit
    depth: np.ndarray  # int16: condition blocks open after the op
    syndrome_qubits: np.ndarray  # qubits appended to the S registers

    @classmethod
    def from_lines(cls, lines: Iterator[str]) -> CircuitColumns:
        op: list[int] = []
        qubits: list[tuple[int, ...]] = []
        cond: list[int] = []
        registers: dict[str, str] = {}
        appended: list[tuple[int, str]] = []
        for line in lines:
            if not line.strip():
                continue
            name, *args = line.split()
            parsed = Op(name, tuple(args), 0)
            op.append(OPCODE[name])
            qubits.append(parsed.qubits)
            cond.append(int(bit[1:]) if (bit := parsed.condition_bit) else NONE)
            if name == "APPEND_TO_REGISTER" and parsed.qubits:
                appended.append((parsed.qubits[0], args[1]))
            elif name == "REGISTER":
                registers[args[0]] = line.split('"')[1]
        codes = np.array(op, np.uint8)
        table = np.full((len(qubits), 3), NONE, np.int16)
        for i, qs in enumerate(qubits):
            table[i, : len(qs)] = qs
        depth = np.cumsum(
            (codes == OPCODE["PUSH_CONDITION"]).astype(np.int16) - (codes == OPCODE["POP_CONDITION"]),
            dtype=np.int16,
        )
        syndrome = sorted(q for q, r in appended if SYNDROME_REGISTER.fullmatch(registers[r]))
        return cls(codes, table, np.array(cond, np.int32), depth, np.array(syndrome, np.int16))

    @classmethod
    def from_text(cls, text: str) -> CircuitColumns:
        return cls.from_lines(iter(text.splitlines()))

    @staticmethod
    def text_dump(circuit_path: Path, work_dir: Path) -> Path:
        """`work_dir/<stem>.kmx`, kickmix's text form of a `.kmb`, written on first use.

        It is 1 GB for dqi, so it goes in a working directory, never beside the circuit.
        """
        text_path = work_dir / f"{circuit_path.stem}.kmx"
        if not text_path.exists():
            work_dir.mkdir(parents=True, exist_ok=True)
            km.Circuit.from_file(str(circuit_path)).to_file(str(text_path))
        return text_path

    @classmethod
    def load(cls, text_path: Path) -> CircuitColumns:
        """The columns for a `.kmx` text dump, from its `.npz` cache when that is current."""
        cache = text_path.with_suffix(".npz")
        if cache.exists() and cache.stat().st_mtime >= text_path.stat().st_mtime:
            with np.load(cache, allow_pickle=False) as data:
                return cls(**{k: data[k] for k in data.files})
        with text_path.open() as fh:
            columns = cls.from_lines(fh)
        np.savez(cache, **columns.__dict__)
        return columns

    @cached_property
    def candidate_wires(self) -> dict[int, np.ndarray]:
        """Non-syndrome qubit -> sorted indices of every instruction naming it."""
        rows, slots = np.nonzero(self.qubits != NONE)
        wires = self.qubits[rows, slots]
        order = np.argsort(wires, kind="stable")  # rows are already ascending
        rows, wires = rows[order], wires[order]
        cuts = np.flatnonzero(np.diff(wires)) + 1
        syndrome = set(self.syndrome_qubits.tolist())
        return {
            int(w[0]): r
            for w, r in zip(np.split(wires, cuts), np.split(rows, cuts), strict=True)
            if int(w[0]) not in syndrome
        }

    def written(self, i: int) -> list[int]:
        qs = [q for q in self.qubits[i].tolist() if q != NONE]
        if self.op[i] in WRITES_ALL:
            return qs
        if self.op[i] == OPCODE["CSWAP"]:
            return qs[1:]
        return qs[-1:]

    def recomputations(self) -> list[Recomputation]:
        """Every CCX whose control pair is unchanged since that pair's previous CCX."""
        last_write: dict[int, int] = {}
        last_pair: dict[tuple[int, int], int] = {}
        found = []
        for i in np.flatnonzero(np.isin(self.op, WRITERS)).tolist():
            if self.op[i] == OPCODE["CCX"] and self.depth[i] == 0 and self.cond[i] == NONE:
                a, b = sorted(self.qubits[i, :2].tolist())
                prev = last_pair.get((a, b))
                if prev is not None and last_write.get(a, -1) < prev and last_write.get(b, -1) < prev:
                    found.append(Recomputation(prev, i, (a, b)))
                last_pair[(a, b)] = i
            for q in self.written(i):
                last_write[q] = i
        return found

    def idle_wires(self, start: int, stop: int) -> Iterator[int]:
        """Qubits clean before `start` and untouched over `[start, stop)`."""
        for wire, rows in self.candidate_wires.items():
            k = np.searchsorted(rows, start)
            if (k == len(rows) or rows[k] >= stop) and (k == 0 or self.op[rows[k - 1]] in CLEANING):
                yield wire

    def hold_assignment(self, found: list[Recomputation]) -> dict[Recomputation, int]:
        """A spare wire per recomputation where one exists, shortest gaps first, no overlaps."""
        busy: dict[int, list[tuple[int, int]]] = {}
        held = {}
        for r in sorted(found, key=lambda r: r.gap):
            for wire in self.idle_wires(r.first, r.second + 1):
                if all(r.second < s or r.first > e for s, e in busy.get(wire, [])):
                    busy.setdefault(wire, []).append((r.first, r.second))
                    held[r] = wire
                    break
        return held
