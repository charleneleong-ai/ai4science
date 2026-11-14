"""A submission's files: the layout the leaderboard expects, and what the circuits contain.

The layout, from the Evaluation page: `invsqrt` and `dqi` as single circuits, and a
`dfthcsos/` folder holding `angle_k_0` .. `angle_k_<n-1>` numbered from zero with no gaps
plus `submission.npz` carrying the factorization and the Givens tree. `SubmissionAudit`
reads only names and the npz index -- it parses no circuits, so auditing a 905 MB submission
is instant, which is what lets the gate run before anything expensive does. `Resources` is
what a parse buys once one is worth paying for: the measured width and non-Clifford content
of a circuit, taken from kickmix's own accessors so that kickmix owns the semantics -- the
Clifford-hierarchy level of a `Z_POW`, and what counts as Toffoli-like or T-like.

Nothing here knows the competition's cost model, not even the 1/4 that turns T gates into
Toffoli-equivalents -- that lives in `cost`, which consumes these counts. Nor can the counts
be turned into a cost: the scored Toffoli figure is a simulated per-shot mean, so what is
measured here bounds a candidate's cost from above and never gives it.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import kickmix as km
import numpy as np

from .cost import require_task

CIRCUIT_SUFFIXES = (".kmb", ".kmx")
ANGLE_PATTERN = re.compile(r"^angle_k_(\d+)$")
CHEMISTRY_DIR = "dfthcsos"
FACTORIZATION_FILE = "submission.npz"
# The arrays submission.npz must carry: the factorization and the Givens rotation tree.
NPZ_ARRAYS = ("U", "W", "H2_I", "H1_I", "Hsym_lower", "tree_edges")
# What reading the factorization *without* building its circuits additionally needs. `check`
# does not require these, because the verifier recovers the words by simulating the ANGLE
# circuits rather than reading them; `chemistry.ChemistryScreen` reads them directly, which is
# the whole of why it is fast and why it is a screen rather than a verification.
NPZ_SCREEN_ARRAYS = ("angle_words", "brot")

# A `Z_POW` past the T level is one Clifford+T cannot express exactly, so it needs synthesis.
# The level marks *what needs* synthesising; it does not price it -- the metric charges every
# such rotation at one circuit-wide synthesis error, regardless of its level.
T_LEVEL = 3


@dataclass(frozen=True)
class Resources:
    """Measured width and non-Clifford content of one circuit.

    Deliberately small and hashable: a search loop keeps thousands of these while the
    parsed `km.Circuit` each came from -- up to 1.46 GB for a public dqi decoder -- is
    discarded.
    """

    num_qubits: int
    num_toffoli_like: int
    num_t_like: int
    num_fine_rotations: int

    @classmethod
    def from_circuit(cls, circuit: km.Circuit) -> Resources:
        stats = circuit.operation_statistics()
        by_level = stats.z_pow_clifford_hierarchy_counts
        return cls(
            num_qubits=circuit.num_qubits,
            num_toffoli_like=stats.num_toffoli_like,
            num_t_like=stats.num_t_like,
            num_fine_rotations=sum(n for lvl, n in by_level.items() if lvl > T_LEVEL),
        )

    @classmethod
    def from_file(cls, path: Path | str) -> Resources:
        return cls.from_circuit(km.Circuit.from_file(path))


@dataclass(frozen=True)
class Op:
    """One instruction of a circuit's text form, the only op-level view kickmix exposes."""

    name: str
    args: tuple[str, ...]
    depth: int  # condition blocks open after this op

    DECLARATIONS = ("APPEND_TO_REGISTER", "REGISTER")

    @property
    def conditioned(self) -> bool:
        """Inside a `PUSH_CONDITION` block (its own PUSH and POP included) or carrying an `if`."""
        return self.depth > 0 or self.name == "POP_CONDITION" or "if" in self.args

    @property
    def qubits(self) -> tuple[int, ...]:
        return tuple(int(a[1:]) for a in self.args if a.startswith("q"))

    @property
    def is_fine_rotation(self) -> bool:
        """A `Z_POW` whose exponent is not a multiple of T's (Z^0.25), so it needs synthesis."""
        return self.name == "Z_POW" and not (4 * float(self.args[1])).is_integer()

    @property
    def text(self) -> str:
        return " ".join((self.name, *self.args))

    @staticmethod
    def parse(circuit: km.Circuit) -> list[Op]:
        """Every instruction, indexed like the lines of `str(circuit)`."""
        ops, depth = [], 0
        for line in str(circuit).splitlines():
            if not line.strip():
                continue
            name, *args = line.split()
            depth += (name == "PUSH_CONDITION") - (name == "POP_CONDITION")
            ops.append(Op(name, tuple(args), depth))
        return ops

    @staticmethod
    def circuit(ops: list[Op]) -> km.Circuit:
        return km.Circuit("\n".join(op.text for op in ops))

    @staticmethod
    def body_start(ops: list[Op]) -> int:
        """The first instruction after the register declarations, which come first."""
        return sum(op.name in Op.DECLARATIONS for op in ops)

    @staticmethod
    def cut_after(ops: list[Op], at: int) -> int:
        """The first cut at or after `at` that neither splits a condition block nor falls among
        the declarations."""
        at = max(at, Op.body_start(ops))
        while at > 0 and ops[at - 1].depth:
            at += 1
        return at

    @staticmethod
    def steps(ops: list[Op], start: int, stop: int) -> list[tuple[int, int]]:
        """`[first, last)` op ranges covering `start..stop`: one op, or one whole condition block.
        `start` must not fall inside a block."""
        steps, at = [], start
        while at < stop:
            end = at + 1
            while ops[end - 1].depth:
                end += 1
            steps.append((at, end))
            at = end
        return steps


@dataclass(frozen=True)
class SubmissionAudit:
    """What a submission directory holds, and every way its layout is wrong."""

    root: Path
    invsqrt: Path | None
    dqi: Path | None
    angle_circuits: tuple[Path, ...]
    problems: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not self.problems

    @property
    def circuits(self) -> dict[str, Path]:
        """Every circuit the audit found, keyed by the name it was found under."""
        single = {"invsqrt": self.invsqrt, "dqi": self.dqi}
        return {name: p for name, p in single.items() if p is not None} | {
            path.stem: path for path in self.angle_circuits
        }

    @property
    def chemistry(self) -> Path | None:
        """The `dfthcsos/` directory, which is what that task's verifier is pointed at."""
        chemistry = self.root / CHEMISTRY_DIR
        return chemistry if chemistry.is_dir() else None

    def target(self, task: str) -> Path | None:
        """What `task`'s verifier takes: a circuit for two of them, a directory for the third.

        Layout knowledge stays in this module, so a driver never has to rebuild a path.
        """
        if task == "dfthcsos":
            return self.chemistry
        return self.circuits.get(require_task(task))

    @staticmethod
    def task_circuit(directory: Path, stem: str) -> Path | None:
        """The `.kmb` or `.kmx` circuit named `stem` in `directory`, if exactly one exists."""
        found = [p for s in CIRCUIT_SUFFIXES if (p := directory / f"{stem}{s}").is_file()]
        return found[0] if len(found) == 1 else None

    @staticmethod
    def indexed_angle_circuits(directory: Path) -> dict[int, Path]:
        circuits: dict[int, Path] = {}
        for path in sorted(directory.iterdir()) if directory.is_dir() else []:
            match = ANGLE_PATTERN.match(path.stem)
            if match and path.suffix in CIRCUIT_SUFFIXES:
                circuits[int(match.group(1))] = path
        return circuits

    @staticmethod
    def npz_problems(path: Path, extra: Sequence[str] = ()) -> list[str]:
        """Why `path` is not a readable factorization, checking `extra` arrays beyond the
        required ones so a caller with stricter needs does not hand-roll a second reader."""
        if not path.is_file():
            return [f"missing {CHEMISTRY_DIR}/{FACTORIZATION_FILE}"]
        try:
            # allow_pickle stays off: a submission npz is data, never code.
            with np.load(path, allow_pickle=False) as data:
                present = set(data.files)
        except (OSError, ValueError) as exc:
            return [f"{path.name} could not be read as an npz: {exc}"]
        if absent := [name for name in (*NPZ_ARRAYS, *extra) if name not in present]:
            return [f"{path.name} is missing array(s) {absent}"]
        return []

    @classmethod
    def from_root(cls, root: Path | str) -> SubmissionAudit:
        root = Path(root)
        if not root.is_dir():
            return cls(root, None, None, (), (f"{root} is not a directory",))

        single = {stem: cls.task_circuit(root, stem) for stem in ("invsqrt", "dqi")}
        problems = [
            f"missing {stem}{' or '.join(CIRCUIT_SUFFIXES)}"
            for stem, path in single.items()
            if path is None
        ]

        chemistry = root / CHEMISTRY_DIR
        indexed = cls.indexed_angle_circuits(chemistry)
        if not chemistry.is_dir():
            problems.append(f"missing {CHEMISTRY_DIR}/ directory")
        elif not indexed:
            problems.append(f"{CHEMISTRY_DIR}/ holds no angle_k_<n> circuits")
        elif missing := sorted(set(range(max(indexed) + 1)) - set(indexed)):
            problems.append(
                f"{CHEMISTRY_DIR}/ angle circuits must run from 0 with no gaps; "
                f"missing index(es) {missing}"
            )
        problems += cls.npz_problems(chemistry / FACTORIZATION_FILE)

        return cls(
            root=root,
            invsqrt=single["invsqrt"],
            dqi=single["dqi"],
            angle_circuits=tuple(indexed[i] for i in sorted(indexed)),
            problems=tuple(problems),
        )
