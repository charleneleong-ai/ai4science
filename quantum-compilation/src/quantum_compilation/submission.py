"""Structural audit of a submission directory.

The layout the leaderboard expects, from the Evaluation page: `invsqrt` and `dqi` as single
circuits, and a `dfthcsos/` folder holding `angle_k_0` .. `angle_k_<n-1>` numbered from zero
with no gaps plus `submission.npz` carrying the factorization and the Givens tree.

This module only looks at names and the npz index -- it parses no circuits, so auditing a
905 MB submission is instant. Measuring the circuits is `resources.Resources`; pricing them
is `cost`. Keeping the three apart is what lets the gate run before anything expensive does.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .cost import require_task

CIRCUIT_SUFFIXES = (".kmb", ".kmx")
ANGLE_PATTERN = re.compile(r"^angle_k_(\d+)$")
CHEMISTRY_DIR = "dfthcsos"
FACTORIZATION_FILE = "submission.npz"
# The arrays submission.npz must carry: the factorization and the Givens rotation tree.
NPZ_ARRAYS = ("U", "W", "H2_I", "H1_I", "Hsym_lower", "tree_edges")


def task_circuit(directory: Path, stem: str) -> Path | None:
    """The `.kmb` or `.kmx` circuit named `stem` in `directory`, if exactly one exists."""
    found = [directory / f"{stem}{suffix}" for suffix in CIRCUIT_SUFFIXES]
    present = [path for path in found if path.is_file()]
    return present[0] if len(present) == 1 else None


def indexed_angle_circuits(directory: Path) -> dict[int, Path]:
    circuits: dict[int, Path] = {}
    for path in sorted(directory.iterdir()) if directory.is_dir() else []:
        match = ANGLE_PATTERN.match(path.stem)
        if match and path.suffix in CIRCUIT_SUFFIXES:
            circuits[int(match.group(1))] = path
    return circuits


def npz_problems(path: Path) -> list[str]:
    if not path.is_file():
        return [f"missing {CHEMISTRY_DIR}/{FACTORIZATION_FILE}"]
    try:
        # allow_pickle stays off: a submission npz is data, never code.
        with np.load(path, allow_pickle=False) as data:
            present = set(data.files)
    except (OSError, ValueError) as exc:
        return [f"{path.name} could not be read as an npz: {exc}"]
    if absent := [name for name in NPZ_ARRAYS if name not in present]:
        return [f"{path.name} is missing array(s) {absent}"]
    return []


@dataclass(frozen=True)
class SubmissionAudit:
    root: Path
    invsqrt: Path | None
    dqi: Path | None
    angle_circuits: tuple[Path, ...]
    factorization: Path | None
    problems: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not self.problems

    @property
    def circuits(self) -> dict[str, Path]:
        """Every circuit the audit found, keyed by the name it was found under."""
        named = {
            name: path
            for name, path in (("invsqrt", self.invsqrt), ("dqi", self.dqi))
            if path is not None
        }
        return named | {path.stem: path for path in self.angle_circuits}

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


def audit(root: Path | str) -> SubmissionAudit:
    root = Path(root)
    problems: list[str] = []
    if not root.is_dir():
        return SubmissionAudit(
            root, None, None, (), None, (f"{root} is not a directory",)
        )

    single = {stem: task_circuit(root, stem) for stem in ("invsqrt", "dqi")}
    problems += [
        f"missing {stem}{' or '.join(CIRCUIT_SUFFIXES)}"
        for stem, path in single.items()
        if path is None
    ]

    chemistry = root / CHEMISTRY_DIR
    indexed = indexed_angle_circuits(chemistry)
    if not chemistry.is_dir():
        problems.append(f"missing {CHEMISTRY_DIR}/ directory")
    elif not indexed:
        problems.append(f"{CHEMISTRY_DIR}/ holds no angle_k_<n> circuits")
    elif missing := sorted(set(range(max(indexed) + 1)) - set(indexed)):
        problems.append(
            f"{CHEMISTRY_DIR}/ angle circuits must run from 0 with no gaps; "
            f"missing index(es) {missing}"
        )

    factorization = chemistry / FACTORIZATION_FILE
    problems += npz_problems(factorization)

    return SubmissionAudit(
        root=root,
        invsqrt=single["invsqrt"],
        dqi=single["dqi"],
        angle_circuits=tuple(indexed[i] for i in sorted(indexed)),
        factorization=factorization if factorization.is_file() else None,
        problems=tuple(problems),
    )
