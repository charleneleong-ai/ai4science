"""Where a lost dqi shot's decoding goes wrong: rare flags, unrestored steps and zero words.

Re-runs each batch that lost a shot in `dqi_draws.py`'s records, with the verifier's own seeds,
one decoder stage at a time (stages from a public `vqcc-regions-sidecar` map), and compares each
lost shot with the batch's passing shots. It reports, per lost shot:

- `rare_bits`: classical bits whose value at a stage's end at most `--rare` passing shots share,
  a branch the decoder almost never takes (measurement outcomes are near 50/50, never rare);
- `unrestored`: each `*_prefix*` stage whose mirrored `*_return*` stage leaves a qubit that every
  passing shot restores, outermost first, so the last entry is the innermost step that failed;
- `zero_words`: aligned field-width words of qubits that are zero at a stage's start in the lost
  shot and in at most `--rare` passing shots.

    uv run python scripts/dqi_flags.py <circuit.kmb> <work-dir> <regions.json> \\
        <draws.jsonl> <out.jsonl> --verifiers <verifiers/>

(not `-I`: it imports `dqi_draws` from beside it). Boundaries inside a condition block are skipped,
merging the stages either side of them.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Any, NamedTuple

import kickmix as km
import numpy as np
import typer
from dqi_draws import DrawRunner

from quantum_compilation.dqi_ir import OPCODE, CircuitColumns
from quantum_compilation.submission import Op
from quantum_compilation.verify import load_verifier

Job = tuple[tuple[int, int, int], tuple[int, ...]]  # a batch's verifier task, and its lost shot offsets


class Stage(NamedTuple):
    name: str
    circuit: km.Circuit


def stages(text_path: Path, regions_path: Path) -> list[Stage]:
    """The circuit cut at every stage boundary that is outside a condition block."""
    columns = CircuitColumns.load(text_path)
    offset = int(np.isin(columns.op, [OPCODE[n] for n in Op.DECLARATIONS]).sum())  # they lead the text
    cuts = sorted({
        (offset + r["end"], r["name"])
        for r in json.loads(regions_path.read_text())["regions"]
        if r["layer"] == "stage" and columns.depth[offset + r["end"] - 1] == 0
    }) + [(len(columns.op), "(tail)")]
    out, lines, cut = [], [], iter(cuts)
    end, name = next(cut)
    with text_path.open() as fh:  # streamed: the text is 1 GB
        for i, line in enumerate(line for line in fh if line.strip()):
            if i == end:
                out.append(Stage(name, km.Circuit("".join(lines))))
                lines = []
                end, name = next(cut)
            lines.append(line)
    out.append(Stage(name, km.Circuit("".join(lines))))
    return out


@dataclass
class FlagProbe:
    runner: DrawRunner
    stages: list[Stage]
    rare: int
    pairs: list[tuple[int, int]] = field(init=False)  # (prefix stage, its return stage), outermost first

    def __post_init__(self) -> None:
        index = {stage.name: i for i, stage in enumerate(self.stages)}
        found = [(i, index[s.name.replace("_prefix", "_return")]) for i, s in enumerate(self.stages)
                 if "_prefix" in s.name and s.name.replace("_prefix", "_return") in index]
        self.pairs = sorted(found, key=lambda p: p[0] - p[1])  # widest span first
        self.bits = km.array([km.b(i) for i in range(self.runner.circuit.num_bits)])
        self.qubits = km.array([km.q(i) for i in range(self.runner.circuit.num_qubits)])
        self.word = self.runner.spec.degree

    def run(self, job: Job) -> list[dict[str, Any]]:
        task, lost = job
        b = self.runner.start_batch(task)
        passing = np.setdiff1d(np.arange(b.sim.batch_size), lost)
        rare: dict[int, dict[int, dict[str, Any]]] = {s: {} for s in lost}
        snapshots = [b.sim.read_across_shots(self.qubits)]  # snapshots[i]: qubits before stage i
        for index, stage in enumerate(self.stages):
            b.sim.do(stage.circuit)
            snapshots.append(b.sim.read_across_shots(self.qubits))
            self.note_rare_bits(rare, index, b.sim.read_across_shots(self.bits), passing)
        val_err, phase_err = self.runner.errors(b)
        failed = set(np.flatnonzero(val_err | phase_err).tolist())
        if failed != set(lost):
            raise AssertionError(f"batch {task[0]}: lost {sorted(failed)}, records say {sorted(lost)}")
        return [
            {"shot": task[0] + s, "passing": len(passing), "rare_bits": list(rare[s].values()),
             "unrestored": self.unrestored(snapshots, passing, s),
             "zero_words": self.zero_words(snapshots, passing, s)}
            for s in rare
        ]

    def note_rare_bits(self, rare: dict[int, dict[int, dict[str, Any]]], index: int,
                       values: np.ndarray, passing: np.ndarray) -> None:
        ones = values[:, passing].sum(axis=1)
        for s, seen in rare.items():
            same = np.where(values[:, s], ones, len(passing) - ones)
            for bit in np.flatnonzero(same <= self.rare).tolist():
                seen.setdefault(bit, {"stage": index, "name": self.stages[index].name, "bit": bit,
                                      "value": bool(values[bit, s]), "passing_alike": int(same[bit])})

    def unrestored(self, snapshots: list[np.ndarray], passing: np.ndarray, s: int) -> list[dict[str, Any]]:
        out = []
        for forward, back in self.pairs:
            before, after = snapshots[forward], snapshots[back + 1]
            restored = (before[:, passing] == after[:, passing]).all(axis=1)
            if wrong := np.flatnonzero(restored & (before[:, s] != after[:, s])).tolist():
                out.append({"forward": self.stages[forward].name, "return": self.stages[back].name,
                            "qubits": wrong})
        return out

    def zero_words(self, snapshots: list[np.ndarray], passing: np.ndarray, s: int) -> list[dict[str, Any]]:
        found: dict[int, dict[str, Any]] = {}
        for index, state in enumerate(snapshots[:-1]):
            words = state[: self.word * (len(state) // self.word)].reshape(-1, self.word, state.shape[1])
            zero = ~words.any(axis=1)  # (word, shot)
            for word in np.flatnonzero(zero[:, s] & (zero[:, passing].sum(axis=1) <= self.rare)).tolist():
                found.setdefault(word, {"stage": index, "name": self.stages[index].name, "word": word,
                                        "qubits": [self.word * word, self.word * (word + 1) - 1],
                                        "passing_zero": int(zero[word, passing].sum())})
        return list(found.values())


def failed_batches(runner: DrawRunner, draws_path: Path) -> list[Job]:
    batch = runner.dqi.DEFAULT_BATCH_SIZE
    lost: dict[tuple[int, int, int], set[int]] = {}
    for line in draws_path.read_text().splitlines():
        row = json.loads(line)
        tasks = runner.tasks(row["salt"])
        for f in row["failed_shots"]:
            lost.setdefault(tasks[f["shot"] // batch], set()).add(f["shot"] % batch)
    return [(task, tuple(sorted(shots))) for task, shots in sorted(lost.items())]


def main(
    circuit_path: Path,
    work_dir: Path,
    regions_path: Path,
    draws_path: Path,
    out: Path,
    verifiers: Annotated[Path, typer.Option(envvar="VQCC_VERIFIERS")],
    rare: int = 2,
    max_workers: int = 8,
) -> None:
    runner = DrawRunner(load_verifier(verifiers, "dqi"), km.Circuit.from_file(str(circuit_path)), max_workers)
    probe = FlagProbe(runner, stages(CircuitColumns.text_dump(circuit_path, work_dir), regions_path), rare)
    tasks = failed_batches(runner, draws_path)
    typer.echo(f"{len(probe.stages)} stages, {len(probe.pairs)} prefix/return pairs; "
               f"re-running {len(tasks)} batches")
    shots = [s for found in runner.dqi._run_parallel(probe.run, tasks, max_workers) for s in found]
    out.write_text("".join(json.dumps(s) + "\n" for s in shots))
    innermost = Counter(s["unrestored"][-1]["forward"] if s["unrestored"] else "(none)" for s in shots)
    typer.echo(f"{len(shots)} lost shots; innermost unrestored forward stage: {dict(innermost)}")


if __name__ == "__main__":
    typer.run(main)
