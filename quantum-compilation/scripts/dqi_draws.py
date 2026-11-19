"""Re-verify one dqi circuit under many salts and keep every failed shot's input: gate G2's draws.

The official `verify()` reports only how many shots failed, and its batch loop is a closure that
cannot be imported, so `run_batch` and `draw` are a copy of it: the seed derivation, simulator
setup, error masks, the 1/16 rejection and the score. The helpers it calls (`sample_shots`,
`validate_registers`, `repetition_overhead`, ...) are the shipped module's own. Re-check the copy
whenever the verifier changes. The copy additionally keeps, for each failed shot, the error word
it was asked to decode: positions, values, the MBUC words at those positions, and which check
failed. The error word is recovered by replaying the batch's RNG in
`sample_shots`' own draw order, and each replayed word is checked against the verifier's
syndromes, so a replay that drifted from the verifier cannot produce a wrong record.

The empty salt reproduces the verifier's default draw, so its cost must equal the official
`verify()` cost exactly; `--check-default` asserts that before any other draw is trusted.

    uv run python -I scripts/dqi_draws.py <circuit.kmb> <out.jsonl> --draws 48 \\
        --verifiers <verifiers/>
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Annotated, Any

import kickmix as km
import numpy as np
import typer

from quantum_compilation.verify import load_verifier

SALT_PREFIX = "g2-draw-"


@dataclass(frozen=True)
class FailedShot:
    shot: int
    kind: str  # "value" (registers or ancillas not restored) or "phase"
    positions: list[int]
    values: list[int]
    mbuc_at_positions: list[int]
    syndrome_zero_words: list[int]  # indices i with S_i == 0


@dataclass
class Batch:
    sim: km.Simulator
    syndromes: np.ndarray
    mbuc: np.ndarray
    phases: np.ndarray
    inputs: np.ndarray  # the syndromes as the bits written to the S registers


@dataclass
class DrawRunner:
    """One circuit, one verifier, any number of salts."""

    dqi: ModuleType
    circuit: km.Circuit
    max_workers: int

    def __post_init__(self) -> None:
        self.spec = self.dqi.ProblemSpecDQI.benchmark_instance()
        self.tables = self.dqi.FieldTables(km.GF2Field(self.spec.degree))
        self.weight_probs = (
            self.dqi.optimal_dqi_weights(self.spec.m, self.spec.max_errors, self.spec.degree, self.spec.r)
            ** 2
        )
        syndromes, mbuc, self.ancillas = self.dqi.validate_registers(self.circuit, self.spec)
        self.syndrome_array, self.mbuc_array = km.array(syndromes), km.array(mbuc)
        self.circuit_hash = self.dqi._circuit_sha256(self.circuit)

    def error_words(self, rng_seed: int, batch: int) -> dict[int, tuple[np.ndarray, np.ndarray]]:
        """`shot -> (positions, values)`, replaying `sample_shots`' draws in its order."""
        rng = np.random.default_rng(rng_seed)
        max_value = (1 << self.tables.degree) - 1
        weights = rng.choice(len(self.weight_probs), size=batch, p=self.weight_probs)
        rng.integers(0, max_value, size=(batch, self.spec.m), endpoint=True, dtype=np.int64)
        words = {}
        for shot in np.flatnonzero(weights):
            pos = rng.choice(self.spec.m, size=int(weights[shot]), replace=False)
            vals = rng.integers(1, max_value, size=len(pos), endpoint=True, dtype=np.int64)
            words[int(shot)] = (pos, vals)
        return words

    def syndromes_of(self, pos: np.ndarray, vals: np.ndarray) -> np.ndarray:
        rows = np.arange(1, self.spec.n + 1, dtype=np.int64)[:, None]
        exps = (self.tables.log[vals][None, :] + rows * pos[None, :]) % self.tables.order
        return np.bitwise_xor.reduce(self.tables.antilog[exps], axis=1)

    def failures(self, start: int, rng_seed: int, syndromes: np.ndarray, mbuc: np.ndarray,
                 val_err: np.ndarray, phase_err: np.ndarray) -> list[FailedShot]:
        failed = np.flatnonzero(val_err | phase_err)
        if not len(failed):
            return []
        words = self.error_words(rng_seed, len(val_err))
        out = []
        for shot in failed.tolist():
            pos, vals = words.get(shot, (np.zeros(0, np.int64), np.zeros(0, np.int64)))
            if not np.array_equal(self.syndromes_of(pos, vals), syndromes[shot]):
                raise AssertionError(f"replayed error word for shot {start + shot} disagrees")
            out.append(FailedShot(
                shot=start + shot,
                kind="value" if val_err[shot] else "phase",
                positions=pos.tolist(),
                values=vals.tolist(),
                mbuc_at_positions=mbuc[shot, pos].tolist(),
                syndrome_zero_words=np.flatnonzero(syndromes[shot] == 0).tolist(),
            ))
        return out

    def start_batch(self, task: tuple[int, int, int]) -> Batch:
        """A simulator loaded with one batch's inputs, before any instruction runs."""
        _, rng_seed, sim_seed = task
        batch = self.dqi.DEFAULT_BATCH_SIZE
        rng = np.random.default_rng(rng_seed)
        syndromes, mbuc, phases = self.dqi.sample_shots(self.spec, self.tables, rng, batch, self.weight_probs)
        sim = km.Simulator(batch_size=batch, debug_print=False, counting_operations=True, seed=sim_seed)
        sim.use_same_registers_as(self.circuit)
        inputs = self.dqi.to_bits(syndromes, self.spec.degree)
        sim.write_across_shots(self.syndrome_array, inputs)
        sim.write_across_shots(self.mbuc_array, self.dqi.to_bits(mbuc, self.spec.degree))
        return Batch(sim, syndromes, mbuc, phases, inputs)

    def errors(self, b: Batch) -> tuple[np.ndarray, np.ndarray]:
        """`(value error, phase error)` per shot, by the verifier's rule, once the circuit has run."""
        batch = b.sim.batch_size
        junk = b.sim.read_across_shots(self.ancillas) if len(self.ancillas) else np.zeros((0, batch), bool)
        val_err = (b.sim.read_across_shots(self.syndrome_array) ^ b.inputs).any(axis=0) | junk.any(axis=0)
        return val_err, (b.sim.read_phase_flipped_across_shots() ^ b.phases) & ~val_err

    def run_batch(self, task: tuple[int, int, int]) -> tuple[int, list[FailedShot]]:
        b = self.start_batch(task)
        b.sim.do(self.circuit)
        val_err, phase_err = self.errors(b)
        toffolis = int(b.sim.operation_statistics().num_toffoli_like)
        return toffolis, self.failures(task[0], task[1], b.syndromes, b.mbuc, val_err, phase_err)

    def tasks(self, salt: str) -> list[tuple[int, int, int]]:
        """`(first shot, rng seed, simulator seed)` per batch, as `verify()` derives them."""
        starts = range(0, self.spec.num_shots, self.dqi.DEFAULT_BATCH_SIZE)
        base_seed = self.dqi._initial_seed(self.circuit_hash, salt)
        seeds = np.random.SeedSequence(base_seed, pool_size=8).generate_state(
            8 * len(starts), dtype=np.uint64
        ).reshape(-1, 2, 4)
        return [
            (s, int.from_bytes(seeds[i, 0], "little"), int.from_bytes(seeds[i, 1], "little"))
            for i, s in enumerate(starts)
        ]

    def draw(self, salt: str) -> dict[str, Any]:
        num_shots = self.spec.num_shots
        results = self.dqi._run_parallel(self.run_batch, self.tasks(salt), self.max_workers)
        failed = [f for _, batch_failures in results for f in batch_failures]
        eps_plus = self.dqi.clopper_pearson_upper_bound(len(failed), num_shots,
                                                        log2_delta=self.dqi.DEFAULT_LOG2_DELTA)
        if eps_plus >= 1 / 16:
            raise AssertionError(f"salt {salt!r}: effective error rate {eps_plus:.4f} >= 1/16")
        toffolis = sum(t for t, _ in results) / num_shots
        cost = self.dqi.repetition_overhead(eps_plus) * self.circuit.num_qubits ** 2 * toffolis
        return {
            "salt": salt,
            "circuit_sha256": self.circuit_hash,
            "shots": num_shots,
            "failures": len(failed),
            "toffoli_count": toffolis,
            "cost": cost,
            "failed_shots": [f.__dict__ for f in failed],
        }


def main(
    circuit_path: Path,
    out: Path,
    verifiers: Annotated[Path, typer.Option(envvar="VQCC_VERIFIERS")],
    draws: int = 48,
    first: int = 0,
    max_workers: int = 8,
    check_default: Annotated[float | None, typer.Option(help="official default-salt cost")] = None,
) -> None:
    runner = DrawRunner(load_verifier(verifiers, "dqi"), km.Circuit.from_file(str(circuit_path)), max_workers)
    if check_default is not None:
        default = runner.draw("")
        if default["cost"] != check_default:
            raise typer.Exit(f"default salt cost {default['cost']!r} != official {check_default!r}")
        typer.echo(f"default salt reproduces the official cost {check_default!r}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("a") as fh:
        for index in range(first, first + draws):
            t0 = time.time()
            row = runner.draw(f"{SALT_PREFIX}{index}")
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            stamp = time.strftime("%H:%M:%S")
            typer.echo(f"[{stamp}] draw {index}: {row['failures']} failed, cost {row['cost']:.6e} "
                       f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    typer.run(main)
