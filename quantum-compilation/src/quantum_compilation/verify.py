"""Driving the three official verifiers -- the only scorer whose numbers are the real ones.

This competition has no train/dev split: there is one hidden test, re-drawn each time the
leaderboard scores a submission. The dev set is therefore the **official verifier run
locally**: the same code the leaderboard runs, over the same circuits. Submitting to learn a
score is a 5-a-day ration; running this is free and unlimited, so every candidate goes
through here first and only a winner costs a submission.

**Comparing candidates needs an explicit `--seed`.** Every verifier derives its simulator
seed from the *candidate's own bytes* -- `int(circuit_sha256, 16)`, or
`sha256(circuit_sha256 + salt)` when a salt is given -- so two different circuits draw two
different seeds no matter what salt is set, and a re-run of one circuit is bit-identical.
Left alone, that puts an unmatched draw inside every local comparison and makes its noise
simultaneously present and unmeasurable. All three verifiers accept `--seed`, which
overrides the hash entirely; pass one (or several) to hold the draw fixed across candidates.
The salt only reproduces the published default derivation.

All three verifiers share one contract, which is what lets a single driver cover them: a
JSON object on stdout with `status` of PASS or FAIL, `score` on a pass, `error` on a fail,
and a non-zero exit on a fail. This module owns the subprocess and that JSON, and nothing
else: the cost model, including which fields it reads, lives in `cost`.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .cost import (
    MODEL_AGREEMENT_REL,
    TASKS,
    MetricCalibration,
    cost_from_report,
    require_task,
)
from .submission import SubmissionAudit

# Where each verifier lives under the extracted `verifiers/` tree.
VERIFIER_SCRIPTS = {
    "invsqrt": "invsqrt/verify_submission_invsqrt.py",
    "dqi": "dqi/verify_submission_dqi.py",
    "dfthcsos": "dfthcsos/verify.py",
}

# What each verifier calls the hash of the thing it scored. Recorded with every verdict, so
# a cost can always be traced back to the exact bytes that produced it.
REPORT_HASH_KEYS = ("circuit_sha256", "submission_sha256")

# The empty salt reproduces the verifiers' published default seed derivation.
DEFAULT_SALT = ""


@dataclass(frozen=True)
class Verification:
    """One verifier's verdict on one task, with the invocation that produced it.

    The invocation is part of the result: a cost is a simulated mean at a particular seed
    and shot budget, so a verdict that does not record them cannot be compared with another
    one later.
    """

    task: str
    status: str
    score: float | None
    seconds: float
    salt: str = DEFAULT_SALT
    seed: int | None = None
    detail: dict[str, Any] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return self.status == "PASS" and self.score is not None

    @property
    def error(self) -> str | None:
        return None if self.passed else str(self.detail.get("error", self.status))

    @property
    def candidate_sha256(self) -> str | None:
        """The hash of what was scored, under whichever name this verifier reports it."""
        return next(
            (self.detail[key] for key in REPORT_HASH_KEYS if key in self.detail), None
        )

    def row(self) -> dict[str, Any]:
        """A flat record for a results log: the verdict, its invocation and the model's cost."""
        return {
            "task": self.task,
            "status": self.status,
            "score": self.score,
            "model_cost": recompute_cost(self),
            "seconds": self.seconds,
            "salt": self.salt,
            "seed": self.seed,
            "candidate_sha256": self.candidate_sha256,
            "detail": self.detail,
        }


def verifier_script(verifiers_root: Path | str, task: str) -> Path:
    script = Path(verifiers_root) / VERIFIER_SCRIPTS[require_task(task)]
    if not script.is_file():
        raise FileNotFoundError(
            f"no {task} verifier at {script}; point --verifiers at the extracted "
            "`verifiers/` directory from the competition download"
        )
    return script


def parse_verifier_output(
    task: str, stdout: str, seconds: float, **provenance: Any
) -> Verification:
    """The verifiers' shared JSON contract, with their own `status` taken at face value."""
    try:
        detail = json.loads(stdout)
    except json.JSONDecodeError:
        tail = stdout.strip()[-400:] or "(no output)"
        return Verification(task, "UNPARSEABLE", None, seconds, detail={"error": tail}, **provenance)
    if not isinstance(detail, dict):
        return Verification(
            task, "UNPARSEABLE", None, seconds, detail={"error": repr(detail)}, **provenance
        )
    score = detail.get("score")
    return Verification(
        task=task,
        status=str(detail.get("status", "UNKNOWN")),
        score=float(score) if isinstance(score, (int, float)) else None,
        seconds=seconds,
        detail=detail,
        **provenance,
    )


def run_verifier(
    verifiers_root: Path | str,
    task: str,
    target: Path,
    *,
    salt: str = DEFAULT_SALT,
    seed: int | None = None,
    timeout: float | None = None,
) -> Verification:
    """Run one verifier to completion and return its verdict.

    `-I` keeps the interpreter from putting the target's own directory on `sys.path`: a
    submission is data, and a circuit directory that happens to contain a `json.py` must not
    be able to shadow the stdlib. The verifiers each insert their own package root, so none
    of them needs that path anyway.
    """
    command = [
        sys.executable,
        "-I",
        str(verifier_script(verifiers_root, task)),
        str(target),
        "--salt",
        salt,
        *(() if seed is None else ("--seed", str(seed))),
    ]
    provenance = {"salt": salt, "seed": seed}
    started = time.monotonic()
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired:
        return Verification(
            task,
            "TIMEOUT",
            None,
            time.monotonic() - started,
            detail={"error": f"over {timeout}s"},
            **provenance,
        )
    seconds = time.monotonic() - started
    result = parse_verifier_output(task, completed.stdout, seconds, **provenance)
    if result.status == "UNPARSEABLE" and completed.stderr.strip():
        # A crash before the verifier's own error handler leaves nothing on stdout.
        return Verification(
            task,
            "CRASH",
            None,
            seconds,
            detail={"error": completed.stderr.strip()[-400:]},
            **provenance,
        )
    return result


def verify_submission(
    verifiers_root: Path | str,
    audit: SubmissionAudit,
    tasks: Sequence[str] = TASKS,
    *,
    salt: str = DEFAULT_SALT,
    seed: int | None = None,
    timeout: float | None = None,
) -> Iterator[Verification]:
    """Each task's verdict as it lands, so a caller can report an hour-long run as it goes."""
    for task in tasks:
        target = audit.target(task)
        if target is None:
            yield Verification(
                task, "MISSING", None, 0.0, detail={"error": "not in the submission"}
            )
            continue
        yield run_verifier(
            verifiers_root, task, target, salt=salt, seed=seed, timeout=timeout
        )


def recompute_cost(verification: Verification) -> float | None:
    """Re-derive a passing verifier's own score from the fields it reported.

    A **validation** of the transcription, not a predictor: the inputs the formulas need --
    the per-shot Toffoli mean, `rotation_eps`, `eps_plus` -- are each obtainable only by
    simulating, so this cannot price a candidate that has not been run. What it does is keep
    `cost` honest against the metric, permanently, over the three committed reference
    outputs as well as every live run.

    None when the verifier failed, or did not report the parts the formula needs.
    """
    if not verification.passed:
        return None
    try:
        return cost_from_report(verification.task, verification.detail)
    except (KeyError, TypeError, ValueError):
        return None


def model_agrees(verification: Verification) -> bool | None:
    """Whether `cost` reproduces this verdict's score; None if it could not be evaluated.

    Three states, deliberately: agreeing, disagreeing, and not evaluable are different
    facts, and collapsing the last two hides a model that never ran.
    """
    ours = recompute_cost(verification)
    if ours is None or verification.score is None:
        return None
    return math.isclose(ours, verification.score, rel_tol=MODEL_AGREEMENT_REL)


def passing_costs(verifications: Sequence[Verification]) -> dict[str, float]:
    """Just the costs that came back from a passing run, keyed by task."""
    return {v.task: v.score for v in verifications if v.passed and v.score is not None}


def measured_calibration(
    verifications: Sequence[Verification], provenance: str
) -> MetricCalibration:
    """Turn three passing reference verdicts into the calibration -- gate G0's output.

    The hashes come along so the result records *which* circuits it measured: calibrating
    against a compressed candidate instead of the reference would otherwise make that
    candidate score 1.0 and every later comparison silently relative to it.
    """
    costs = passing_costs(verifications)
    if missing := [t for t in TASKS if t not in costs]:
        failed = {v.task: v.error for v in verifications if not v.passed}
        raise ValueError(
            f"every task must pass to calibrate; missing {missing}"
            + (f", failures: {failed}" if failed else "")
        )
    hashes = {
        v.task: v.candidate_sha256
        for v in verifications
        if v.passed and v.candidate_sha256 is not None
    }
    return MetricCalibration(
        reference_costs=costs,
        provenance=provenance,
        measured=True,
        reference_sha256=hashes,
    )


def candidate_hashes(verifications: Sequence[Verification]) -> dict[str, str]:
    """What each passing run actually scored, keyed by task."""
    return {
        v.task: v.candidate_sha256
        for v in verifications
        if v.passed and v.candidate_sha256 is not None
    }


def write_rows(path: Path, verifications: Sequence[Verification]) -> None:
    """Rewrite the results file from scratch -- called after each verdict, not at the end.

    A full run is ~12 minutes of simulation; writing only at the end means a Ctrl-C or an
    OOM kill at minute 11 throws away everything already measured.
    """
    path.write_text(
        "\n".join(json.dumps(v.row()) for v in verifications) + "\n" if verifications else ""
    )
