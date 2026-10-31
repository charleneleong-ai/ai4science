# quantum-compilation — design

**Date:** 2026-10-31 · **Status:** scaffold landed, task choice open

## Intent

Enter the [Google Quantum AI Verifiable Quantum Circuit Compilation](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation)
competition (closes 2027-01-20) with a verifier-in-the-loop search for smaller circuits.

**Why this competition, in this repo.** It is the `ai4science` thesis with the generator and
verifier handed over pre-separated: an exact classical verifier (`kickmix`) that accepts or
rejects a circuit in seconds, and an unbounded generator space. There is no model to train
and no label noise to argue about — the only question is how well a search loop is built
around a ground truth that cannot be gamed. That is the pattern
[`touchstone/`](../../../touchstone) argues for and [`casmi/`](../../../casmi) learned the
hard way when a local estimate of 0.664 scored 0.115.

## Goals

1. A local compile-and-verify loop that reproduces a known public score before trusting any
   number it produces.
2. Depth on **one** task, aimed at that task's $10,000 single-task prize, with the other two
   riding the public Apache-2.0 frontier.
3. Every result traceable: a funnel with denominators, a verdict per experiment, and the
   public/private score gap stated rather than smoothed.

## Non-goals

- **Not** chasing the overall leaderboard top-6 from a standing start. 1,729 teams, the
  leader at 26.76 against a public frontier of 23.14, and ~55% of unexplained cost reduction
  held privately. A single-task prize is a reachable target; sixth place overall is not the
  thing to optimise for.
- **Not** the mid-competition milestone prize (2026-11-20). It needs the highest *public*
  score 20 days from an empty repo.
- **Not** hand circuit-surgery on `invsqrt`. Six-plus strong competitors are doing
  single-Toffoli surgery on one shared lineage with frontier coding models; marginal return
  for a late entrant is the worst of the three.

## Anti-goals — what makes this a failure even if it scores

- **A local score we believe over the leaderboard.** The exact failure in `casmi`. The
  verifier is seeded from the submission hash plus a secret salt and **final scoring
  re-draws**, so any number measured at the default salt is one sample, not the score.
- **An unattributed derivative.** The public circuits are Apache-2.0 and reuse is permitted,
  but there are live licensing threads
  ([753660](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753660),
  [753740](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753740)).
  Anything derived carries its lineage and its author's name.
- **A reported gain below submission noise.** 5 submissions/day, ~45 min to score. One
  change per submission, and a measured re-run to establish noise before any delta is
  called real.

## Constraints

- **Measurement is rationed**: 5 submissions/day. The external score is the only
  trustworthy one.
- **CPU-only is sufficient** — `kickmix==1.0.0` ships a cp313 arm64 macOS wheel and
  simulates tens of millions of gates in seconds. No GPU, no quantum hardware.
- **No Kaggle API token on this machine.** The MCP `authorize` tool fails in this client
  (mcp-remote schema mismatch), so the competition data needs a token installed by hand
  before anything can be measured locally. This blocks step 1 below and nothing else.
- Hard feasibility on `dfthcsos`: `Λ > 49.0`, `E_gap > 3.914898` Ha, Hamiltonian within
  20 mHa, every ANGLE circuit ≤ 128 qubits. A violation is a rejection, not a bad score.

## The decided how

**Scoring is transcribed, not inferred.** [`cost.py`](../../src/quantum_compilation/cost.py)
holds the three cost functions and the geometric mean from the official Evaluation page.
Its only validation so far is that it reproduces the published public blend (23.131 against
a reported 23.13587) from three third-party ratios — which is enough to confirm the shape of
the metric and nothing more.

**Resources are measured, never estimated.** [`resources.py`](../../src/quantum_compilation/resources.py)
reads qubits, Toffoli-likes and rotation counts straight off the circuit via kickmix, taking
kickmix's own accessors so kickmix owns the Clifford-hierarchy semantics. It holds no metric
coefficients at all — not even the `/4` — so there is exactly one place a cost constant can
hide, and it is `cost.py`.

**The metric's unknowns travel together.** Each task's reference cost and the T price of a
rotation Clifford+T cannot express exactly are the same kind of unknown: a number the
official metric owns and we have not read. Both sit in one
[`MetricCalibration`](../../src/quantum_compilation/cost.py) with a `provenance` string and
a `measured` flag, loadable from JSON — so the measured numbers arrive by **filling the
object in**, not by editing code or deleting a code path, and every `vqcc` run warns until
`measured` is true. A `Callable` pricing hook was the first attempt and was wrong for
exactly one reason: it cannot be loaded from the file the plan says will replace it.

**The gate runs before the expensive work.** [`submission.py`](../../src/quantum_compilation/submission.py)
audits a submission directory by filenames and the npz index alone — no circuit parsing — so
a structurally broken 905 MB submission is rejected instantly. `vqcc check` then measures
and prices, and reports per task which inputs only the official verifier can supply rather
than assuming any of them.

**Repo placement — my call.** `ai4science/quantum-compilation/`, as a self-contained
subproject matching [`casmi/`](../../../casmi)'s layout (own `pyproject.toml`, `uv.lock`,
`src/`, `tests/`, `docs/`, `reports/`, `research_notes/`). Rationale: it is a Kaggle
competition and `casmi` is the established pattern for one here, and the verifier-first
thesis is the repo's through-line. It is not in `quantum-ml-research` — that repo is a QRC
benchmark with a different question and a different claim to defend. Easy to relocate if
you'd rather it stood alone.

**Python 3.13+, uv, typer, rich, pytest** — repo convention. `ruff` complexity gate comes
from the repo-root `.pre-commit-config.yaml`.

## Left to the builder

**The task choice is open, and it is the decision that matters.** The read from the
[landscape report](../../reports/2026-10-31-vqcc-landscape.md):

| Task | Case for | Case against |
| --- | --- | --- |
| `dqi` | `Q²` makes a qubit worth double. 888 → 800 is ~23% off the task, ~7% on score. An independent 1,797-qubit design already passes the verifier, so the space is not closed. | Reed–Solomon decoder over GF(2¹²) is real domain depth. `R(ε₊)` couples cost to a re-drawn failure rate, so gains must be shown across many draws. |
| `dfthcsos` | Best fit for an optimization/ML background — a joint search over tensor factorization, Givens tree and QROM tables, not circuit surgery. Public gains have collapsed to 0.001–0.2% weight tweaks on one shared word table, which is what a local optimum of one parameterization looks like. Two new structures appeared on 2026-10-30, neither at the frontier. | Hard feasibility constraints reject rather than penalise. Cost spans the whole block encoding, so the submitted ANGLE circuits are only part of what is scored. |
| `invsqrt` | Smallest, most legible circuit; fastest iteration. | The crowded one. See non-goals. |

My recommendation is **`dfthcsos`**, on skill fit and on the evidence that the public line
is stuck in one parameterization — but this is a judgement call on where you want to spend
three months, not something the data settles.

Then, in order:

1. **Install a Kaggle token, download the data, measure `sample_submission/`.** It must come
   out at exactly 1.0 — if it doesn't, the cost model is wrong and nothing downstream is
   worth running. That is gate G0.
2. **Fill the calibration**: the measured reference costs, and the rotation→T price read off
   [the published metric](https://www.kaggle.com/code/metric/google-quantum-ai-metric).
   `vqcc calibration --write` emits the file to edit; set `measured: true` only once both
   are real, because that flag is what silences the warning on every run.
3. **Reproduce the public frontier blend locally** from the Apache-2.0 circuits dataset, with
   attribution, and submit it once. That calibrates the local loop against an external score
   and establishes the baseline every later claim is measured from.
4. **Measure submission noise** — one identical re-run — before calling any delta real.
5. Only then start searching, on one task.

Standing obligations from the repo rules: a weekly delta sweep (next due **2026-11-07**), a
verdict row in [`docs/experiments/verdicts.yaml`](../experiments/verdicts.yaml) when each
experiment closes, and the W&B project / artifact / PR description refreshed together at
each result that moves the headline.
