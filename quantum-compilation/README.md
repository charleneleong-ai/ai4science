# quantum-compilation

Verifier-in-the-loop circuit compression for the [Google Quantum AI — Verifiable Quantum
Circuit Compilation](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation)
Kaggle competition ($150k, closes 2027-01-20).

Compress three fault-tolerant quantum subroutines — a 1/√x phasing circuit for fusion
plasma models, a Reed–Solomon syndrome decoder for Decoded Quantum Interferometry, and a
FeMoCo Hamiltonian block encoding. No quantum hardware, no training data, no GPU: every
candidate is checked exactly by [kickmix](https://github.com/quantumlib/Qualtran/tree/main/kickmix),
which simulates reversible-classical + diagonal-phase circuits of tens of millions of gates
in seconds. Paper: [arXiv:2610.23832](https://arxiv.org/abs/2610.23832).

This is the `ai4science` thesis with the two halves handed over pre-separated — an exact
verifier that cannot be gamed, and an unbounded generator space. The only question is how
well the search loop around it is built.

## Status — 2026-10-31

**Scaffold only. Nothing has been measured locally and nothing has been submitted.** The
competition data is not yet downloaded: there is no Kaggle API token on this machine and the
Kaggle MCP `authorize` tool fails in this client, so a token has to be installed by hand.

| | |
| --- | --- |
| Our best score | — (no submission) |
| Public frontier blend | **23.14** (*third-party*, 2026-10-30) |
| Leader | **26.76** (*third-party*, 2026-10-30) |
| Teams | 1,729 |

The grounded entry sweep is [`reports/2026-10-31-vqcc-landscape.md`](reports/2026-10-31-vqcc-landscape.md);
the plan and the one open decision are in
[`docs/specs/2026-10-31-quantum-compilation-design.md`](docs/specs/2026-10-31-quantum-compilation-design.md).

## The metric

Score is the geometric mean of three cost ratios against a reference submission that scores
exactly 1.0:

| Task | Cost |
| --- | --- |
| `invsqrt` | `Q × C`, with `C = Toffolis + T/4` |
| `dqi` | `R(ε₊) × Q² × T`, with `R(ε) = ln2 / −ln(½ + 2√ε)` |
| `dfthcsos` | `λ_eff × T × Q`, with `λ_eff = √(E_gap(2Λ − E_gap))` |

Transcribed in [`cost.py`](src/quantum_compilation/cost.py) from the
[Evaluation page](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/overview/evaluation)
and the [published metric](https://www.kaggle.com/code/metric/google-quantum-ai-metric). Two
consequences worth holding onto:

- **Equal weight.** A 1% cost cut in any task buys the same 0.33% of score.
- **`dqi` pays twice for width.** Its `Q²` term makes a qubit saved worth double.

The only validation so far is that [`score`](src/quantum_compilation/cost.py) reproduces the
published public blend — 23.131 against a reported 23.13587, from three ratios quoted to four
significant figures ([the test](tests/test_cost.py)). That confirms the *shape* of the
metric and nothing about its absolute calibration.

## Run it

```bash
uv sync --extra dev
uv run pytest            # 56 tests
uv run vqcc --help
```

Four verbs. `check` is the one that matters — it is the submission gate:

```bash
uv run vqcc check path/to/submission/      # layout, ANGLE widths, every cost its inputs allow
uv run vqcc resources path/to/circuits/    # qubits, Toffoli-likes, rotations by level
uv run vqcc blend --invsqrt 100941 --dqi 4.0334e11 --dfthcsos 3.3780e8
uv run vqcc calibration                    # the metric constants in force, and their provenance
uv run vqcc diagram circuit.kmb --flame    # kickmix's own SVG flame chart
```

`check` audits the layout the Evaluation page specifies — `invsqrt`, `dqi`, and a
`dfthcsos/` folder of `angle_k_0 … angle_k_<n-1>` with no gaps plus `submission.npz` — and
does it by filename alone, so auditing a 905 MB submission is instant and the gate runs
before anything expensive does. It then measures the circuits and prices what it can.

## Three deliberate refusals

Each exists because a plausible default would silently mis-rank every candidate.

**The metric's unpublished constants travel together, with provenance.** Each task's
reference cost and the price in T gates of a rotation Clifford+T cannot express exactly are
the same kind of unknown — a number the official metric owns and we have not read. Both live
in one loadable [`MetricCalibration`](src/quantum_compilation/cost.py) carrying a
`provenance` string and a `measured` flag, so they arrive by filling a JSON in rather than
by editing code, and every `vqcc` run prints a warning until `measured` is true. The shipped
one is back-derived from third-party ratios at ~0.1% slop. Write a copy to edit with
`vqcc calibration --write`.

**Pricing an unpriced rotation raises.** `toffoli_equivalents` refuses when the circuit holds
`Z_POW` past the T level and the calibration has no price for that level, naming the levels
it is missing.

**A gate checks what it was given, not what it assumed.** `dfthcsos_rejections` takes the
Hamiltonian error, `E_gap` and `Λ` as optional, because only the ANGLE widths are knowable
before the official verifier runs. Defaulting them to zero would invent rejections.

## What is not established

- No local measurement of anything, and no submission. The leaderboard is the only
  trustworthy score and we have no reading from it.
- The calibration is unmeasured — see `vqcc calibration` for exactly what is assumed.
- **The leaders hold ~55% more total cost reduction than the public frontier** (26.76/23.14
  on score is ×1.55 on the product of ratios) and none of it is visible. Which task it sits
  in is the most valuable unknown in the competition.
- Which task to commit to is **open** — see the design doc. The recommendation is
  `dfthcsos`, on skill fit and on evidence that the public line is stuck in one
  parameterization.

## Attribution

The public circuits are Apache-2.0 and reuse is permitted with attribution; there are live
threads on exactly this ([753660](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753660),
[753740](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753740)).
Anything derived here carries its lineage and its author's name. The landscape numbers above
come from [manwithacat's briefings](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754174),
a competitor's self-reported data, labelled as third-party wherever quoted.
