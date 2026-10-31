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

**The measurement loop is validated against the leaderboard; no compression has been
attempted yet.** All three official verifiers run locally over the competition's own
`sample_submission/`, every cost function here reproduces the verifier's own score to 1e-12
from the fields it reports, and submitting that unmodified reference scored **1.00000** —
exactly what the local verifier said. The local dev set can therefore be trusted to rank
candidates without spending submissions.

| | |
| --- | --- |
| Our best score | **1.00000** (reference submission `57824925`, local verifier 1.00000) |
| Public frontier blend | **23.14** (*third-party*, 2026-10-30) |
| Leader | **26.76** (*third-party*, 2026-10-30) |
| Teams | 1,729 |

The grounded entry sweep is [`reports/2026-10-31-vqcc-landscape.md`](reports/2026-10-31-vqcc-landscape.md);
the plan and the one open decision are in
[`docs/specs/2026-10-31-quantum-compilation-design.md`](docs/specs/2026-10-31-quantum-compilation-design.md);
verdicts are in [`docs/experiments/verdicts.yaml`](docs/experiments/verdicts.yaml).

## The metric

Score is the geometric mean of three cost ratios against a reference submission that scores
exactly 1.0. The reference costs below are **measured**, by running the shipped verifiers
over `sample_submission/` at the default salt:

| Task | Cost | Reference cost | Reference shape |
| --- | --- | --- | --- |
| `invsqrt` | `Q × C`, with `C = Toffolis + T/4` | 2,076,774.23 | Q=260, C=7,987.59 |
| `dqi` | `R(ε₊) × Q² × T`, with `R(ε) = ln2 / −ln(½ + 2√ε)` | 2.4439463583673 × 10¹³ | Q=1,893, T=5,725,108.4, R=1.1913 |
| `dfthcsos` | `λ_eff × T × Q`, with `λ_eff = √(E_gap(2Λ − E_gap))` | 3,355,599,588.32 | Q=730, T=93,237, λ_eff=49.3014 |

[`cost.py`](src/quantum_compilation/cost.py) is transcribed from the **shipped verifiers**,
not the [Evaluation page](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/overview/evaluation),
which states the formulas but not the constants. Three things it omits, each of which
changes the number:

1. **A rotation costs `0.56·log₂(1/ε) + 5.3` T gates** ([mixed fallback synthesis](https://arxiv.org/abs/2203.10064)),
   where `ε` is a *single circuit-wide* synthesis error — the largest that still fits the
   accuracy budget the simulated error leaves over — not a per-Clifford-level constant. The
   reference's 37 rotations land at `ε = 1.53 × 10⁻¹²`, so they cost 27.3 T each.
2. **`T` is the simulated mean per shot, not the static gate count.** Measurement-based
   uncomputation makes it data-dependent; the static count is only an upper bound. The dqi
   reference executes 5,725,108.4 of a possible 5,725,461.
3. **`ε₊` is a one-sided Clopper–Pearson upper bound** at `δ = 2⁻¹²⁸`, not the observed
   failure fraction. A decoder that fails **0 of 102,400** shots is still bounded at
   8.66 × 10⁻⁴ and pays **R = 1.1913** — a 19% floor nobody can optimise away.

Two consequences worth holding onto:

- **Equal weight.** A 1% cost cut in any task buys the same 0.33% of score.
- **`dqi` pays twice for width.** Its `Q²` term makes a qubit saved worth double.

## Run it

```bash
uv sync --extra dev
uv run pytest            # 111 tests
uv run vqcc --help
```

`verify` is the one that matters — it is the dev set:

```bash
export VQCC_VERIFIERS=path/to/verifiers      # the extracted competition `verifiers/` tree
uv run vqcc verify path/to/submission/ --seed 1   # the real score, free and unlimited
uv run vqcc check path/to/submission/        # layout + the 128-qubit cap, instant
uv run vqcc resources path/to/circuits/      # static gate counts — a triage bound
uv run vqcc blend --invsqrt 100941 --dqi 4.0334e11 --dfthcsos 3.3780e8
uv run vqcc calibration                      # the reference costs in force, and their provenance
uv run vqcc diagram circuit.kmb --flame      # kickmix's own SVG flame chart
```

**There is no train/dev split in this competition**, so the dev set is the *official verifier
run locally*. The leaderboard runs this same code and rations measurements to five a day;
running it locally is free, so every candidate goes through `verify` and only a winner costs
a submission. It is trustworthy because the reference submission scores 1.00000 both ways.

**Pass `--seed` to compare candidates.** Each verifier derives its simulator seed from the
*candidate's own bytes* — `int(circuit_sha256, 16)`, or `sha256(circuit_sha256 + salt)` when
a salt is given — so two different circuits draw two different seeds whatever the salt is
set to, and a re-run of one circuit is bit-identical. Left alone that puts an unmatched draw
inside every comparison and makes its noise simultaneously present and unmeasurable. All
three verifiers accept `--seed`, which overrides the hash; the salt only reproduces the
published default derivation. Reference run times on 10 CPU cores: invsqrt 5.7 min
(200,540,928 shots), dqi 2.4 min, dfthcsos 1 min.

`check` is the cheap pre-filter: it audits the layout the Evaluation page specifies —
`invsqrt`, `dqi`, and a `dfthcsos/` folder of `angle_k_0 … angle_k_<n-1>` with no gaps plus
`submission.npz` — from filenames and the npz index alone, so a structurally broken 905 MB
submission is rejected before anything is parsed. It then checks the one rejection knowable
without simulating: the 128-qubit ANGLE cap.

## Three deliberate refusals

Each exists because a plausible default would silently mis-rank every candidate.

**No cost is computed without the verifier — invsqrt included.** The scored Toffoli count is
a simulated mean and the rotation price depends on the accuracy margin the simulation leaves
over, so neither is knowable from a static reading of a circuit. `resources` reports static
counts and calls them a *bound*; [`verify`](src/quantum_compilation/verify.py) produces
costs. An earlier draft priced invsqrt statically and would have ranked candidates on a
number the leaderboard does not use.

**Every formula is checked against the verifier that owns it.**
[`recompute_cost`](src/quantum_compilation/verify.py) re-derives each verifier's own score
from the fields that verifier reported, and `vqcc verify` prints the comparison per task —
agreeing, disagreeing and *not evaluated* as three distinct states, so a model that never ran
cannot read as one that agreed. The three real reference outputs are committed as
[test data](tests/data/) so the check is a permanent regression test, not a one-off. It is a
*validation* of the transcription, not a predictor: the inputs it needs — the per-shot
Toffoli mean, `rotation_eps`, `eps_plus` — are obtainable only by simulating, so it cannot
price a candidate that has not been run.

**A gate checks what it was given, not what it assumed.** `dfthcsos_rejections` takes the two
accuracy errors, `E_gap` and `Λ` as optional, because only the ANGLE widths are knowable
before the official verifier runs. Defaulting them to zero would invent rejections.

## What is not established

- **No compression has been attempted.** Everything so far is measurement apparatus, and
  1.00000 is the reference's score, not an achievement.
- **Run-to-run noise is unmeasured.** A re-run at the same seed is bit-identical, so the
  only noise that exists is across seeds, and nothing has yet measured its size. Until it
  is, no local difference has a threshold to clear (gate G1).
- **One leaderboard reading, on an input designed to score 1.0.** That the local verifier
  and the leaderboard agree there does not prove they agree on a compressed candidate.
- **The leaders hold ~55% more total cost reduction than the public frontier** (26.76/23.14
  on score is ×1.55 on the product of ratios) and none of it is visible. Which task it sits
  in is the most valuable unknown in the competition.
- Which task to commit to is **open** — see the design doc. The recommendation is
  `dfthcsos`, on skill fit and on evidence that the public line is stuck in one
  parameterization. Two facts from the reference run bear on it: the ANGLE circuits are
  85% of that task's Toffoli cost (79,480 of 93,237) but only 119 of its 730 qubits, so
  width is governed by the state-prep workspace rather than by the circuits a submission
  ships; and the reference sits close to two cliffs, `error_2` at 0.01864 against a 0.020
  cap and `E_gap` at 3.9655 against a 3.914898 floor.

## Attribution

The public circuits are Apache-2.0 and reuse is permitted with attribution; there are live
threads on exactly this ([753660](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753660),
[753740](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753740)).
Anything derived here carries its lineage and its author's name. The landscape numbers above
come from [manwithacat's briefings](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754174),
a competitor's self-reported data, labelled as third-party wherever quoted. Those briefings
held up well: their implied reference costs came in 0.02%, 0.006% and 0.04% below what the
verifiers measured.
