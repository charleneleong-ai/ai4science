# VQCC landscape — grounded 2026-10-31

Entry sweep for the [Google Quantum AI — Verifiable Quantum Circuit Compilation](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation)
competition. Raw captures in [`research_notes/2026-10-31-vqcc-landscape/`](../research_notes/2026-10-31-vqcc-landscape).

**Provenance labels used throughout.** *Official* = the competition's own pages or the
published metric. *Third-party* = another competitor's forum post, re-verified by them and
not by us. Nothing here has been measured by us yet — the competition data is not yet
downloaded (no Kaggle API token on this machine).

## The task

Compress three fault-tolerant quantum subroutines. No quantum hardware, no training data,
no GPU: circuits are *verified classically* by [kickmix](https://github.com/quantumlib/Qualtran/tree/main/kickmix),
which simulates reversible-classical + diagonal-phase + measurement-uncomputation circuits
of tens of millions of gates in seconds. The companion paper is
[arXiv:2610.23832](https://arxiv.org/abs/2610.23832).

Score is the geometric mean of three cost ratios against a reference submission that scores
1.0 (*official*):

| Task | Domain | Cost |
| --- | --- | --- |
| `invsqrt` | fusion — 1/√x phasing for Coulomb interactions | `Q × C`, `C = Toffolis + T/4` |
| `dqi` | combinatorial optimization — Reed–Solomon syndrome decoder | `R(ε₊) × Q² × T` |
| `dfthcsos` | chemistry — FeMoCo Hamiltonian block encoding | `λ_eff × T × Q` |

with `R(ε) = ln2 / −ln(½ + 2√ε)` and `λ_eff = √(E_gap(2Λ − E_gap))`.

Two structural facts that set strategy:

- **Equal weight.** A geometric mean means a 1% cost cut in *any* task buys the same 0.33%
  of score. There is no task worth more than another per unit of relative progress.
- **`dqi` pays twice for width.** Its `Q²` term makes a qubit saved worth double what it is
  worth in the other two tasks, which are linear in `Q`.

## Where the field is (2026-10-30, *third-party*)

From [manwithacat's briefing, topic 754174](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754174)
— a Claude-run research-ops agent publishing daily leaderboard and per-task frontier
snapshots. Easily the best public intel in the competition; it is also a competitor, so
treat its numbers as unverified by us.

| | Value |
| --- | --- |
| Top score | **26.76409** (Breaking1800) |
| Gold line (≤ rank 13) | 24.74 |
| Silver line (≤ rank 84) | 23.137 |
| Median | 22.486 |
| Teams | 1,696 (1,729 as of today) |

Public per-task frontier, and the blend it implies:

| Task | Qubits | Cost | Ratio |
| --- | --- | --- | --- |
| `invsqrt` | 62 | 100,941 | 20.57 |
| `dqi` | 888 | 4.0334 × 10¹¹ | 60.59 |
| `dfthcsos` | 430 | 3.3780 × 10⁸ | 9.93 |
| **blend** | | | **23.136** |

Our [`score`](../src/quantum_compilation/cost.py) reproduces 23.131 from those three ratios —
0.02% off the briefing's 23.13587, consistent with its ratios being quoted to four
significant figures. That is the only check we have so far that the metric is transcribed
correctly, and it is [a test](../tests/test_cost.py).

**The leaders hold ~55% more total cost reduction than the public frontier** (26.76/23.14 =
1.157 on score, so 1.157³ = 1.55 on the product of ratios), and none of it is visible. Where
that 55% sits — one task or spread across three — is the most valuable unknown in the
competition.

## Read on each task

**`invsqrt` — crowded, hand-optimized, integral.** At 62 qubits a single qubit is ~1.6% of
the task. Public progress is now single-Toffoli surgery on one shared lineage (mrmeows/team
selenium → busyaprime → yoyobar's 62-qubit → Özer's −1 Toffoli), much of it done with
frontier coding models. Lowest marginal return for a late entrant; highest density of people
already doing exactly this.

**`dqi` — the one structural lever in plain sight.** 888 qubits with a `Q²` term: 888 → 800
would be `(888/800)² = 1.23`, i.e. 23% off the task and ~7% on score. Leo Tang's
independent 1,797-qubit decoder (ratio ~10 vs 60) proves non-lineage designs pass the
verifier, so the design space is not exhausted. The subtlety is `R(ε₊)`: failures are drawn
from the file's own hash, so any byte change re-rolls them, and **final scoring re-draws**.
Judge a decoder on expected failures across many draws, never the draw you got.

**`dfthcsos` — the best fit for an ML/optimization background.** This is not circuit
surgery: it is a joint classical optimization over a tensor factorization, a Givens-rotation
tree and QROM lookup tables, scored by `λ_eff × T × Q` subject to hard feasibility
(`Λ > 49.0`, `E_gap > 3.914898` Ha, Hamiltonian within 20 mHa, every ANGLE circuit ≤ 128
qubits). Public gains have collapsed to 0.001–0.2% weight tweaks on one shared word table,
which is what a local optimum of a single parameterization looks like. Two new structures
appeared on 2026-10-30 (Den Pugovkin at 423 qubits, Leo Tang at 742) and neither is at the
frontier yet — the structure search is open.

## Prizes and dates (*official*)

$150k total. Leaderboard 1st–6th ($50k/25k/10k/5k/5k/5k), plus **three $10,000 single-task
prizes** — one per task, for the best selected submission on that task alone.

| Date | What |
| --- | --- |
| 2026-11-20 | Mid-competition milestone ($10k) — highest public score, notebook must be public |
| 2027-01-04 | Longest-leader prize period ends; qualifying notebooks must be public |
| 2027-01-13 | Entry and team-merger deadline |
| 2027-01-20 | Final submission |

5 submissions/day, max team size 5, scoring takes ~45 minutes.

**The single-task prizes are the realistic target.** Beating 26.76 overall from a standing
start against 1,729 teams is a different proposition from being best on one task, and the
geometric mean means a single-task specialist still scores respectably overall if the other
two ride the public frontier.

## What is not established

- **Nothing has been measured locally.** No competition data on this machine; the Kaggle MCP
  `authorize` tool fails in this client, so a token has to be installed by hand.
- **The metric is uncalibrated.** Both unpublished constants — each task's reference cost and
  the T price of a rotation Clifford+T cannot express exactly — carry their own provenance in
  [`MetricCalibration`](../src/quantum_compilation/cost.py); run `vqcc calibration` for the
  current state rather than trusting a copy of it written here.
- **Leader composition is unknown**, and that is the 55% that matters.
- **Licensing is live.** Public circuits are Apache-2.0 and reuse is permitted with
  attribution, but there are open threads on exactly this
  ([753660](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753660),
  [753740](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753740)).
  Read them before shipping anything derived.

## Next sweep

Weekly delta per the standing rule, and before any decision to commit to one task. Next due
**2026-11-07**; the questions to re-ask are whether the leader gap has narrowed, whether a
second independent `dqi` design has reached the frontier, and whether any `dfthcsos`
structure has displaced the shared word table.
