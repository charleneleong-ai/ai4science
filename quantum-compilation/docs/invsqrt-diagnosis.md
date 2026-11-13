# invsqrt — cost anatomy and levers (2026-11-12)

The [frontier sweep](../reports/2026-11-12-vqcc-frontier-delta.md) picked invsqrt as the first task to
attack: it is the only one with reported headroom above public (mikelou1 claims 25.1x in
[thread 756998](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756998),
method not given). This breaks down the circuit we submitted, `6c6a734d` (manwithacat HW-25, ratio
22.1255), to see where 25.1x could come from. Numbers are reproduced by
[`scripts/invsqrt_anatomy.py`](../scripts/invsqrt_anatomy.py) from the circuit and our verifier row.

## How the cost splits

invsqrt cost is `Q x C` with `C = mean Toffolis per shot + T/4`
([verifier](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/data)).
For `6c6a734d`: **Q 62 x C 1,513.93 = 93,863.4**.

| part | C | share |
| --- | ---: | ---: |
| Toffoli-like, mean per shot (1,380 static: 1,242 CCX, 92 CCZ, 46 CSWAP) | 1,218.83 | 80.5% |
| 41 arbitrary rotations, 28.76 T each at the verifier's eps 2.43e-13 | 294.84 | 19.5% |

**The rotations are a bit-serial phase kickback.** All 41 act on one qubit (`q28`) as `Z^a` with
`a` = 2^-3 .. 2^-35 (a phase of `a/2` turns): each bit of the computed 1/sqrt(x) is brought to `q28` in
turn and rotated by its weight.
Exponents 2^-10 .. 2^-17 appear **twice**, so 8 of the 41 rotations kick a second term of the same
weight.

**The rotation price is set by the accuracy margin, and the margin is spent.** The verifier picks the
largest `eps` with `n_rot * 2 eps` inside the tolerance the simulated error leaves over; the worst
input (`x = 161`) already uses 99.932% of its tolerance, so `eps` is 2.4e-13 and each rotation costs
28.8 T (7.19 C, 0.475% of the task cost).

Static stage map (instruction index; Toffoli-like counts are static, not per-shot means):

| instructions | what is there | Toffoli-like |
| --- | --- | ---: |
| 0 – 1.5k | normalising shifter: all 46 CSWAPs (x0..x20, control q48) + 54 CCX | ~100 |
| 1.5k – 7.5k | steady arithmetic, no rotations | ~132 |
| 7.5k – 17k | arithmetic with heavy HMR uncompute; first 25 rotations | ~470 |
| 17k – 25.5k | conditioned multiplies (CCZ, HMR); last 19 rotations | ~570 |
| 25.5k – 33k | Clifford only | ~1 |
| 33k – 37k | tail uncompute, mostly conditioned | ~145 |

## Qubit liveness

[`scripts/invsqrt_liveness.py`](../scripts/invsqrt_liveness.py) steps the circuit one gate (or
condition block) at a time over 16,384 inputs and records which qubits are nonzero. Each qubit is a
definite bit per shot, so its occupied stretches are intervals and the peak overlap is the qubit
count reachable by relabelling alone. A qubit is *parked* where it is occupied with no
unconditioned use; HMR uncompute and its phase fixups touch qubits without using them.

**The peak is the declared 62, held in 12 runs of 50+ steps; relabelling frees nothing.** The long
run, instructions 1,663–7,493 (5,831 steps), writes a Boolean function of the top 8 normalised bits
(`q14`–`q21`) into 39 outputs (`q22`–`q61`) with 5,627 CX and 129 CCX; the outputs are first used
together at 8,034–8,074, so they cannot be loaded lazily. The short runs sit in kick pass 1
(8.6k–10.8k) and pass 2 (16.5k–16.7k).

**What fills it is information, not slack.** Through the long run and the five pass-1 runs at
9.8k–10.8k the same 15 qubits are parked: the low input bits `q0`–`q13` (last used ~700, next
~8,570) and `q23`. `q23` starts as `NOT(x20 XOR x21)` and steers the shifter, and by 1,663 it is
not recomputable: 23 states of the other 61 qubits occur with both of its values. So the long run
holds 14 input bits, 1 shift bit, 8 address bits and 39 lookup outputs, all needed at once.

## Levers, priced

25.1x needs cost 82,740: **a 11.85% cut**. Each 1% of invsqrt cost is +0.336% on the blend.

| lever | measured | verdict |
| --- | --- | --- |
| free accuracy margin so `eps` grows | **0.44% at most** (fixing the 1,000 tightest inputs frees 6.6 C) | closed |
| merge the 8 repeated-weight rotations | saves 57.5 C, but the 8 bits would cross 62-qubit stretches between the passes | closed |
| one fewer qubit (62 → 61) | 1.6%, worth up to ~25 extra Toffolis; no free qubit exists (above), so all 12 peak runs must drop one | needs a narrower lookup |
| shifter from 46 CSWAPs to the high-30s counting bound ([thread 757030](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757030)) | 0.5% | open, small |

**The margin is spent everywhere, not at one input.** [`scripts/invsqrt_error_profile.py`](../scripts/invsqrt_error_profile.py)
reruns the verifier's own shot loop over all 3,133,452 inputs and keeps each signed error. It reproduces
the verifier's `eps` (2.434e-13) exactly, and the error is deterministic per input (no input varies with
measurement outcomes). The smallest margin, 3.177e-12 turns, is shared by at least 30 inputs; 168,010
inputs use more than half their tolerance and 184 more than 90%. Freeing margin means more accurate
arithmetic everywhere, which costs Toffolis.

**The repeated weights are two kick passes, not one doubled bit.** Pass 1 (instructions 8.6k–9.0k) kicks
exponents 2^-17 .. 2^-3; pass 2 (16.4k–16.7k) kicks 2^-35 .. 2^-10. Each kick loads a bit into `q28` with a
Toffoli (`CX q28 q4; CCX q21 q4 q28; CX q28 q4`, cycling a 4-qubit window q4–q7), rotates, and erases it
with `HMR`. So the phase is the sum of two terms computed ~7,800 instructions apart.

**Conclusion: trimming this design is closed.** The open levers (rows 3-4) total ~2%, against 11.85%
for 25.1x. That gap needs a structural change — fewer arithmetic Toffolis in the two ~500-Toffoli sections, or
several fewer qubits, which only a narrower lookup can give — not tuning of the kick or the margin.

## What is not established

- Per-shot (not static) Toffoli counts per stage: needs a stage-tagged simulation.
- What the 39 lookup outputs encode (likely piecewise-polynomial coefficients over 256 segments), so a narrower lookup can be priced.
- What the two phase terms are mathematically (likely a coarse estimate plus a correction).
- Whether mikelou1's 25.1x is real; it is a forum claim with no method.
