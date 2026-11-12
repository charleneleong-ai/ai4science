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

## Levers, priced

25.1x needs cost 82,740: **a 11.85% cut**. Each 1% of invsqrt cost is +0.336% on the blend.

| lever | size if it works | share of cost | status |
| --- | ---: | ---: | --- |
| merge the 8 repeated-weight rotations (add the two terms before kicking) | -57.5 C, plus the adder | up to 3.8% | untested; the adder costs Toffolis and maybe qubits |
| free accuracy margin so `eps` grows (each doubling saves 5.74 C) | 9.5 doublings to reach 50% slack at x=161 | up to 3.6% | upper bound: other inputs may bind first, and accuracy costs Toffolis |
| one fewer qubit (62 → 61) | 1/62 | 1.6% | worth up to ~25 extra Toffolis |
| shifter from 46 CSWAPs to the high-30s counting bound ([thread 757030](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757030)) | -8 C | 0.5% | SAT did not resolve w=21 publicly |

Together these are about 9.5% and **not additive** — the first two pull against each other through the
accuracy budget. They do not reach 11.85% on their own, so 25.1x most likely comes from a structural
change (fewer arithmetic Toffolis or several fewer qubits), not from trimming this design.

## What is not established

- Per-shot (not static) Toffoli counts per stage: needs a stage-tagged simulation.
- The error profile over all `x`, which decides how much the margin lever is really worth.
- What the repeated-weight terms are, and whether adding them first fits in 62 qubits.
- Whether mikelou1's 25.1x is real; it is a forum claim with no method.
