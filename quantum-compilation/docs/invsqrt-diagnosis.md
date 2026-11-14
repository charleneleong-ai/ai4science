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

Static stage map (instruction index; Toffoli-like counts are static, not per-shot means; per-shot
counts are under "Measured constants"):

| instructions | what is there | Toffoli-like |
| --- | --- | ---: |
| 0 – 1.5k | normalising shifter: all 46 CSWAPs (x0..x20, control q48) + 54 CCX | ~100 |
| 1.5k – 7.5k | steady arithmetic, no rotations | ~132 |
| 7.5k – 17k | arithmetic with heavy HMR uncompute; first 25 rotations | ~470 |
| 17k – 25.5k | conditioned multiplies (CCZ, HMR); last 19 rotations | ~570 |
| 25.5k – 33k | Clifford only | ~1 |
| 33k – 37k | tail uncompute, mostly conditioned | ~145 |

## What the circuit computes

[`scripts/invsqrt_decompose.py`](../scripts/invsqrt_decompose.py) reads the lookup's outputs and
the bit on `q28` before each fine rotation over 16,384 inputs. Write `x = 2^e · M` with `M` in
[1, 2):

1. **Shifter** (0 – 1.5k): normalises `x` by a power of two; `q23` keeps the parity of `e`, which
   is why it cannot be recomputed later.
2. **Lookup** (1.7k – 7.5k): a 256-entry table on the 7 mantissa bits after the leading one plus
   that parity. It has 38 output bits, no collisions at 7 bits, 128 at 6; the outputs are not
   a linear (XOR) encoding of samples of `1/sqrt(M)`, so the coefficients are not decoded bit by bit.
3. **Pass 1** (8.6k – 9.0k): kicks a seed, about 8 bits of `1/sqrt(x)` (x = 100: 0.099609 of
   0.1), scaled by `2^-floor(e/2)`.
4. **Arithmetic and pass 2** (to 16.7k): kicks the correction. Per segment it is a **quadratic in
   the 14 low mantissa bits**: a degree-2 fit leaves 0.69x tolerance at worst, degree 1 leaves
   70x, and degree 3 gains nothing (0.70x), so the remainder is fixed-point rounding.

So this is a seed-plus-quadratic design: 256 segments, two multiplies by a 14-bit `u`, and the
seed kicked and erased before the multiplies need the room.

## Qubit liveness

[`scripts/invsqrt_liveness.py`](../scripts/invsqrt_liveness.py) steps the circuit one gate (or
condition block) at a time over 16,384 inputs and records which qubits are nonzero. Each qubit is a
definite bit per shot, so its occupied stretches are intervals and the peak overlap is the qubit
count reachable by relabelling alone. A qubit is *parked* where it is occupied with no
unconditioned use; HMR uncompute and its phase fixups touch qubits without using them.

**The peak is the declared 62, held in 12 runs of 50+ steps; relabelling frees nothing.** The long
run, instructions 1,663–7,493 (5,831 steps), is the bulk of the lookup (step 2 above, which starts
at 688): the 8 address qubits `q14`–`q21` and 38 outputs in `q22`–`q61`, written with CX and 167
CCX and a scratch `q41` erased at 7,493; the outputs are first used
together at 8,034–8,074, so they cannot be loaded lazily. The short runs sit in kick pass 1
(8.6k–10.8k) and pass 2 (16.5k–16.7k).

**What fills it is information, not slack.** Through the long run and the five pass-1 runs at
9.8k–10.8k the same 15 qubits are parked: the low input bits `q0`–`q13` (last used ~700, next
~8,570) and `q23`. `q23` starts as `NOT(x20 XOR x21)` and steers the shifter, and by 1,663 it is
not recomputable: 23 states of the other 61 qubits occur with both of its values. So the long run
holds 14 input bits, the parity bit, 8 address bits, 38 lookup outputs and 1 scratch bit: 62.

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

## Redesign cost model

[`invsqrt_model.py`](../src/quantum_compilation/invsqrt_model.py) models the seed-plus-series
family bit-exactly: seed `y0` from a table on `A` mantissa bits plus parity, `eps = 1 - y0^2 v`,
then `y = y0 (1 + eps/2 + 3 eps^2/8 + ...)` with every product floored at `2^-F` and the phase
floored at `2^-K` turns. Each variant is scored on all 3,133,452 inputs against the verifier's own
target and tolerance, so feasibility is exact. Toffolis are counted as partial products of
truncated multipliers plus adder widths, at two per-unit constants measured on `6c6a734d` (below).
The model's nearest feasible twin of `6c6a734d` (A 7, 8-bit seed, 3 terms, F 27, K 36) calibrates
the rest: the measured 1,218.83 Toffolis and the 62-qubit arithmetic peak are matched by offsets
applied to every variant.

### Measured constants

[`scripts/invsqrt_stage_counts.py`](../scripts/invsqrt_stage_counts.py) simulates each stage of
`6c6a734d` on 131,072 uniform inputs and counts the Toffolis that fire per shot (conditioned gates
counted only when their condition holds). The stages sum to the verifier's 1,218.83.

| stage | Toffolis per shot |
| --- | ---: |
| shifter | 65.00 |
| lookup | 167.00 |
| lookup cleanup | 7.00 |
| pass 1 kick | 51.00 |
| arithmetic | 371.99 |
| pass 2 kick | 78.00 |
| uncompute | 408.84 |
| Clifford kicks | 0.50 |
| tail | 69.50 |

- **Lookup: 0.652 Toffolis per entry** (167 for 256 entries; see "How the lookup works" below).
- **Uncompute costs more than the forward arithmetic** (409 vs 372). It scales with the
  arithmetic, so it belongs in the per-product constant, not in the fixed offset.
- **1.44 Toffolis per partial product**: pass 1, arithmetic, pass 2 and uncompute (910), less the
  44 kicks and 105 adds, over the twin's 530 partial products
  ([`CostModel.fitted`](../src/quantum_compilation/invsqrt_model.py), pinned by a test). If the adds are uncomputed too, the
  figure is 1.24. Forward-only it is 0.66.

### Search

[`scripts/invsqrt_design_search.py`](../scripts/invsqrt_design_search.py) searches 300
configurations (A 5–9, seed 7–11 bits, 2–3 terms, three table contents, Taylor or fitted
coefficients) over precision and floor: 2,527 variants pass on every input.

The lookup is priced two ways: at `6c6a734d`'s synthesised 0.652 per entry, or as kickmix's
`init_lookup` (2^a − 2 for 2^a entries). The calibration offset leaves the measured lookup out, so
`6c6a734d` keeps its own 167 under both.

| Toffolis per partial product | lookup | best variant | vs reference |
| ---: | --- | --- | ---: |
| 1.00 | synthesised | A 8, 10-bit seed, 2 terms, F 26, K 36; Q 58 | −3.7% |
| 1.24 | synthesised | A 8, 10-bit seed, 2 terms | −5.1% |
| **1.44** | **synthesised** | **A 8, 10-bit seed, 2 terms; Q 58, C 1,530** | **−6.3%** |
| 2.00 | synthesised | A 8, 10-bit seed, 2 terms | −9.8% |
| 1.00 | `init_lookup` | A 6, 8-bit seed, 3 terms, F 27, K 36; Q 62 | −0.2% |
| 1.44 | `init_lookup` | A 6, 8-bit seed, 3 terms; Q 62, C 1,537 | +0.6% |
| 2.00 | `init_lookup` | A 8, 10-bit seed, 2 terms; Q 58 | +1.0% |

**At the measured constants the redesign models at −6.3%, about half the −11.85% that 25.1x
needs.** That prices its 512-entry table at 0.652 per entry (334 Toffolis), which only a
synthesis like `6c6a734d`'s could reach (below). Built with `init_lookup` no variant beats
`6c6a734d` at the measured constants.
The winner is a 512-entry table with a 10-bit seed, which shrinks `eps` to ~2^-10 so two Taylor
terms suffice. That drops the `eps^3` multiply and narrows the `eps` and sum registers (modelled
Q 58, from 62), paid for by a lookup twice the size. At −6.3% the invsqrt ratio would rise from
22.13 to about 23.6, lifting the overall score by about 2% (the score is a geometric mean over three
tasks). Closing the rest needs a cheaper multiplier or uncompute than `6c6a734d`'s, not a different
series.

## Where the uncompute goes

The multipliers already use measurement-based uncomputation (temporary ANDs undone by `HMR` plus a
conditioned `CZ`, [arXiv:1709.06648](https://arxiv.org/abs/1709.06648)): the arithmetic and
uncompute hold 588 `HMR` and ~2,900 conditioned `CZ`.

[`scripts/invsqrt_uncompute_rewind.py`](../scripts/invsqrt_uncompute_rewind.py) fingerprints every
qubit's value across 1,024 shots after every step, and tracks the latest forward value the
uncompute still holds. Per-qubit attribution fails here (qubits 14–61 are each rewritten in four or
more stages), as does whole-state matching (no uncompute state but the first equals a forward one).

| uncompute Toffolis spent | rewound to forward op | forward Toffolis undone |
| ---: | ---: | ---: |
| 0 | 16,417 (pass 2 undone free) | 81 |
| 161 | 15,024 | 229 |
| 323 | 13,268 | 344 |
| 404 | 12,689 | 407 |
| 409 | 8,577 (all of the arithmetic) | 502 |

- **It is a full reversal.** The uncompute walks back steadily to op 12,689, then releases the
  arithmetic's first 3,700 ops by measurement for ~5 Toffolis.
- **The reversal is a mirror, within 3%.** Undoing a measurement-based circuit recomputes each AND
  the forward measured away, so a mirror costs the forward `HMR` count over the span it rewinds,
  not its Toffolis: 332.0 unconditioned uncompute Toffolis against 323 `HMR` in ops 8,577–16,417.
- **The only non-mirror cost is 76.8 conditioned Toffolis** (5% of C): phase fixups for bits
  measured out in the X basis whose value is not a single AND.
  [`scripts/invsqrt_fixups.py`](../scripts/invsqrt_fixups.py) costs each of the 65 blocks and
  traces its condition bit back to the measurement outcomes it depends on:

  | condition bit | blocks | fired per shot | fires |
  | --- | ---: | ---: | ---: |
  | one outcome | 45 | 60.1 | ~half |
  | parity of 4–8 outcomes | 8 | 4.5 | ~half |
  | OR of 2–29 outcomes | 12 | 12.4 | 76–100% |

  The single-outcome blocks are ripple-carry chains (`CX, CX, CCX` majority steps with `Z` kicks):
  a measured-out sum bit depends on the carry from every lower bit, so fixing its phase rebuilds
  that carry. A fair coin decides whether each fires, which is the halving measurement-based
  uncomputation buys. The OR-gated blocks, mostly one CCX each, fire as 1 − 2^-n for n outcomes,
  so their condition saves almost nothing.

So there is no slack in the uncompute as built. Cutting it means computing fewer measured-away
ANDs in the forward, which the redesign's narrower products already model.

## How the lookup works

[`scripts/invsqrt_lookup_decode.py`](../scripts/invsqrt_lookup_decode.py) runs one input per entry
through instructions 688–7,494.

- **Address:** at 688, `q14`–`q20` hold mantissa bits 14–20 and `q21` the exponent parity, plainly
  (degree 1, all 256 values distinct); the 65 Toffolis before it are the shifter's.
- **Not a QROM:** 167 CCX, of which 166 write directly into the 38 output qubits and 1 into the
  scratch `q41`; no ancilla tree and no measurement-based uncompute (one HMR, erasing `q41`). CX
  networks rewrite the address register between them (it reaches degree 8 mid-lookup) and restore
  it by the end. No CCX is conditioned, so every shot costs 167.
- **The table is not simple:** 20 outputs have algebraic degree 7 and 18 degree 8 in the address,
  and all 247 nonlinear monomials appear across them. So 4.4 ANDs per output bit is bought by
  sharing ANDs across outputs in place, a multiplicative-complexity-minimised XOR-AND circuit
  that a logic-synthesis tool would emit, not by structure in the seeds.

A generic `init_lookup` would cost 254 here, so the synthesis saves 87 Toffolis (7% of 1,218.83). It
cannot be reused for another table; a redesign must be synthesised afresh.

## What is not established

- What a synthesis like `6c6a734d`'s would cost on the redesign's 512 × 29-bit table: the
  0.652-per-entry rate is one measurement on one table, not a scaling law.
- Whether the 12.4 OR-gated fixup Toffolis (under 1% of C) can be avoided. Even if all of them
  could, the redesign would still be ~2.5% short of 11.85%.
- Whether a built 512-entry design matches its modelled gain: the per-product constant is
  measured on `6c6a734d`'s multipliers, whose layout a new design may not share.
- How the 38 lookup bits split between the seed and the correction, and why `6c6a734d` meets the
  tolerance with an 8-bit seed and a quadratic where the model needs 3 Taylor terms (likely seeds
  or coefficients tuned per segment); so the model's reference is a twin, not the circuit.
- The model's qubit count beyond the lookup peak is calibrated, not derived from a layout.
- Whether mikelou1's 25.1x is real; it is a forum claim with no method.
