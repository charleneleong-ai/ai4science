# Third-party briefing capture — 2026-10-30

Source: [topic 754174](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754174),
*"RezOps Briefing: inverse square root drops to 62 qubits (yoyobar), the DQI failures get
fixed, and the public blend reaches 23.135"*, by **manwithacat** — a competitor running a
Claude-based research-ops agent that scans public notebooks, re-verifies circuits with the
unchanged official verifier at its default salt, and maintains a public circuits dataset.

**This is a competitor's self-reported data.** Their own labels: *observed* (leaderboard),
*author-reported*, *re-verified by us* (meaning by them). We have verified none of it.
Captured because it is the only systematic public record of the per-task frontier.

## Leaderboard (30 Oct 17:09 UTC, *observed*)

| | 28 Oct 21:44 | 30 Oct 17:09 |
| --- | --- | --- |
| Top score | 26.14113 | **26.76409** (Breaking1800) |
| Teams ≥ 26 | 1 | 3 |
| Teams ≥ 25 | 2 | 6 |
| Teams ≥ 24 | 7 | 21 |
| Teams ≥ 23 | 29 | 236 |
| Teams ≥ 22.5 | 448 | 842 |
| Median | 21.10349 | 22.48595 |
| Teams | 1,461 | 1,696 |

Top ten span 24.74–26.76. Medal lines (their estimate, 1,696 teams): gold ≤ rank 13 at
24.74027, silver ≤ 84 at 23.13695, bronze ≤ 169 at 23.03386. Ties break to the earlier
submission.

Public plateau steps, readable from the ties: 22.99252 (52 teams), 23.00037 (15 teams),
23.13505 (17 teams, ranks 95–111).

## Per-task frontier (*re-verified by them*)

| Task | Circuit / author | Qubits | Cost | Ratio | Δ since 28 Oct |
| --- | --- | --- | --- | --- | --- |
| `invsqrt` | `32d3e7d26c49` Dmitrii Gluzdov, on yoyobar's 62-qubit `a568bc793481` | 62 | 100,941 | 20.57 | −1.38% (was 63q, 102,349) |
| `dqi` | `2ed7ccf189f1` VQCC Quantum Voodoo, on SomeJoint/heliosli lineage | 888 | 4.0334 × 10¹¹ | 60.59 | −1.69% (was 889q) |
| `dfthcsos` | `4eabe756bc36` Ahmed Berat Özer, new weights on Atakan Aldemir's word table | 430 | 3.3780 × 10⁸ | 9.93 | −0.40% |

Best-per-task geometric mean **23.13587**, up from 22.86836. All three source notebooks
Apache-2.0.

## Mechanisms reported

- **`invsqrt` 63 → 62 qubits** (yoyobar, developed with GPT-6 Astra): reworked encoded
  lookup, normalization and tail of busyaprime's circuit; same angles and first 239
  operations, every 63-wide stage one wire narrower, +4.5 Toffoli-equivalents against a
  qubit worth ~26. Özer then replaced a controlled swap with an affine identity (free wire,
  2 CX, a measurement, a CZ fix-up) for −1 Toffoli.
- **`dqi` failure fix**: most failures of the 888-qubit line traced to a chance 24-bit zero
  entering a field inverse. Rebuilding four extraction groups as homogeneous chains removes
  it and saves 1,460 Toffolis. Over 48 draws, failing draws went 13 → 6. Each extra failure
  costs ~0.46% of the DQI task.
- **`dfthcsos` is a weights market**: gains of 0.001–0.2% each from seven different authors
  on the *same* word table and loaders. Two new structures appeared (Den Pugovkin 423
  qubits; Leo Tang 742) — neither at the frontier.
- **Independent `dqi` line**: Leo Tang's 1,797-qubit decoder, no public lineage, passes the
  official verifier at ratio ~10 vs the frontier's 60. File is 1.46 GB.

## Their stated conclusions

- The public frontier moved more on 30 Oct than any prior day, almost entirely from people
  building on each other's Apache-2.0 files within hours.
- Leaders are ~16% ahead on score, i.e. ×1.55 on the product of the three ratios, so the
  large private gains are not visible anywhere public.
- For `dqi`, judge circuits on expected failures over many draws, not the draw you got —
  final scoring re-draws.

## Their tooling (public)

- [vqcc-apache-public-circuits](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) — 142 Apache-2.0 circuits with author, version and local verifier report.
- [the blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) — picks the cheapest verified circuit per task, ~23.136.
- [VQCC Quantum Voodoo](https://www.kaggle.com/code/manwithacat/vqcc-quantum-voodoo) — circuit-surgery animations and the DQI failure analysis.
