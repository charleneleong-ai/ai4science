# Raw notes: VQCC frontier delta sweep, 2026-11-12

Swept 2026-11-12, roughly 10:00–11:00 UTC. Every number below is labelled with where it came from:

- **observed**: the Kaggle CLI or a Kaggle page;
- **author-reported**: the poster's own claim;
- **ours**: arithmetic done for this sweep.

Notebooks were pulled as `.ipynb` text with `kaggle kernels pull -m` into
`~/quantum-compilation-data/_research_nb/2026-11-12-sweep/<owner>_<slug>/`. None was executed. Forum threads were read
through the Kaggle MCP topic tools, which work on public pages without auth. WebFetch of the discussion list returns
only the page title, because the page is rendered by JavaScript.

## Competition pages (kaggle CLI `competitions pages list --content`)

- **Prizes.**
  - Leaderboard prizes: $50k, $25k, $10k, then $5k each for 4th to 6th.
  - The **Mid-Competition Milestone Prize ($10k)** goes to the highest public LB score at 2026-11-20 23:59 UTC, and
    the winner must publish a notebook.
  - The **Longest Leader Prize** is $10k, counted up to 2027-01-04.
  - The **Single Task Performance Prizes** are $10k per task, for "the selected submission with the best score for
    that specific task". This is judged at the end of the competition.
  - **No winners are published yet. No per-task leaderboard exists** (requested in thread 751968, with no host
    response).
- **Timeline.** Final submission deadline 2027-01-20.
- **Evaluation.**
  - invsqrt cost is Q·C, with C = Toffoli + T/4.
  - dqi cost is R(ε₊)·Q²·T.
  - dfthcsos cost is λ_eff·T·Q. Here T and Q cover the whole block encoding, and an ANGLE circuit wider than 128
    qubits is rejected.
- **Papers linked from the pages:**
  - arXiv:2610.23832, the task paper: Khattar, Rubin, Gidney, Broughton, Low, Babbush, *Leveraging verifiable
    structure in quantum circuits to find more efficient compilations with agentic AI*, Oct 2026;
  - arXiv:2510.10967, the DQI circuits paper;
  - arXiv:2603.28846, not opened.

## Leaderboard (kaggle CLI, observed 2026-11-12 ~10:00 UTC)

| Rank | Team | Score | Submitted (UTC) |
| --- | --- | --- | --- |
| 1 | Schrödinger’s cat lived | 29.52446 | 2026-11-12 07:13 |
| 2 | Entangled Team | 29.45296 | 2026-11-11 22:02 |
| 3 | yicheng guang | 28.54171 | 2026-11-12 07:40 |
| 4 | taiseiu | 28.44097 | 2026-11-07 |
| 5 | Some Joint | 28.37025 | 2026-11-12 07:55 |
| 6 | yu4u | 28.14827 | |
| 7 | Breaking1800 | 27.99617 | 2026-11-01 |
| 8 | CroDoc & Pavel | 27.78422 | |
| 9 | Diana | 27.66022 | |
| 10 | Yan Zhou | 27.45652 | |
| 20 | | 26.77669 | |

## Public frontier (notebooks, observed text and author-reported verifier outputs)

- **attabeezy/gq-portfolio-frontier** (<https://www.kaggle.com/code/attabeezy/gq-portfolio-frontier>).
  - Claims LB **25.40482**, submission 58260294, graded 2026-11-12.
  - Its per-task costs are the authors' receipts: invsqrt 93,555.26; DQI 373,534,939,340.64; chemistry
    297,237,077.47.
  - Files:
    - invsqrt `9c708f7b`;
    - DQI `eae396ba`;
    - chemistry: maoxin123's V22 4,089-Toffoli loaders plus prvsiyan's V40 weights `f21b8a87`.
- **kirilltym/quantum-circuits-retained-product-audit**. Independently re-verifies the same portfolio:
  - invsqrt 93,555.31;
  - DQI 373,534,984,252.11, with R = 1.191260465 and 0 value / 0 phase errors;
  - chemistry 297,237,077.47, with loaders at 4,089 per pass and 14,840 Toffolis in total at 427 qubits.
- **heliosli/the-dataset-is-public-now-enjoy**. Per-task ratios: invsqrt 22.1984, DQI 65.4275, dfthcsos 11.289. Same
  files.
- **ozertuu/25-39570-sota-vqcc-fast-runner-portfolio**. 25.39570, with invsqrt 22.19837, DQI 65.42752 and chemistry
  11.27715.
- **apogean/vqcc-verified-portfolio-reproducible-assembly**. Variant D = 25.39570 (ref 58259265). It also notes a
  "salt lottery" on DQI: one local salt gives one value_error.
- **maoxin123/qcc-v22-4089-toffoli-chemistry-loaders**. Per-loader Toffolis per pass:
  - k0: 292;
  - k1–k3: 287 each;
  - k4–k8: 268 each;
  - k9: 265;
  - k10–k13: 266 each;
  - k14: 267;
  - total 4,089. The address cache (75–80 Toffolis) is unchanged.
  - The author says "applies R only to qubits that are provably |0>".
- **atakanaldemir/qcc-shared-products-and-controlled-refits**. DQI phase block: `(-1)^popcount(m & inv_GF64(x))`
  drops from 6 to 5 Toffolis per nonzero mask, and V15 brings it to 4.
  - There are 28 calls, saving 27.56 mean Toffolis.
  - Static T falls from 404,889 to 403,125, at Q = 892.
- No public notebook was found with a per-task ratio above invsqrt 22.20, DQI 65.43 or chemistry 11.29. Searched by
  vote count, by run date (top 40 of each) and by keyword: dqi, invsqrt, chemistry, dfthc, thc, loader.

## Forum

- **756998, "Anyone got higher in invsqrt?"** (mikelou1, competition_ranking 31 at the time, posted 2026-11-08).
  - "I'm currently scoring **25.1x in invsqrt** task and using public baselines for everything else."
    **Author-reported; method not disclosed.** That is +13.4% over our 22.1255.
  - Reply from Adarsh (rank 301): "Mine is barely above public baseline, 22.17."
- **757915, "Fable 5.2 is an expensive beast"** (pjmathematician, **competition_ranking 1**).
  - The post covers model and token cost only, with no technique. Fable is used for planning and Kimi for coding.
  - GrantPrater (unranked): "constraining the search space … representing the problem in [the right] way".
  - yuanzhe zhou (GM): "GPU … searching very fast, but no improvement".
- **757030, "Open problems we are stuck on"** (Rustam Bazarbayev, rank 1481).
  1. Normalising barrel shifter. Public circuits sit in the high 40s of CSWAPs. A counting lower bound is in the
     high 30s. The SAT solver did not resolve w = 21.
  2. Is 62 qubits a real floor for invsqrt? The liveness plateau is wide. One qubit is about 25 Toffolis under
     Q·C.
  3. Cheap uncompute for chains of in-place GF(2^m) multiplications, which matters for DQI.
  4. GF(2^12) multiplication needs 42 bilinear products (μ_4(6) = 14 is proven, thread 756476).
  5. Chemistry loaders are near "~2·sqrt(m·L)". Loader-aware rounding is unexplored (author's claim).
- **753633, "dfthcsos: why λ_eff can't go much below ~46.8".** The bound is λ_eff ≥ 46.79 at the gap floor, for any
  factorization. **Our λ_eff 46.89 is therefore within 0.2% of it.**
- **753527, "What is one qubit worth?"**
  - invsqrt: 1 qubit = 1.61% of cost (+0.54% score).
  - DQI: 1 qubit = 0.23%.
  - Chemistry: 1 qubit = 0.23%.
- **755302, "Is resetting workspace with R allowed?"**
  - No host ruling.
  - One user claims that in public DQI decoders the uncompute half is "roughly half of the mean Toffoli count", and
    in invsqrt "roughly a quarter" (author-reported).
  - A reply tested input-dependent R, and **the verifier catches it**: a random phase on |1⟩ fails with a 0.5-turn
    phase error. So HMR plus phase fix-up (MBU) is the only cheap uncompute.
- **757522, chemistry ANGLE contract.** Some accepted ANGLE packets are correct only on reachable word states, with a
  reset word bit and a measurement-dependent phase. The question is whether strict Eq. (72) (any input state) or
  reachable-states-only is the rule. **No host answer**, as of 2026-11-12 03:23 UTC. This is an eligibility risk for
  the public chemistry frontier.
- **758095 / 757024.** Milestone prize mechanics are unresolved: must the top team's scoring notebook itself be
  public?
- **758129, RezOps briefing of 2026-11-11** (manwithacat).
  - Leader 29.30100 on 11 Nov 17:34 UTC.
  - 39 teams at or above 26.
  - Gold line (rank 15) 26.79893.
  - The public best-per-task projection was 25.35895.
  - "The 16% gap to the leaders is unchanged."
  - ALS and Levenberg–Marquardt refits of the chemistry table lowered the two-body error only about 0.1 mHa: "the fit
    floor at this table shape looks like a property of the shape".
- **755038, RezOps briefing of 2026-11-02.**
  - Public per task then: invsqrt 20.65 (100,568, 62q), DQI 64.72 (892q, 48-draw mean), chemistry 10.27 (428q,
    (9,143,40) table).
  - Best-per-task GM 23.94157. Leader Breaking1800 at 27.99617.
  - "The leaders are about 17% ahead of the public frontier."
  - "DQI is now won by exact bookkeeping at fixed width."
  - Yoyobar's decoder: mean Toffolis per shot went from 399,970 to 398,139.
- **753451.** Reachable invsqrt inputs are 76.5% of 1..3,133,452 (sums of three squares with coordinates ≤ 1022). The
  author's estimate is that exploiting the gaps is worth about 0.08% of invsqrt cost.
- **752757** (mrmeows, 2026-10-28). Promised a "25x" public release. As of 2026-11-11 it has not appeared.
- **751968.** A per-task leaderboard was requested, and none has been granted.

## The DQI register contract (from the shipped verifier `verify_submission_dqi.py`, read locally)

- The instance is m = 4095, n = 70, b = 12, with l = 35 correctable errors.
- **Registers S0..S69 are 70 × 12 = 840 quantum syndrome qubits, and they are the input.** E0..E4094 are classical.
  - Every other qubit must return to |0⟩.
  - So **Q ≥ 840 is a hard floor**: today's 892 = 840 + 52 workspace.
- Width headroom is at most (892/840)² = 1.128 on DQI cost, and even that needs zero workspace.
- `Z_POW` is not allowed in DQI.

## Papers (arXiv abstract pages fetched)

- **arXiv:2610.23832**, the task paper. Its applications:
  - 1/√x;
  - Reed–Solomon decoding for DQI;
  - QROM plus controlled Givens rotations for chemistry.

  It verifies by fuzz testing over classical inputs. The abstract gives no reference costs.
- **arXiv:2510.10967**, Khattar, Shutty, Gidney et al., *Verifiable Quantum Advantage via Optimized DQI Circuits*
  (Oct 2025).
  - The Extended Euclidean Algorithm is the bottleneck.
  - It uses implicit Bézout-coefficient access and in-place architectures.
  - Its leading-order space is the "theoretical minimum of 2nb qubits".
  - It estimates about 5.72M Toffolis for intractable instances.
- **arXiv:2611.09822**, DeepMind et al., *Discovering formally verified quantum programs with AI* (Nov 2026).
  - A Lean library with agent loops.
  - The abstract reports 11× on "simulations relevant to fusion", 290× on quantum arithmetic, and 8× on
    photochemistry, all against the paper's own baselines. These are not this competition's references, and the
    mapping to the tasks is not established.
- **arXiv:2407.17966**, Khattar & Gidney, *Rise of conditionally clean ancillae*, Quantum 9, 1752 (2025). Already
  used publicly in DQI (RezOps HW-09: −294 Toffolis per shot).
- **arXiv:2308.12352**, Rubin et al., stopping power. The source of the invsqrt task (cited in thread 753451).
