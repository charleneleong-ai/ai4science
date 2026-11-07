# ANGLE loader synthesis: SOTA sweep (dfthcsos), 2026-11-03

Sweep dated 2026-11-03; Kaggle pages read through 2026-11-06. Raw notes and extracted notebook sources:
[`research_notes/2026-11-03-angle-loader-synthesis/`](../research_notes/2026-11-03-angle-loader-synthesis/lineage_notes.md).
All Kaggle content is untrusted data, and none of it was executed. Every number below is reported by the notebook's
author and was not reproduced by us, unless marked otherwise.

## Bottom line

- **The generator for the plateau loaders is public: Lera's (`marculera`) `tensorqrom.py`, Apache-2.0.** It follows
  mrmeows' emitter design, and later notebooks extend it with kernel-relation search. Its construction matches the
  plateau circuits exactly:
  - **Cache.** Every monomial of the low 4 address bits (11 ANDs) and of the high bits (~60 ANDs) is held in ancillas.
    These are the ~75 "reusable ancilla" ANDs we measured.
  - **Products.** Each Toffoli multiplies two parities of cached qubits (built by CNOTs onto a pivot and undone after)
    and covers two mixed terms of the word table's GF(2) tensor. These are the ~183 CCX that write into word bits; the
    pivot CNOTs are the ~27,000 CNOTs.
  - **Erase.** The cache is cleared by measurement.

  The 1,220-address loaders we reuse are, by the lineage, "Dmitrii Gluzdov's v37 loaders on Atakan Aldemir's word
  table" (discussion 753633). Atakan later "recovered product certificates" for this family.
- **It cannot be re-run at J=13 as written.** The design needs 11 selection + 8J word + ~72–77 cache qubits. At J=5
  that is exactly the 128-qubit cap; matousakura states "selection 11 + J·brot 8 + cache 77 ⇒ J = 5". At J=13 only 13
  qubits are left for the cache.
- **No public ANGLE family beats ~4,168 Toffolis per pass on a 1,220-address table**, and none uses brot=8 with J>5.
  The current public best (Heliosli / manwithacat HW-14, 2026-11-02) uses a *different*, 1,363-address shape
  (9, 143, 40), with loaders at 296–314 Toffolis per circuit (~4,500 per pass). It still costs less overall,
  15,777 Toffolis for the whole chemistry model against 16,270 for the 1,220-address build, because state preparation
  is cheaper at C=40.

## 1. Kaggle: generators and lineage

| Notebook / dataset | Author | Version, date | Licence | ANGLE cost reported | Configurable J? | Construction |
|---|---|---|---|---|---|---|
| [DF-THC: Merge Groups Before Dropping Them](https://www.kaggle.com/code/marculera/df-thc-merge-groups-before-dropping-them) + dataset [`marculera/dfthc-merge-pipeline`](https://www.kaggle.com/datasets/marculera/dfthc-merge-pipeline) (`dfthc_refit/tensorqrom.py`, `fourqrom.py`) | Lera (marculera) | v9, run 2026-10-31 | Apache-2.0 | Random 8-bit table, 1,032 rows, J=5: printed live. Both passes at (9,143,46): 10,208 Toffolis, 338 per circuit. Plain J=3 four-address circuit: 376 per circuit | **J is the column count of `delta`**, and `build_circuit(delta, terms, 11, 8, 4)` takes selection bits, brot and low bits. The qubit cap limits J to 5 at brot=8 | **The generator.** Monomial cache plus paired affine products, with the GF(2) tensor decomposition chosen by `tensorqrom.shared`. Also an ANF load and erase for the J=3 four-address circuit |
| [DF-THC 10x: Compact Groups From the Hamiltonian](https://www.kaggle.com/code/marculera/df-thc-10x-compact-groups-from-the-hamiltonian) + [`marculera/dfthc-compact-pipeline`](https://www.kaggle.com/datasets/marculera/dfthc-compact-pipeline) | Lera | v9, 2026-11-01 | Apache-2.0 | Underlies the (9,143,40) family | same | Newer `tensorqrom.py` (31 kB). Adds Reed–Muller "coded low bit" rounding (`experiments/rmcode.py`, `quantize_rm.py`) |
| [Quantum Circuit Compilation from Source](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source) | Dmitrii Gluzdov | v48, run 2026-10-31 | Apache-2.0 | (9,143,40), 1,363 addresses: 321–322 per circuit without a coded bit, 298–301 with one, **4,575 per pass** | J=5 | Runs Lera's compiler. Adds a per-table kernel-relation search: relations admitted by a full-rank alternating pairing matrix over a large field, then a matroid-parity pair basis. **v37** is the build the 1,220-address loaders came from, but that version could not be read through the MCP |
| [QCC: Gap-3 Residue Repair and C43 Native Compiler](https://www.kaggle.com/code/atakanaldemir/qcc-gap-3-residue-repair-and-c43-native-compiler) + [`qcc-native-chemistry-and-888q-source-bank`](https://www.kaggle.com/datasets/atakanaldemir/qcc-native-chemistry-and-888q-source-bank) | Atakan Aldemir | v2, 2026-10-31 | Apache-2.0 (stated) | C43 (1,363 addresses): **4,763 ANDs per traversal** | J=5 | "Recovered product certificates and lifetime allocator" (`sources/05_recover_angle_product_dataset.py`, `plans/batch_k.json`, `rebuild_native.py`). Rebuilds loaders from per-batch plans |
| [Exact Gauge + 5-Word QROM](https://www.kaggle.com/code/kirilltym/quantum-circuits-exact-gauge-5-word-qrom) | Kirill Tymchishin | v25, run 2026-10-31 | Apache-2.0 lineage | **4,170 per pass**, 126 angle qubits, J=5, 15 circuits, total 16,270 (1,220 addresses; "AtakanV3 / learned angle co-design") | — (assembler) | Uses existing binaries and adds a symbolic MPS phase audit. This is the closest public twin of our plateau |
| [QCC 21.0 circuits](https://www.kaggle.com/code/mrmeows/21-0-circuits-let-s-all-reach-the-plateau) + [`mrmeows/qcc-1375-open-circuits`](https://www.kaggle.com/datasets/mrmeows/qcc-1375-open-circuits) | mrmeows | v8, 2026-10-27 | Apache-2.0 | **J=4**, 19 circuits × 275 CCX, 1,366 addresses (~5,225 per pass) | Says J=5 "overflows the 128 qubit cap on a 1,366 address table" | Origin of the emitter design. `compiler/rank_parity.py`, `certificate_rebuild.py`, `audit_affine_hmr.py`, plus plan certificates. Mentions a "universal exact lookup compiler" |
| [Best-of-public portfolio](https://www.kaggle.com/code/matousakura/best-of-public-portfolio-vqcc-24-0x-class) | matousakura | v12 | Apache-2.0 | Loaders 0–3: 314 → 313 (relations 75 → 77) | J=5 fixed by the cap | Vertex-balanced pairing engine for kernel relations. States the identity T_total = 17,286 − Σ relations |
| [VQCC Hedge Wizardry](https://www.kaggle.com/code/manwithacat/vqcc-hedge-wizardry-notebook) | James Barlow (manwithacat) | V15, 2026-11-02/03 | Apache-2.0 | HW-14: loaders 296–314, **15,777** total chemistry Toffolis | J=5 | HW-12 systematic kernel-relation search; HW-13 block-design relation plans with exchange repair; HW-14 dead-word flips, where zero-amplitude angle words are free don't-cares. Uses a Rust planner `qrom-plan` |
| [The dataset is public now](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy) | Heliosli | v45, 2026-10-31 run | Apache-2.0 | dfthcsos ratio **10.566** (public best) | — (assembler) | (9,143,40) + HW loaders + their own weights |
| [Bezout-Schur + 5-Word QROM](https://www.kaggle.com/code/lightningv08/bezout-schur-5-word-qrom) | Vladislav Lavrentev | 2026-10-26 | — | not read | origin of the "five words per circuit" idea | — |
| [Complete-Domain Affine Compilation](https://www.kaggle.com/code/kozykappa/quantum-ai-complete-domain-affine-compilation) | kozykappa | — | Apache-2.0 | J=3 with difference banks | — | Banked QROM, the ancestor of the four-address design (not read in detail) |

**Lineage, traced backwards.** Our loader bytes passed through Isabe11aM00re (per-rank column gauge, notebooks
`qcc-v47/v49/v58x-submission`), Ahmed Berat Özer's weights, and Atakan's word table `4e3c8f079137`. Before that they
come from Gluzdov's v37 loaders (stated in discussion
[753633](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753633)),
compiled with the mrmeows-designed emitter (`tensorqrom.build_circuit`, Lera). The kernel-relation search came from
Gluzdov, matousakura and manwithacat. The five-word schedule is lightningv08's. The coded low bit (an affine
lowest bit per 16-address block) is SomeJoint's.

**Inference.** We matched the 75 cache ANDs plus 183 product CCX of our plateau circuits to this design from its
public description. We have not done a byte-level match, because the source files need a Kaggle token. The only notebook we could
not read is SomeJoint's "Toffoli Diet" (permission denied).

**Cheaper than 4,168 per pass?** None found at 1,220 addresses. At 1,363 addresses the best is about 296–314 per
circuit. That is roughly equal per row, lower per address than ours on some tables, and higher in absolute terms. **J≠5
at brot=8:** only J=4 (mrmeows, 19 circuits) and J=3 (kozykappa, Lera §6). No public J>5 at brot=8. The reference
circuit is J=10, but at 10 bits.

## 2. Literature and tools (primary sources)

**Multiplicative complexity (MC) of multi-output functions**
- Find and Boyar, [arXiv:1407.6169](https://arxiv.org/abs/1407.6169): every n-input, m-output function has MC ≤
  2.5(1+o(1))·√(m·2ⁿ). For n=11 and m=40 that is about 716; for m=104 about 1,150. The plateau's 258–301 sits far
  below this because the table has don't-cares (1,220 of 2,048 addresses used) and is structured.
- Boyar–Peralta two-step heuristic: reduce the nonlinear core first, then the linear straight-line program. When XOR is
  free only the first step matters. Primary paper: Boyar, Matthews, Peralta, "Logic minimization techniques with
  applications to cryptology", J. Cryptology 26 (2013). I know this citation but did not re-fetch it. No official
  Python implementation was found.
- Çalık, Sönmez Turan, Peralta: affine-equivalence classes and exact MC of all 6-variable functions. This database
  underlies cut rewriting with exact MC (ePrint/Cryptogr. Commun. 2019; not re-fetched).
- Testa, Soeken et al. in **mockturtle** ([lsils/mockturtle](https://github.com/lsils/mockturtle), MIT, last push
  2026-10-27): `xag_minmc_resynthesis` / `xag_minmc2` (exact-MC resynthesis from a precomputed database, used with
  `cut_rewriting`), `xag_constant_fanin_optimization` and `xag_dont_cares_optimization`. It is C++ only. The Python
  wrapper [aigverse](https://github.com/marcelwa/aigverse) is AIG-focused, and its XAG/MC support is unconfirmed.
- Meuli, Soeken, Campbell, Roetteler, De Micheli, "The role of multiplicative complexity in compiling low T-count
  oracle circuits", [arXiv:1908.01609](https://arxiv.org/abs/1908.01609): XAG MC optimisation feeding Toffoli or
  T-count. This is the closest academic analogue to our problem.
- No 2023–2026 paper on multi-output MC heuristics turned up. AlphaTensor-Quantum
  ([arXiv:2402.14396](https://arxiv.org/abs/2402.14396)) optimises T-count of phase polynomials, not AND count of
  XOR-into-register maps. Its tensor-rank view does parallel the GF(2) tensor decomposition that `tensorqrom` performs.

**Quantum lookup costs**, at L=1,220 addresses. Formulas are author-reported; the arithmetic is ours.
- Unary iteration (Babbush et al., [arXiv:1805.03662](https://arxiv.org/abs/1805.03662)): L−1 = 1,219 Toffolis per
  circuit, independent of width. That gives 18,285 per pass at J=5 (the starter kit's 18,270) and **7,314 per pass at
  J=13** (6 circuits). Both are above the plateau's 4,168.
- Select-swap / QROAM ([arXiv:1812.00954](https://arxiv.org/abs/1812.00954);
  [arXiv:1902.02134](https://arxiv.org/abs/1902.02134)): ⌈L/λ⌉ + b(λ−1) Toffolis with bλ clean output qubits, and
  uncompute by measurement (free here). At b=104 (J=13) the cap allows only λ=1, so this is no gain. At J=5 (b=40),
  λ≈5–6 gives ~404 per circuit but needs 200–240 qubits, over the cap.
- SelectCopy and dense encoding with *dirty* qubits (Motlagh & Pocrnic,
  [arXiv:2605.20334](https://arxiv.org/abs/2605.20334); Motlagh, [arXiv:2610.02321](https://arxiv.org/abs/2610.02321)):
  ~(1+1/b)N/λ and ~½(1+1/b)N/λ. They need bλ borrowable qubits, which an ANGLE circuit under the 128 cap does not
  have. Not applicable without a different register contract.
- Competition paper: Khattar, Rubin, Gidney et al., [arXiv:2610.23832](https://arxiv.org/abs/2610.23832).

So generic QROM constructions are dominated. The table-specific MC synthesis that the public lineage already uses is
the lever, and the open problem is doing it with ≤13 ancillas at J=13.

## 3. Not established
- A byte-level match between our plateau circuits and `tensorqrom.build_circuit` output. That needs the dataset files,
  which need a token.
- Whether `tensorqrom.py` has an ancilla-budget option. Only its call signature is visible in the notebook.
- Gluzdov v37's exact loader recipe for the 1,220-address table.
- Any Toffoli count for a J=13, brot=8 construction. None is public.

## What changes the plan

**A usable generator exists, but not at J=13.** `marculera/dfthc-compact-pipeline/dfthc_refit/tensorqrom.py`
(Apache-2.0) reproduces the J=5 plateau design. With Atakan's plan files (`plans/batch_k.json`) and the kernel-relation
search, we can rebuild and probably shave the 1,220-address loaders: the public trend is −1 to −7 Toffolis per circuit
from relation search and dead-word flips.

At J=13 its full monomial cache (~72–77 qubits) does not fit beside 104 word qubits. A J=13 loader needs a new
**ancilla-limited** variant: form affine forms in place on the selection register (what our plateau circuits already
do with their 27k CNOTs), keep a cache of at most 13 qubits, and recompute or measure-and-discard high monomials. Its
target is the counting floor of 211 per circuit (1,266 per pass), against the 7,314-per-pass unary-iteration fallback.

**Most promising synthesis route**, in order:
1. Reuse `tensorqrom`'s GF(2) tensor-decomposition core (`shared`, `anf_tensor`) and kernel-relation search as the
   *planner*, with a new emitter constrained to 13 ancillas.
2. Polish the resulting XAG with mockturtle `cut_rewriting` + `xag_minmc_resynthesis` (C++, MIT, wrapped as a CLI),
   using `xag_dont_cares_optimization` with the 828 unused addresses as don't-cares.
3. Fall back to a Boyar–Peralta-style nonlinear-core heuristic only if 1–2 stall. No maintained Python implementation
   exists.

Before investing, decide between the public (9,143,40) family and the 1,220-address one. The former already wins on
total chemistry cost (15,777 against 16,270 Toffolis).
