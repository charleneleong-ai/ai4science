# CASMI 2026 spike — findings and whether to commit

**Verdict:** worth entering *only* if the goal is a demonstrable methodology piece. It is not a
realistic medal push at ~9 weeks part-time against 2,520 teams, because the one thing standing
between a working ranker and a real submission — a candidate pool that covers the test
chemistry — is unsolved and is the largest remaining cost.

Grounded 2026-10-07. [Competition](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra):
predict up to 25 ranked SMILES per molecule from LC-MS/MS, scored MRR@25 on the
tautomer-canonical InChIKey connectivity block. Entry closes 2026-12-07.

## What was measured

| | recall | MRR@25 (mass floor) | MRR@25 (fragments) | median candidates |
|---|---|---|---|---|
| train pool (276k) | **0.9959** | 0.4119 | **0.7011** | 46 |
| COCONUT (474k) | **0.0806** | 0.0363 | — | 10 |

500 held-out molecules, split on `inchikey14`, 5 ppm window. Recall is the ceiling: ranking
cannot exceed it.

## Three findings, in order of consequence

### 1. The public leaderboard measures an exact join, not elucidation

Every one of the 400 public test molecules has a train spectrum at peak-cosine ≥0.95, median
exactly **1.000**. Confirmed twice by independent methods: full-scan cosine (identical peak
counts and m/z, all from `enveda-180`), and a hash join on `(adduct, precursor, n_peaks)` that
resolved **400/400**.

So the public test spectra are duplicated training rows. This explains the ~0.87 MRR figures in
public repos and inverts their meaning — they are evidence the join works, not that a model
does. **Every public-LB comparison is uninformative**, which is why evaluation lives in
[`harness.py`](../src/casmi/harness.py) instead.

This is the one genuine competitive edge the spike found: teams tuning against the public LB are
optimising a lookup.

### 2. No public natural-product database covers the test chemistry

COCONUT (CC0, 474k structures after collapsing to `inchikey14`) covers **0 of 400** test
molecules, and only 8% of held-out train molecules.

The reason is that the test set is not natural products. All 400 resolve to a single library,
`enveda-180`; **36% contain a halogen**; and the structures are plainly medicinal chemistry —
brominated cyclopropyls, triazoles, spiro-morpholines. Despite the CASMI name and Enveda's
natural-products business, the public test molecules are synthetic drug-like compounds.

**Consequence:** a real submission needs PubChem-scale retrieval, not COCONUT. That is the
dominant unsolved cost, and it will move the numbers — windows grow from ~46 to plausibly
thousands.

### 3. Ranking is the constraint, and combinatorics gets surprisingly far

Recall is ~1.0 *when the pool contains the answer*, so retrieval is not the bottleneck. The real
problem is that a 5 ppm window holds ~46 candidates of which **~29 share the truth's exact
molecular formula** — mass accuracy is blind among them.

Fragment explainability ([`fragments.py`](../src/casmi/fragments.py)) takes MRR@25 from the
**0.4119** mass floor to **0.7011**, with no training, no spectral library and no model weights.

## What the 0.70 is and is not

It is ranking performance *within train's own chemical space*, on a pool that holds the answers
because it is the same file. It is **not** comparable to any leaderboard number, and it will fall
against a pool that covers `enveda-180`-like space at PubChem scale.

## Premises this spike overturned

Three things believed at the start proved wrong on contact with data, which is the main argument
for having run it at all:

| believed | measured |
|---|---|
| retrieval + rerank is the research frontier | de-novo generation is; graph diffusion broke 0% exact-match on MassSpecGym first |
| ICEBERG is right *because* CASMI is natural products | the test set is synthetic drug-like; the natural-product argument does not apply |
| fragment enumeration cannot separate branching isomers | it can — that was a hydrogen-accounting bug in RDKit's fragmenter |

Two were my own errors, found by a result being implausibly good or a fixture failing.

## If continuing, in order

1. **PubChem-scale pool.** Everything else is uninterpretable without it; the 0.70 has no meaning
   against a pool that cannot contain the answers.
2. **Then** a learned reranker. ICEBERG remains plausible but its stated natural-product
   advantage is now irrelevant here, and its pretrained weights are described as NIST'23 — a
   licensed library, so usability needs checking before any plan depends on it.
3. Private-set composition is unknowable directly, but is presumably drawn like the public set.

## Reusable regardless of the decision

The harness is the asset: molecule-disjoint splitting on the metric's own equivalence class, and
recall reported **separately** from ranking so it is visible which one binds. Both leakage
findings came from that separation — the first because recall was impossibly perfect, the second
because a ranker using no spectral information scored 0.94.
