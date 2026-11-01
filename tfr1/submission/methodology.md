# A pH-conditional, mouse-cross-reactive TfR1 binder

20 single-chain *de novo* miniproteins, 46–49 aa, 2 scaffolds × 10 histidine variants each.
Nothing is measured in a wet lab. The three objectives are ranked pH-selectivity, mouse
cross-reactivity, affinity, and **two design rules were measured and found to be wrong** partway
through — those two corrections are the most useful thing in this entry. The first is about
chemistry: the switch is set by burial, not by pairing a histidine with a receptor carboxylate.
The second is about where the designs actually sit, and it cost this entry its best epitope.

## The epitope was chosen by measurement, and the recommended one lost

Four candidate epitopes were scored on human/mouse conservation, solvent exposure, N-glycan
clearance and distance from the transferrin site, from structure and sequence rather than from
the brief:

| epitope | human/mouse identity | carboxylates facing it | transferrin clearance |
| --- | --- | --- | --- |
| apical — the brief's recommendation | 47% (VP1u core 22%) | E350, E369 | 46.3 Å |
| 529 | **90%** (only A529→I) | D525, E533 | 12.3 Å |
| 605 | 83% | 4: D604 E606 E664 D667 | 6.3 Å |
| 647 | **100%** | D757, E759 | **0.0 Å — it is the transferrin site** |

The most conserved surface on this receptor is the one transferrin already occupies.
Conservation and function coincide here, which is the structural reason cross-reactive TfR1
binders are hard. Three epitopes were carried rather than betting the entry on one reading of
that trade-off: apical (the recommendation), 529 (`conserved`) and 605 (`carboxylate`).

**The recommended epitope had the lowest interface-confidence ceiling of the three** — ipTM
median 0.130 and best 0.188 over its 100 designs, against 0.153 / 0.366 for 529 and 0.110 /
0.387 for 605 — and contributed none of the 8 scaffolds that reached the switch stage.

The scan's two binding gates (contact patch ≥ 10 residues, transferrin clearance ≥ 6 Å) were
**loosened after seeing the first result** and are declared here as a post-hoc change. The size
gate was landing on the mode of the patch-size distribution, so it was deciding the campaign by
two residues.

### Calibrated against the organisers' own figures first

The scan reproduces two numbers stated in the brief, independently from structure and from
sequence: the nine apical VP1u contacts in **11RN** come out as exactly 202, 204, 208, 209, 210,
211, 212, 348, 371, of which two are conserved; and human/mouse identity over the ectodomain
121–760 comes out at **78.0%**. Four of the brief's six named divergent residues land exactly;
210 and 211 differ only in where the gap falls inside the 210–214 indel, which is
alignment-ambiguous and leaves both non-conserved either way.

## The switch is set by burial, not by carboxylate pairing

A histidine that titrates in the right pH window does nothing to binding if its pKa is the same
bound and free — an equilibrium identical in both states cancels. The quantity that matters is
the **bound-vs-free linkage, Δ = pKa(complex) − pKa(free)**, which must be negative for the
binder to release at pH 6.0.

Each epitope here was originally chosen partly for the receptor carboxylates facing it, on the
reasoning that an engineered histidine titrates against them at the interface. **That reasoning
was measured and it is wrong.** Over 114 readings on 44 designs, Δ tracks how deeply binding
buries the histidine — Spearman ρ **+0.74** (PROPKA) and **+0.51** (pKAI) against distance to the
target — and not proximity to a carboxylate, where ρ is **+0.12** and **−0.27** with the two
methods disagreeing on the sign. Buried within 4 Å the median Δ is −0.98; further out, −0.04.
Within 6 Å of an acid it is −0.225 against −0.200 further away, which is no difference.

The mechanism is desolvation, not pairing. Burial penalises protonation in the complex but not in
the free binder, lowering the complex pKa and giving Δ < 0. A histidine beside a receptor
carboxylate is the opposite case: the acid stabilises the protonated form in the complex, raising
its pKa, so Δ > 0 and the binder holds on **harder** at pH 6.0. Three separate positions in this
campaign read Δ > 0 for exactly that reason.

Reproduced independently on this target over 32 histidine variants of the 3DS18 scaffold, each on
three crystallographically independent copies of the complex: Δ runs −3.0 at 0.9 Å of burial,
−2.9 at 1.9 Å, −1.1 at 2.3 Å, −0.2 at 2.9 Å, and **−0.02 for the native His69 at 11 Å in all 32
variants** — an internal null that the method returns correctly.

Substitution positions in the submitted designs are therefore chosen on buried interface area.
The epitope named `carboxylate` keeps that name only because the generation specs, design names
and recorded results are all written in terms of it; it is retained on its merits as the
second-most conserved epitope (83%), not on its carboxylate count.

## Two pKa methods, and a position counts only when both agree

Every candidate position is read by **PROPKA 3.5.1** and **pKAI 1.2.0**, which share no method.
A position contributes to a design's switch only when **both** methods put Δ below the coupling
floor, and it then contributes the **smaller** of the two estimates. The per-design switch is the
sum over surviving positions, reported as log₁₀ fold weaker binding at pH 6.0.

The submitted twenty carry 1–4 such positions (median 2) and switches of **0.42–1.46 log units**,
i.e. 2.6–29× weaker at pH 6.0. None reaches the 2 log units that "no detectable binding at pH
6.0" would require.

## Epitope engagement is a property of the scaffold, and it was measured late

Every design was named for the epitope it was generated against. Whether it *touches* that
epitope is a separate question, and it was asked only after the first selection had been made.
Each design's own complex was re-read and the binder's contacts intersected with the five
hotspot positions nominated for its epitope, over all 256 histidine variants:

| scaffold | epitope | variants | hotspots touched, of 5 | distance to the nearest one |
| --- | --- | --- | --- | --- |
| `carboxylate_s_37` | 605 | 32 | mean **3.50**, max 4, never 0 | 2.2–3.0 Å |
| `carboxylate_s_10` | 605 | 32 | mean **2.00**, max 2, never 0 | 1.2–3.4 Å |
| six `conserved` scaffolds | 529 | 192 | **0 on every one of the 192** | **5.7–18.8 Å** |

**192 of 256 variants — 75% — touch none of the hotspots they were generated against**, and the
distance says they are not near-misses: the closest `conserved` scaffold sits 5.7 Å off its
nominated patch and the furthest 18.8 Å, while both `carboxylate` scaffolds are in contact. The
split is exact: every variant of the two `carboxylate` scaffolds touches at least two, every
variant of the six `conserved` scaffolds touches none. Engagement does not vary within a
scaffold, because the histidine substitutions do not move the backbone. It is a property of the
backbone the generator produced, not of the variant built on it.

The consequence for this entry is direct. The 529 patch was chosen for its 90% human/mouse
identity, the best of any epitope clear of the transferrin site, and the generator produced
nothing that lands on it. Ten designs in the first selection came from those scaffolds. They
were replaced.

**This gate is post-hoc and is declared as such**: a minimum of one nominated hotspot, applied
after the engagement measurement was read, and off by default in the selection code so that the
pre-registered selection stays reproducible. It is a filter chosen after seeing the designs, and
the funnel below shows what it removed.

### An independent model, and what it says about the scaffold that was dropped

Each of the 26 designs in the pre-registered selection was re-folded with **Chai-1**, a different
model family from the BoltzGen/Boltz-2 weights that both generated and ranked them, at 5 samples
per design — 130 poses. Because that selection predates the engagement gate, **only 10 of the 20
submitted designs carry an independent co-fold read**; the other 10 are ranked on the generator's
pose alone.

**The co-fold is reported as calibration, not used as a gate**, and the first thing it calibrates
is the score itself. The published 3DS18 binder — a real crystallographic TfR1 binder — was run
through the identical pipeline: ipTM **0.296** (0.205–0.326 over 5 samples) against the
639-residue ectodomain, and **0.190** (0.186–0.201) against a 90-residue carved patch of the kind
the designs are scored on. The usual 0.8 threshold is unreachable here by a binder that
demonstrably works. The bar on this target is 0.19–0.30.

Against that bar the 130 design poses span **0.081–0.388**, and the spread is not where the
switch measurement would predict. The three highest-scoring designs in the whole campaign —
0.388, 0.315 and 0.285, two of them above the native control's mean and all three above the
carved control — come from `conserved_18`, one of the scaffolds the engagement gate removed.
They also have by far the most reproducible poses: mean pairwise interface overlap 0.51–0.67
across Chai's 5 samples, against 0.12–0.24 for everything else.

That looked at first like the gate discarding the best designs. Measuring the distance from each
pose to the nearest nominated hotspot says otherwise:

| scaffold | co-fold ipTM | distance to nearest nominated hotspot | pose agreement |
| --- | --- | --- | --- |
| `conserved_18` (dropped) | **0.135–0.388** | **8.5–14.9 Å** | 0.51–0.67 |
| `carboxylate_37` / `_10` (submitted) | 0.103–0.123 | 3.1–10.8 Å | 0.12–0.24 |
| `conserved_45` (dropped) | 0.088–0.102 | 4.7–7.2 Å | 0.16–0.20 |

`conserved_18` scores well because Chai is confident and consistent about an interface **on a
different face of the carved patch**, 8.5–14.9 Å from the epitope these designs were generated
against. A high interface score on the wrong surface is not evidence for the objective, so the
gate's verdict stands and the score is reported rather than acted on.

The same table is the honest limit of the gate. On the generator's own poses the two
`carboxylate` scaffolds touch 2.0 and 3.5 of 5 hotspots; on Chai's poses the submitted designs
touch 0.2–1.4 and sit 3.1–10.8 Å away. **The engagement separation is clear in the model that
generated the designs and much weaker in the independent one.** What survives both readings is
only the ordering — the submitted scaffolds sit closer to the nominated patch than the dropped
ones in both — not the claim that they land on it.

## Mouse cross-reactivity is carried by the epitope, and this entry lost its best one

The second-ranked objective is carried entirely by where the binder sits. All twenty submitted
designs bind the 605 patch, **83%** identical between human and mouse. No design is individually
optimised for mouse, and none is predicted against a mouse structure, which does not exist for
this receptor.

The entry is weaker on this objective than it was designed to be. Half the slots were meant for
the 529 patch at 90% identity; the engagement measurement removed them. 83% is still well above
the brief's recommended apical site (47%, and 22% over the VP1u core), but the hedge across two
epitopes is gone and every design now shares one surface — so a single wrong call about that
surface takes the whole entry with it.

## Funnel

```
300 designs   3 epitopes x 2 length bands x 50
226 pass composition triage (75%)   <- degeneracy, free cysteine, Glu-fraction
  8 taken forward as scaffolds
256 histidine variants of those 8 (32 each)
192 clear the novelty floor          <- 64 dropped: 2 scaffolds score novelty 2
                                        (fold TM 0.94 / 0.90 at 95-100% coverage)
192 parent passes composition triage <- 0 dropped, already applied at the scaffold stage
128 parent has no unpaired cysteine  <- 64 dropped: 2 scaffolds
113 both pKa methods agree on release <- 15 dropped
 56 touch a nominated hotspot        <- 57 dropped: the 2 surviving `conserved` scaffolds,
                                        in full (POST-HOC gate, declared above)
 20 submitted (cap), 10 from each of 2 scaffolds
```

Slots are filled round-robin by scaffold rather than greedily by switch, because a scaffold the
novelty filter rejects scores nothing whatever its switch is. Filling greedily instead was
measured and changes nothing here: the same median switch, 0.81, and a 12/8 scaffold split
rather than 10/10, because the diversity rule below is what limits how deep either fill reaches.

No two submitted designs differ by fewer than two substitutions. Within a scaffold, mutual
identity is 0.91–0.96; between the two scaffolds it is 0.041–0.061.

Two scaffolds carrying an unpaired cysteine were dropped by default. That selection, and the
pre-registered selection that predates the engagement gate and hedges across both epitopes, are
both recorded alongside this entry.

## The strongest switch in the campaign is not submitted

Re-engineering the published **3DS18** TfR1 binder produced the largest measured switch of the
whole campaign, 1.44–1.57 log units, with the mechanism cleanly established on three independent
crystal copies. It is **excluded as not novel**: its sequence is **100.0% identical** to its own
crystal chain `6wrv_C` (e-value 4e-45, full coverage), and the best histidine variant of it is
still 96.7%. It is kept as the positive control that establishes the burial rule and calibrates
the co-fold score, not as an entry.

## Novelty

Sequence: MMseqs2 against SwissProt and PDB seqres. **All 20 designs have no significant hit**
at e-value 1e-3.

Structure: Foldseek against CATH50, scoring the **design chain only** so the target does not
inflate the match. Both submitted scaffolds fall in the **moderate** band (TM 0.64 and 0.75,
above 95% coverage) and score a predicted 3 of 4. Two scaffolds that scored `high` (TM 0.94 over
100% coverage, and 0.90 over 95%) were excluded, as was the 3DS18 scaffold above.

These scores are a **conservative local proxy** for the organisers' own scorer, not its output.
The fold half of it measures more similarity than theirs appears to, so a `high` verdict reliably
predicts rejection while a `moderate` one may still pass.

## Not established

No binding, expression, or pH response is measured. The switch numbers are predictions from two
pKa calculators, not a measurement of Kd at either pH, and no design reaches the 2 log units that
would mean no detectable binding at pH 6.0.

Touching a nominated hotspot is not binding it, and the engagement measurement is itself
model-dependent: it separates the scaffolds cleanly on the generator's own poses and only
orders them on the independent co-fold, where the submitted designs still sit 3.1–10.8 Å off the
nominated patch. It says nothing about affinity, and the submitted interfaces score below the
level a real TfR1 binder reaches on the same pipeline.

The conserved patch at 529 has **no precedent as an epitope** — no structure shows anything bound
there, so unlike the apical site there is no evidence it is druggable. The 605 patch carried by
this entry has none either. The cut targets expose faces that are interior in the intact
ectodomain. Binding is measured on the homodimer, so avidity may flatter a weak monovalent
design. Mouse cross-reactivity is predicted from a sequence alignment, not from a mouse
structure.

The switch on the *de novo* designs is measured on **one pose per design**, so the replication
axis is the two pKa methods rather than independent structures. The hotspot positions used by the
engagement gate are the five nominated per epitope during target preparation; they were fixed
before any design was generated, but the gate that uses them was not.
