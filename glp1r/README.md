# GLP1R binders — Anthropic × Adaptyv Challenge 3

A de novo protein that **binds human GLP1R and switches it on**. Affinity is ranked first and
activation second, but activation is a gate: the organisers discard a binder that does not
activate. Submissions closed **Monday 26 October 2026, 22:59 Sydney**.

**Status (2026-10-24): 20 designs submitted — extracellular-domain binders, with agonism
explicitly not claimed.** The entry is in
[`submission/glp1r_submission_v2.csv`](submission/glp1r_submission_v2.csv) with the per-design
evidence beside it and the submitted text in
[`submission/methodology.md`](submission/methodology.md). Nothing is wet-lab measured.

## Why the entry does not claim agonism

Measured off PDB 6X18 rather than taken from a review: GLP1R's two ligand sites are **disjoint**.
The 26 receptor residues contacted by GLP-1's activation segment (residues 7–13) are all in the
transmembrane bundle and loops; the 19 ECD residues forming the affinity groove are contacted
only by the peptide's C-terminal half. Just three residues (141, 201, 298) see both halves. The
centroids are **31.6 Å** apart and the peptide N-terminus sits **19.2 Å** below the vestibule rim.

Six independent campaigns have put a protein in that groove — including a 10¹⁰ antibody library
carrying >150,000 GPCR motifs in HCDR3, and a 10,080-design de novo campaign read out in a
CRE-luciferase line. **Zero full agonists**; the best ECD-restricted agonism reported is Emax
<2% of exendin-4. The groove is the affinity trap, not the switch, so these 20 are submitted as
what the evidence supports: binders.

The activation gate was built anyway, and calibrating it is the most useful negative result here
— on *predicted* structures it **rejects but does not select**. See
[`docs/calibration.md`](docs/calibration.md) §5: an inactive 15-mer ranked second, above two full
agonists, and scoring one sample gave native GLP-1 0/6 against its experimental 6/6.

## The funnel

```
200 designs   4 hotspot variants x 50, BoltzGen against the ECD core (29-128)
140 triage    70% — 46 rejected for alanine over-representation, 27 for Ala+Gly
 23 iPTM      also reach BoltzGen iPTM >= 0.5
 21 Chai-1    confirmed by an independent model at >= 0.65
 20 submitted cap
```

| | BoltzGen (generator) | Chai-1 (independent) |
| --- | --- | --- |
| median | 0.535 | **0.802** |
| range | 0.502–0.616 | 0.150–0.862 |

Spearman between them is **+0.377**, so **no fine-grained ranking is claimed** — row order is by
Chai-1, the model that did not generate them, and is labelled weak. The disagreement did real
work: two designs at 0.551 and 0.516 on the generator fell to 0.236 and 0.150 on Chai-1 and were
dropped, and both would have made a top-20 ranked on the generator alone.

Mutual sequence identity among the 20 is **0.011–0.317, median 0.096** — a diverse set, not a
substitution series. (The TNF-α entry's 20 designs sat at 0.934–0.989.)

**Novelty.** All 20 clear the 30% sequence-identity condition on any significant measure: the
best hit for any design has e-value 0.27 against 1.75M sequences, most exceed 10. Structurally,
19 of 20 fall in Foldseek's "moderate" band; `glp1r_20` is **HIGH at TM 0.844** and may be
rejected by the de novo filter, so it is ranked last deliberately. Its best alternative scored
worse (TM 0.869, e-value 0.03), so no swap was available.

## Reproducing

Inputs are public and `data/` is gitignored, so regenerate rather than commit it:

```bash
mkdir -p glp1r/data
curl -sL -o glp1r/data/6x18.cif https://files.rcsb.org/download/6X18.cif
curl -s "https://rest.uniprot.org/uniprotkb/P43220.json?fields=ft_topo_dom,ft_transmem,ft_signal,sequence" \
  -o glp1r/data/p43220.json

uv run --with gemmi --with typer python3 glp1r/scripts/target_sites.py \
  --structure glp1r/data/6x18.cif --out glp1r/data/6x18_contacts.csv

uv run --with gemmi --with typer python3 glp1r/scripts/prepare_target.py \
  --structure glp1r/data/6x18.cif --contacts glp1r/data/6x18_contacts.csv \
  --uniprot glp1r/data/p43220.json --out-dir glp1r/data/target
```

| script | what it does |
| --- | --- |
| [`target_sites.py`](scripts/target_sites.py) | assigns every receptor residue within 4.5 Å of GLP-1 to the groove or the activation segment against the UniProt P43220 topology, and reports span and insertion depth |
| [`prepare_target.py`](scripts/prepare_target.py) | cuts `receptor.pdb`, `ecd.pdb` (24–139) and `ecd_design.pdb` (29–128, the core ordered across all 13 structures resolving it), plus `hotspots.json`; fails loudly if structure and UniProt disagree on numbering |
| [`activation_geometry.py`](scripts/activation_geometry.py) | the six-feature pharmacophore gate plus insertion depth |
| [`calibrate_scores.py`](scripts/calibrate_scores.py) | AUROC of the public predicted scores against wet-lab outcomes; names constant placeholder fields instead of scoring them |
| [`triage_designs.py`](scripts/triage_designs.py) | composition degeneracy, free cysteine and the Glu-fraction expression predictor |
| [`chai_independent.py`](scripts/chai_independent.py) | wraps `run_folding_on_context` to capture the `ranking_data` that `run_inference` discards |
| [`compare_models.py`](scripts/compare_models.py) | Spearman between generator and independent model, binder-to-target only |

## What is not established

No measured binding, expression or activation. Two predictors agreeing raises confidence but is
not a measurement — on the public corpus, ipTM's own within-target AUROC for binding is 0.607 and
nothing tested reaches 0.66. **Agonism is explicitly not claimed.** The novelty thresholds have
not been verified against the organisers' own scorer. The calibration corpus contains **zero
GLP1R rows and no GPCR of any kind**, so every AUROC here is borrowed from soluble targets; the
thermostabilising mutations in the assay construct (Sino Biological 13944-H18H-DD) are
undisclosed, and affinity is scored on purified detergent receptor while activation is scored in
live cells, so those are not the same target.
