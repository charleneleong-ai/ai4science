# RFdiffusion2 vs BoltzGen for Cu²⁺ site design — an A/B through the corrected stack (2026-07-17)

touchstone is generator-agnostic ("the generator is a commodity; the verifier is the asset"), so a
better generator just plugs in and is scored identically. This spec sets up the A/B: **does
RFdiffusion2 produce better Cu²⁺ coordination sites than BoltzGen, measured through the corrected
touchstone stack?**

## Why RFdiffusion2

BoltzGen is a fold-centric binder generator — it infers the metal geometry through folding rather
than building it. RFdiffusion2 (Baker lab, *Nature Methods* 2025; open, single-A100) does **atom-level
theozyme scaffolding**: given the metal + coordinating functional-group atoms, it co-diffuses the
backbone and catalytic residues to satisfy the *exact* coordination geometry. It is the generator
behind the Dec-2025 [computational metallohydrolases](https://pmc.ncbi.nlm.nih.gov/articles/PMC12727532/)
(wet-lab-validated Zn sites). For a designed Cu²⁺ site — where coordination geometry *is* the
objective — conditioning on the metal during backbone generation is the right inductive bias.

Runners-up: **RFdiffusion3** (Dec 2025, 10× faster, but metal-motif input undocumented — confirm
before use); **RFdiffusionAA** (our existing path, superseded by RFD2 for precise geometry). No
model has a *published* Cu²⁺ example — Zn is the demonstrated transition metal — so Cu²⁺ is
extrapolation from a metal-agnostic mechanism. See the model-landscape research
([research agent, 2026-07-17]).

## The baseline is already measured

The 96 existing BoltzGen Cu²⁺ designs (`~/materialhack/boltzgen_cu_out`), scored through the
**corrected** stack (`verify_structure(cif, "Cu2+", selectivity_metals=("Ni2+","Cu2+","Co2+"))`):

| | trust | weak | defer | mean reward |
|---|---|---|---|---|
| **BoltzGen Cu²⁺ (n=96)** | **14** | 25 | 57 | **0.197** |

Per-tier: geometry 62% trust · bond_valence 42%/33% defer · coord_symmetry 62% · coord_geometry
43%/27% defer · precedent 68% · **motif_selectivity 64% trust**. Two reads: (1) the Cu spec *did*
produce Cu-characteristic sites (64% motif_selectivity trust — the donor sets are more Cu- than
Ni/Co-characteristic); (2) the bottleneck is coordination quality — bond_valence and coord_geometry
defer ~30%, i.e. the metal is imperfectly coordinated. **That is exactly what an all-atom,
metal-conditioned generator should improve** — the A/B hypothesis.

## The theozyme

[`examples/cu_type1_theozyme.pdb`](../../examples/cu_type1_theozyme.pdb) — the type-1 "blue copper"
site from plastocyanin **1PLC (1.33 Å)**: His37/His87 (N) + Cys84 (S) + Met92 (S), the **N₂S₂**
donor set. Chosen because the [occupancy analysis](../experiments/2026-07-14-selectivity-from-occupancy.md)
shows N₂S₂ is the most Cu-characteristic donor set (5.5%, ×3.8 over Ni/Co). Geometry verified
physical: Cu–N 1.91/2.06 Å, Cu–S(Cys) 2.07 Å, Cu–S(Met) 2.82 Å (long axial). *(A low-resolution
azurin structure was rejected first — its Cu–S(Cys) came out at 1.79 Å, physically impossible.)*

## The pipeline (mirrors the BoltzGen path)

    RFdiffusion2 (theozyme → scaffold)  →  LigandMPNN (sequence, metal-aware)  →  Chai-1 (fold)  →
    touchstone verify_structure (corrected stack)  →  A/B vs the BoltzGen baseline

Same generate→inverse-fold→fold→score shape as the Ni RLVR loop, so `rlvr_select` and the scoring
harness carry over unchanged.

## To run — RFD2 install is yours (external `git clone`)

I can't clone external repos. Run these on `pi-a100-80gb` (the `!`-prefix works, or a detached
daemon per the long-job convention):

```bash
# 1. clone
git clone https://github.com/RosettaCommons/RFdiffusion2.git ~/RFdiffusion2

# 2. weights (~30 min — "these files are quite large", per their README)
cd ~/RFdiffusion2 && export PYTHONPATH=~/RFdiffusion2 && python setup.py

# 3. Apptainer — RFD2 runs in a container, NOT a conda env
sudo add-apt-repository -y ppa:apptainer/ppa && sudo apt update && sudo apt install -y apptainer
```

## The run config — verified against the repo, not guessed

Read from RFdiffusion2's own configs (`rf_diffusion/benchmark/open_source_demo.json`,
`configs/inference/aa.yml`), so this is the real interface:

| field | value | why |
|---|---|---|
| `inference.contig_as_guidepost` | `True` | **theozyme mode** — the `active_site_unindexed_atomic` pattern: rotamer-free atomic motif scaffolding, exactly what a metal site needs |
| `contigmap.contig_atoms` | `{'A37':'ND1,CE1,CD2','A84':'SG,CB','A87':'ND1,CE1,CD2','A92':'SD,CE,CG'}` | only the coordinating functional groups are constrained — backbone + rotamer stay free |
| `inference.ligand` | `CU` | **verified supported**: RF2AA's `chemical.py` lists `CU`/`CU1`/`CU2` in `METAL_RES_NAMES` (BioLiP metals-in-PDB) and `Cu` in its element vocabulary |
| `contigmap.contigs` | `['10-30,A37-37,10-30,A84-84,10-30,A87-87,10-30,A92-92,10-30']` | variable linkers → scaffold diversity at fixed total length |

Wrapped in [`scripts/rfd2_cu_design.sh`](../../scripts/rfd2_cu_design.sh). Note RFD2 runs via
**Apptainer**, not conda — `apptainer exec --nv .../bakerlab_rf_diffusion_aa.sif`.

**The one real unknown:** the open repo ships **no metal-ion benchmark** — every example ligand is
organic (`LG1`, `NAD`, `OXM`, `PH2`). Zn metallohydrolases are published from this model so metals
plainly work, but the metal path isn't exercised by any shipped config. Hence the script's stage-1
smoke test: **generate 2 designs and confirm the Cu survives into the output** before committing a
full batch.

Then tell me it's installed and I'll:
1. Run the stage-1 smoke test and confirm the Cu is retained + the coordination geometry is sane.
2. Scale to n=96 (matching the BoltzGen pool), then LigandMPNN (`ligmpnn` env) → Chai (`chai` env).
3. Score through the corrected stack and fill in the A/B table below.

## A/B result (pending RFD2 run)

| generator | n | trust | weak | defer | mean reward | motif_selectivity trust | bond_valence defer |
|---|---|---|---|---|---|---|---|
| BoltzGen Cu²⁺ | 96 | 14 | 25 | 57 | 0.197 | 64% | 33% |
| RFdiffusion2 Cu²⁺ | — | — | — | — | — | — | — |

**Hypothesis:** RFD2's metal-conditioned scaffolding lifts the coordination-quality tiers
(bond_valence, coord_geometry) that gate the BoltzGen designs — mean reward and full-consensus trust
rise. If it doesn't, that's also a real result: it says fold-inferred geometry (BoltzGen) is already
competitive for this donor set, and the lever is elsewhere.

## Honest flags

- No published Cu²⁺ RFD2 design exists; Zn is the demonstrated metal. This A/B is partly a test of
  whether the metal-agnostic mechanism transfers to Cu.
- The N₂S₂ type-1 motif buys Cu-over-Ni (×3.8) but **not** Cu-over-Zn (only ×1.12) — see the
  occupancy writeup. A Cu design that must reject Zn needs more than this donor set; that's a
  motif-design question independent of the generator.
- Selectivity as a *thermodynamic* ranking (ΔΔG) is still unavailable (no MLIP passes the
  Irving–Williams gate). The A/B is scored on geometry + coordination + precedent + occupancy, not
  on a binding free energy.
