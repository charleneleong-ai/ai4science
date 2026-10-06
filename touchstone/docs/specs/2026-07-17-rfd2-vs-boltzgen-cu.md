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

**Caveat on the control arm:** this pool is *measured* but not *reproducible from the repo* — there
is no committed Cu run config. `boltzgen_cu_out` appears in the repo only as the path above, and the
only committed BoltzGen reproduce block is the **Ni** motif at `--num_designs 12`
([`boltzgen-metal-design.md`](../boltzgen-metal-design.md)). So pool size, `--protocol`, `--steps`,
and whether the checkpoint was base or RLVR-tuned are all unrecorded for the Cu arm. Commit that
config alongside this spec, or any RFD2-vs-BoltzGen gap is unattributable to the generator.

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

**Score this site at `--cutoff 2.9`, not the 2.8 default.** The axial Cu–S(Met) bond is 2.82 Å, so
the default parse cutoff drops Met92 and the verifier sees **N₂S₁ / CN 3**, not the N₂S₂ / CN 4 the
motif was chosen for. Measured on this exact file:

| cutoff | CN | donors |
|---|---|---|
| 2.80 (default) | 3 | N, S, N |
| **2.85 – 3.00** | **4** | **N, S, N, S** |

It is excluded by 0.02 Å. Three consequences if left at the default: `motif_selectivity` scores the
N₂S₁ enrichment (×3.15) rather than the N₂S₂ one the occupancy analysis cites (×3.8);
`coord_geometry` fits against ideal **CN 3**; and both A/B arms must use the same cutoff or their
coordination numbers aren't comparable. `--cutoff` is exposed on `verify` and `rank`, and the
default stays 2.8 so existing numbers don't move silently.

At 2.9 Å the tiers do pick up the intended motif — `precedent` reports 17 precedents for **Cu-N₂S₂**
and `motif_selectivity` the N₂S₂ enrichment, against Cu-N₂S₁ at the default.

**Calibration note.** At *both* cutoffs this real plastocyanin site lands at **consensus `weak`**,
held there by `bond_valence` (BVS 2.89 vs formal 2 at 2.9 Å; 2.71 at 2.8 Å — the axial S adds
valence, so widening the cutoff makes that tier marginally worse). Worth remembering when reading
the A/B: a *designed* site scoring `weak` is matching experimentally-determined copper, not failing
— and it suggests the BVS tier may be worth re-checking against type-1 Cu specifically.

## The pipeline

    RFdiffusion2 (theozyme → scaffold)  →  LigandMPNN (sequence, metal-aware, REFINE mode)  →
    Chai-1 (fold)  →  touchstone verify_structure (corrected stack)  →  A/B vs BoltzGen

Same generate→inverse-fold→fold→score *shape* as the Ni RLVR loop, but two things do **not** carry
over unchanged, and both affect how the A/B can be read:

- **The arms share no downstream.** The committed BoltzGen run
  ([`boltzgen-metal-design.md`](../boltzgen-metal-design.md)) uses BoltzGen's own
  `--steps design inverse_folding folding` — its inverse-folder and its refold. This arm uses
  LigandMPNN + Chai-1. So the comparison is **pipeline-level, not generator-level**, unless the
  BoltzGen pool is re-run through LigandMPNN + Chai too.
- **LigandMPNN must run in refine mode.** RFD2 with `contig_as_guidepost=True` places its own
  coordinating sidechains, which puts it in the same category as BoltzGen — and the measured cost
  of a full redesign there was **2/24 vs 8/24** with the coordinators fixed
  ([`boltzgen-metal-design.md`](../boltzgen-metal-design.md)). Pass `--fixed_residues` for the motif
  *as numbered in the RFD2 output*; guideposted residues are renumbered in the scaffold, so they are
  not the theozyme's 37/84/87/92.
- **`rlvr_select` is not generator-agnostic.** [`scripts/rlvr_select.py`](../../scripts/rlvr_select.py)
  makes `--npz-dir` a required option and reads BoltzGen `.npz` confidence per design. The
  generator-agnostic claim holds at the `verify_structure` layer and breaks at the harness layer;
  making `--npz-dir` optional is the honest fix before this arm is scored through it.

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
`rf_diffusion/config/inference/aa.yaml`, `rf_diffusion/config/inference/base.yaml`, `setup.py`), so
this is the real interface. Every flag the script passes is listed — nothing is "verified" by
association:

| field | value | why |
|---|---|---|
| `inference.contig_as_guidepost` | `True` | **theozyme mode** — the `active_site_unindexed_atomic` pattern: rotamer-free atomic motif scaffolding, exactly what a metal site needs |
| `contigmap.contig_atoms` | `"{'A37':'ND1,CE1,CD2','A84':'SG,CB','A87':'ND1,CE1,CD2','A92':'SD,CE,CG'}"` | only the coordinating functional groups are constrained — backbone + rotamer stay free. **The enclosing double quotes are required**: Hydra's grammar takes this as a *string*, and a bare `{'A37':…}` dict literal raises `OverrideParseException` |
| `inference.ligand` | `CU` | **verified supported**: RF2AA's `chemical.py` lists `CU`/`CU1`/`CU2` in `METAL_RES_NAMES` (BioLiP metals-in-PDB) and `Cu` in its element vocabulary |
| `contigmap.contigs` | `['10-30,A37-37,10-30,A84-84,10-30,A87-87,10-30,A92-92,10-30']` | variable linkers → scaffold diversity. Passed as a *list*, unlike `contig_atoms` — the asymmetry is Hydra's, not an oversight |
| `--config-name` | `aa` | resolves to `rf_diffusion/config/inference/aa.yaml`, which inherits `base.yaml` |
| `contigmap.length` | `100-140` | real field (`base.yaml` defaults it to `null`); a range is upstream-attested — `open_source_demo.json` passes `contigmap.length=150-150`. Bounds the 54–154 the contig ranges alone permit |
| `inference.ckpt_path` | `…/model_weights/RFD_173.pt` | **a deliberate override, and the right one.** `aa.yaml` ships `RFD_140.pt`, but upstream's docs are explicit: *"there is one recommended model weight file located at `RFdiffusion2/rf_diffusion/model_weights/RFD_173.pt`. This is the set of weights used in the demo in the README"* (`doc/source/usage/usage.rst`). `benchmark/configs/open_source_demo.yaml` passes the same override, so `aa.yaml`'s default is stale and upstream overrides it too. `setup.py` fetches both |

Wrapped in [`scripts/rfd2_cu_design.sh`](../../scripts/rfd2_cu_design.sh). Note RFD2 runs via
**Apptainer**, not conda — `apptainer exec --nv .../bakerlab_rf_diffusion_aa.sif`.

**Which recipe this follows.** Upstream ships *two* ways to scaffold an atomized motif, and they
differ in more than the checkpoint:

| | `open_source_demo.yaml` (what this script follows) | `enzyme_bench_n41.yaml` |
|---|---|---|
| checkpoint | `RFD_173.pt` (the recommended weight) | `RFD_140.pt` |
| motif placement | explicit `contigs` with fixed linker ranges | `contigmap.intersperse='10-100'` + `reintersperse=True` |
| extra | — | `++transforms.configs.CenterPostTransform.center_type='all'` |

Most of the enzyme benchmark's other overrides (`diffuser.T=100`, `idealize_sidechain_outputs`,
`str_self_cond=False`, `rots.sample_schedule=normed_exp`, `guidepost_xyz_as_design_bb`,
`model_runner=NRBStyleSelfCond`) only restate what `aa.yaml` and `base.yaml` already default, so
they are inherited either way — checked, not assumed.

The demo recipe is the right starting point: it is the pattern the `active_site_unindexed_atomic`
mode is documented against, and it uses the recommended weight. But the enzyme benchmark is the
closer *scientific* analogue to a theozyme, and `center_type='all'` is the one flag it sets that
this script does not and that is not an `aa.yaml` default. Worth an arm once stage 1 passes — not
before, since nothing here has been run yet.

**The one real unknown:** the open repo ships **no metal-ion benchmark** — every example ligand is
organic (`LG1`, `NAD`, `OXM`, `PH2`). Zn metallohydrolases are published from this model so metals
plainly work, but the metal path isn't exercised by any shipped config. Hence the script's stage-1
gate, which **exits non-zero** rather than printing advice: every design must carry a `CU` HETATM,
and for a smoke run the first design is put through `touchstone verify --metal Cu2+ --selectivity
Ni2+,Cu2+,Co2+` — the Cu string being present does not mean it is coordinated.

## Stage 1 — passed, 2026-10-06

Installed on `pi-a100-80gb` and run at n=2. The config works end-to-end: RFD2 loaded our
`ckpt_path` override (`RFD_173.pt`, confirmed in its own log), Hydra accepted the `contig_atoms`
string and `contigs` list, diffusion completed, and **the Cu survives into both outputs**.

Two bugs had to be fixed first, neither findable by reading configs:

| | |
|---|---|
| `ModuleNotFoundError: rf_diffusion` | `run_inference.py` imports its own package, so the repo root must be on `PYTHONPATH`. Upstream's `exec/rf_diffusion_aa_shebang.sh` does this, but is unusable off-IPD — it prefers a sif name `setup.py` doesn't install, passes a `--slurm` flag apptainer 1.5.2 rejects, and otherwise falls back to the host python. Fixed with `--env PYTHONPATH`. |
| `ORI HETATM token is required` | The theozyme was missing a **mandatory** input. Added at the Cu coordinates, so the scaffold's centre of mass sits on the metal. |

### The two designs

| | design 0 | design 1 | plastocyanin (reference) |
|---|---|---|---|
| donors | `N, S, S` — **N₁S₂** | `N, N, S` — **N₂S₁** | `N, S, N, S` — N₂S₂ |
| CN | 3 | 3 | 4 |
| geometry | weak (2.0σ) | trust (2.0σ) | trust (0.8σ) |
| bond_valence | **trust (Δ0.15)** | weak (Δ0.64) | weak (Δ0.89) |
| coord_symmetry | weak (0.37) | weak (0.51) | **trust (0.16)** |
| coord_geometry | weak (20.2° vs CN3) | weak (24.6° vs CN3) | trust (14.7° vs CN4) |
| precedent | trust (11× Cu-N1S2) | trust (17× Cu-N2S1) | trust (17× Cu-N2S2) |
| consensus | WEAK | WEAK | WEAK |

**Neither reproduced N₂S₂.** Both came out CN 3, losing a different donor each time — design 0 a
histidine, design 1 the Met thioether. RFD2 scaffolds around the site but drops one donor. That is
the first real signal from this arm, and it is a design-quality result rather than a config fault.

**The sites are under-enclosed**: `coord_symmetry` 0.37 / 0.51 against plastocyanin's 0.16, i.e. the
metal sits one-sided rather than buried. The ORI position is the obvious knob — it was placed on the
Cu precisely to get an enclosed site, and these say it isn't enclosed yet. Now a sweep with a
measurable target.

**Design 0 beat the real protein on bond valence** (Δ0.15, `trust`, vs Δ0.89 `weak`). Together with
the plastocyanin calibration above, that points at the BVS tier being mis-centred for type-1 Cu
rather than the designs failing it.

**Cost:** 939 s for 2 designs (469 s/design) while the GPU was at 100% utilisation from unrelated
jobs. Upper-bound extrapolation: n=24 ≈ 3.1 h, n=96 ≈ 12.5 h.

**What n=2 cannot say:** nothing about trust rate, nothing about the A/B. Two designs from a
contended GPU establish that the pipeline runs and the metal path is real. That is all stage 1 was
for.

### Next
1. ~~Stage-1 smoke test~~ — done.
2. Sweep the ORI placement against `coord_symmetry`, and try the enzyme-benchmark recipe
   (`center_type='all'` + `intersperse`) as a second arm. Note `center_type='all'` *requires* an ORI
   token, which is now present.
3. Scale to a matched n, then LigandMPNN in **refine** mode → Chai. RFD2's `setup.py` ships
   `mlfold.sif` and `chai.sif`, so the downstream can stay containerised.
4. **Blocker:** the BoltzGen Cu control arm has to be rebuilt before any A/B number means anything.

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
