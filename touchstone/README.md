# touchstone

A generator-agnostic **verifier** for designed metal-binding proteins. The generator
(BoltzGen, RFdiffusionAA→LigandMPNN, Chai, …) proposes a metal-coordination site;
touchstone judges whether that *predicted* site is real enough to make — a trust/weak/defer
consensus across independent methods (defense-in-depth, ≥2 per stage):

- **geometry** — z-score vs a CSD/PDB prior · bond-valence sum · Mogul CSD validation
- **precedent** — CSD CrossMiner (is the motif seen in nature?)
- **expression** — ESM-2 pseudo-perplexity · solubility
- **physics / dynamics** — xtb GFN2 · MLIP (MACE) relaxation + MD
- **thermostability** — site MLIP-MD · global Tm (TemStaPro)
- **selectivity** — geometry profile · MLIP ΔE metal-swap
- **cross-verification** — independent co-fold (Boltz-2 / Chai-1 / AllMetal3D): does the
  generator's predicted site reproduce?

Use it to **triage designs to wet-lab** (only `trust` clears the bar) and to **score them
as an RLVR reward** to iterate the generator.

## Use in Claude Code

Install the plugin (bundles the MCP server + the `verify-metal-binder` skill + a
`/verify-binder` command):

```
/plugin marketplace add charleneleong-ai/ai4science
/plugin install touchstone@ai4science
```

Then `/verify-binder design.pdb Ni2+`, or just ask Claude to verify a design — it calls the
`verify_metal_binder` MCP tool. Requires [`uv`](https://docs.astral.sh/uv/) on PATH (the
plugin runs the server via `uv run`).

To register the MCP server directly instead of via the plugin:

```bash
uv tool install "touchstone[mcp] @ git+https://github.com/charleneleong-ai/ai4science.git#subdirectory=touchstone"
claude mcp add -s user touchstone -- touchstone-mcp
```

## CLI

```bash
touchstone verify design.pdb --metal Ni2+     # instant: geometry + bond-valence + CSD
touchstone rank designs/*.pdb --metal Ni2+    # batch, best-first by reward
touchstone verify design.pdb --deep           # + MLIP relax/MD (needs a GPU)
```

Compare a generator's own confidence with the verifier's verdict per design (BoltzGen iPTM/pLDDT/pTM vs touchstone consensus + the reason it deferred), optionally logging to W&B:

```bash
uv run --extra viz python scripts/boltzgen_scores.py \
  --npz-dir <fold_out_npz> --cif-dir <refold_cif> --metal Ni2+ --wandb
```

## Sample output

Each result is a JSON-able dict: per-tier verdicts (`label` / `score` / `reason` + a
machine-readable `metrics` block), a `stack` listing every tier with its `status`
(`ran` / `skipped` / `needs_input`), and the trust/weak/defer `consensus`. A real
LigandMPNN Ni pack, verified **without a GPU** (the default, runs anywhere):

```jsonc
{
  "metal": "Ni2+", "coordination_number": 5, "donors": ["O","O","N","O","N"], "reference": "PDB",
  "verifiers": {
    "geometry":     { "label": "weak",  "score": 0.025, "reason": "strained geometry (2.3σ)",
                      "metrics": { "strain_sigma": 2.32, "cn": 5, "cn_modal": 4 } },
    "bond_valence": { "label": "defer", "score": 0.023, "reason": "BVS 0.90 vs formal 2 (Δ1.10) — defer",
                      "metrics": { "bvs": 0.9, "formal_valence": 2, "delta": 1.1 } }
  },
  "stack": [
    { "stage": "geometry",     "status": "ran" },
    { "stage": "bond_valence", "status": "ran" },
    { "stage": "mogul",        "status": "needs_input", "detail": "a CSD licence (Mogul / CSD Python API)" },
    { "stage": "mlip",         "status": "needs_input", "detail": "pass deep=True (needs a GPU backend)" },
    { "stage": "mlip_md",      "status": "needs_input", "detail": "pass deep=True (needs a GPU backend)" }
    // … cofold, expression, thermostability: needs_input
  ],
  "consensus": "defer"
}
```

With **`--deep`** on a GPU, the two MLIP tiers flip from `needs_input` to `ran` and add
quantitative physics (the rest of the result is identical):

```jsonc
"mlip":    { "label": "defer", "score": 0.088, "reason": "site lost 2 donor(s), drift 1.92 Å, ΔE_bind -3.33 eV — defer",
             "metrics": { "drift_angstrom": 1.92, "cn_before": 5, "cn_after": 3, "interaction_energy_ev": -3.327 } },
"mlip_md": { "label": "defer", "score": 0.059, "reason": "shell survived 6% of 300 K MD — defer",
             "metrics": { "retention": 0.059, "cn_initial": 5, "temperature_k": 300.0 } }
```

Here the GPU physics **confirms** the instant tiers' `defer`: under MACE relaxation the site
drifts ~2 Å and loses 2 of 5 donors, and only 6 % of the first shell survives 300 K MD.
Add `--stress` for a robustness map (`neutral` / `leachate` / `low_pH`) on top of either.

## Remote / deep tiers

The default tiers run anywhere. The deep tiers (MLIP/xtb/ESM) need a GPU, and CSD/Mogul
need a licence — host the server where those live and point clients at it over HTTP:

```bash
touchstone-mcp --http --host 0.0.0.0 --port 8000          # on the GPU box
claude mcp add --transport http touchstone http://<host>:8000/mcp   # on each client
```

## Scope

The trust threshold is grounded in CSD geometry + physics, **not yet calibrated to wet-lab
outcomes** — `trust` means "physically / precedent-plausible," not a calibrated binding
probability. Tiers without their model/licence available report as `not_run` rather than
guessing.

## Develop

```bash
uv sync --extra dev
uv run --extra dev pytest
```
