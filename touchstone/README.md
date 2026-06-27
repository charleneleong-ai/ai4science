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
