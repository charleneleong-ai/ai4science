# GLP1R agonist — Anthropic × Adaptyv Challenge 3

A de novo protein that **binds human GLP1R and switches it on**. Affinity is ranked first and
activation second, but activation is a gate: the organisers discard a binder that does not
activate. Submissions close **Sunday 25 October 23:59 AoE** (Monday 26 October 22:59 Sydney).

**Status (2026-10-16): target characterised, nothing designed yet.** The strategy is waiting on
four open questions listed at the end of [`docs/target.md`](docs/target.md).

## The finding that shapes the entry

Measured off PDB 6X18 rather than taken from a review: GLP1R's two ligand sites are **disjoint**.
The 26 receptor residues contacted by GLP-1's activation segment (residues 7–13) are all in the
transmembrane bundle and loops; the 19 ECD residues forming the affinity groove are contacted
only by the peptide's C-terminal half. Just three residues (141, 201, 298) see both halves, all
at the top of the bundle.

The two site centroids are **31.6 Å** apart and the peptide N-terminus sits **19.2 Å** below the
vestibule rim. So a molecule that both binds and activates has to span 31.6 Å and terminate deep
inside the bundle — which an antibody-derived format cannot do, since it recognises a surface
patch and cannot thread a free terminus down a channel. That makes `protein` (single chain with
a protruding N-terminal arm) the only format with a mechanical route to agonism, and it means a
groove-only binder has no path to switching the receptor on.

## Reproducing the target definition

Inputs are public and the `data/` tree is gitignored, so regenerate rather than commit it:

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

- [`scripts/target_sites.py`](scripts/target_sites.py) — assigns every receptor residue within
  4.5 Å of GLP-1 to the groove or the activation segment, against the UniProt P43220 topology,
  and reports the span and insertion depth.
- [`scripts/prepare_target.py`](scripts/prepare_target.py) — cuts `receptor.pdb` (chain R as
  modelled, 29–423), `ecd.pdb` (24–139), the assay-construct and ECD FASTA, and `hotspots.json`.
  It fails loudly if the structure and UniProt disagree on numbering, which would silently
  offset every hotspot.

## What is not established

Deliberately not yet claimed: whether any de novo design has agonised a class B GPCR, whether
any current tool can place a designed terminus in a deep pocket, whether ipTM means anything for
a 7TM target modelled without a membrane, and whether the undisclosed thermostabilising
mutations in the assay construct (Sino Biological 13944-H18H-DD) narrow the orthosteric pocket.
Affinity is scored on that purified detergent construct while activation is scored in live
cells, so those are not the same target.
