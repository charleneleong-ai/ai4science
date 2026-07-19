#!/usr/bin/env bash
# Generate Cu2+ type-1 (blue copper) sites with RFdiffusion2, for the A/B against BoltzGen.
#
# Config verified against the RFdiffusion2 repo (not guessed):
#   - theozyme mode is `inference.contig_as_guidepost=True` + `contigmap.contig_atoms`, the
#     "active_site_unindexed_atomic" pattern in rf_diffusion/benchmark/open_source_demo.json —
#     rotamer-free atomic motif scaffolding, which is what a metal theozyme needs.
#   - `inference.ligand=CU` is supported: RF2AA's chemical.py lists CU/CU1/CU2 in METAL_RES_NAMES
#     (BioLiP "metals in PDB") and Cu in its element vocabulary.
#
# CAVEAT, stated plainly: the open repo ships NO metal-ion benchmark (only organic ligands — LG1,
# NAD, OXM, PH2). Zn metallohydrolases are published from this model, so metals clearly work, but
# the exact metal path is not exercised by any shipped config. Run STAGE 1 (smoke test, 2 designs)
# and confirm the Cu survives into the output before spending a full batch.
#
# Usage (on the GPU box, after cloning RFdiffusion2 and running its setup.py):
#   RFD2_ROOT=~/RFdiffusion2 ./rfd2_cu_design.sh 2      # smoke test
#   RFD2_ROOT=~/RFdiffusion2 ./rfd2_cu_design.sh 96     # full batch, matches the BoltzGen n
set -euo pipefail

N=${1:-2}
RFD2_ROOT=${RFD2_ROOT:?set RFD2_ROOT to the RFdiffusion2 clone}
THEOZYME=${THEOZYME:-$(cd "$(dirname "$0")/.." && pwd)/examples/cu_type1_theozyme.pdb}
OUT=${OUT:-$HOME/materialhack/rfd2_cu_out}
SIF=${SIF:-$RFD2_ROOT/rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif}

# The type-1 Cu motif: His37 + His87 (imidazole N donors), Cys84 (thiolate S), Met92 (thioether S).
# Only the coordinating functional-group atoms are constrained — backbone and rotamer are free, so
# RFD2 builds a scaffold that satisfies the coordination geometry rather than copying plastocyanin.
CONTIG_ATOMS="{'A37':'ND1,CE1,CD2','A84':'SG,CB','A87':'ND1,CE1,CD2','A92':'SD,CE,CG'}"
# Variable-length linkers between the four motif residues → scaffold diversity at fixed total length.
CONTIGS="['10-30,A37-37,10-30,A84-84,10-30,A87-87,10-30,A92-92,10-30']"

mkdir -p "$OUT"
echo "RFD2 Cu2+ type-1 design: n=$N  theozyme=$THEOZYME  out=$OUT"

apptainer exec --nv "$SIF" python "$RFD2_ROOT/rf_diffusion/run_inference.py" \
  --config-name=aa \
  inference.input_pdb="$THEOZYME" \
  inference.ligand=CU \
  inference.contig_as_guidepost=True \
  contigmap.contigs="$CONTIGS" \
  contigmap.contig_atoms="\"$CONTIG_ATOMS\"" \
  contigmap.length=100-140 \
  inference.num_designs="$N" \
  inference.output_prefix="$OUT/cu_t1" \
  inference.ckpt_path="$RFD2_ROOT/rf_diffusion/model_weights/RFD_173.pt"

echo
echo "wrote $(ls "$OUT"/*.pdb 2>/dev/null | wc -l) designs to $OUT"
echo "CHECK BEFORE SCALING UP: the Cu must survive into the output —"
echo "  grep -h 'CU' $OUT/cu_t1_0.pdb | head"
echo "Then: LigandMPNN (ligmpnn env) -> Chai fold (chai env) -> touchstone verify --metal Cu2+"
