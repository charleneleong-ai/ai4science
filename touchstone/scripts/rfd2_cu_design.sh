#!/usr/bin/env bash
# Generate Cu2+ type-1 (blue copper) sites with RFdiffusion2, for the A/B against BoltzGen.
#
# Config provenance, the per-flag verification table, and the no-metal-benchmark caveat live in
# the spec — it is also where the A/B result gets filled in, so it stays the single source:
#   touchstone/docs/specs/2026-07-17-rfd2-vs-boltzgen-cu.md
#
# Usage (GPU box, after cloning RFdiffusion2 + running its setup.py; RFD2 runs under Apptainer):
#   RFD2_ROOT=~/RFdiffusion2 ./rfd2_cu_design.sh 2     # stage-1 smoke test (foreground)
#
#   # full batch — hours, so run it detached per the long-job convention. `env` is required:
#   # nohup will not accept a leading VAR=value assignment.
#   mkdir -p logs && setsid nohup env RFD2_ROOT=~/RFdiffusion2 ./rfd2_cu_design.sh 96 \
#     </dev/null >>logs/rfd2_cu_$(date -u +%Y%m%dT%H%M%SZ).log 2>&1 & disown
#   ps -o ppid= -p $!   # must print 1
#
# Overridable: RFD2_ROOT (required), THEOZYME, OUT, SIF, CKPT, RUN_ID.
set -euo pipefail

N=${1:-2}
(( N > 0 )) || { echo "n must be > 0" >&2; exit 1; }

RFD2_ROOT=${RFD2_ROOT:?set RFD2_ROOT to the RFdiffusion2 clone}
# Absolute: a relative path handed into the container resolves against the container's cwd.
TOUCHSTONE=$(cd "$(dirname "$0")/.." && pwd)
THEOZYME=${THEOZYME:-$TOUCHSTONE/examples/cu_type1_theozyme.pdb}
SIF=${SIF:-$RFD2_ROOT/rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif}
# RFD_173 is upstream's one recommended weight (doc/source/usage/usage.rst) and what the README
# demo uses; aa.yaml's RFD_140 default is stale, and upstream's own demo config overrides it too.
CKPT=${CKPT:-$RFD2_ROOT/rf_diffusion/model_weights/RFD_173.pt}
# Per-run dir: stage 2 must not overwrite the stage-1 designs that justified scaling up.
RUN_ID=${RUN_ID:-n${N}_$(date -u +%Y%m%dT%H%M%SZ)}
OUT=${OUT:-$HOME/materialhack/rfd2_cu_out/$RUN_ID}

command -v apptainer >/dev/null || { echo "apptainer not on PATH" >&2; exit 1; }
nvidia-smi -L >/dev/null 2>&1 || { echo "no visible GPU" >&2; exit 1; }
for f in "$SIF" "$THEOZYME" "$CKPT"; do
  [[ -s $f ]] || { echo "missing or empty: $f" >&2; exit 1; }
done

# The type-1 Cu motif: His37 + His87 (imidazole N donors), Cys84 (thiolate S), Met92 (thioether S).
# Only the coordinating functional-group atoms are constrained — backbone and rotamer stay free, so
# RFD2 satisfies the coordination geometry rather than copying plastocyanin.
#
# Hydra override grammar, verified against its parser: `contigs` must arrive as a LIST, while
# `contig_atoms` must arrive as a quoted STRING — a bare {'A37':...} dict literal raises
# OverrideParseException. The inner quotes below are load-bearing; do not unify the two forms.
CONTIGS="['10-30,A37-37,10-30,A84-84,10-30,A87-87,10-30,A92-92,10-30']"
CONTIG_ATOMS="\"{'A37':'ND1,CE1,CD2','A84':'SG,CB','A87':'ND1,CE1,CD2','A92':'SD,CE,CG'}\""

mkdir -p "$OUT"
echo "[$(date -u +%FT%TZ)] RFD2 Cu2+ type-1: n=$N out=$OUT ckpt=$(basename "$CKPT")"

# -u: unbuffered, so a tailing monitor sees progress and a crash doesn't swallow the last buffer.
apptainer exec --nv "$SIF" python -u "$RFD2_ROOT/rf_diffusion/run_inference.py" \
  --config-name=aa \
  inference.input_pdb="$THEOZYME" \
  inference.ligand=CU \
  inference.contig_as_guidepost=True \
  contigmap.contigs="$CONTIGS" \
  contigmap.contig_atoms="$CONTIG_ATOMS" \
  contigmap.length=100-140 \
  inference.num_designs="$N" \
  inference.output_prefix="$OUT/cu_t1" \
  inference.ckpt_path="$CKPT"

echo "[$(date -u +%FT%TZ)] inference done in ${SECONDS}s ($((SECONDS / N))s/design)"

# --- stage-1 gate: enforced, not advised -------------------------------------------------------
# No shipped RFD2 config exercises the metal path, so a run can finish clean and carry no Cu at
# all. A 96-batch that silently dropped the Cu would score as "RFD2 is worse than BoltzGen at Cu"
# — a config artefact published as the A/B's headline. Assert it here instead.
shopt -s nullglob
designs=("$OUT"/cu_t1_*.pdb)
(( ${#designs[@]} >= N )) || { echo "FAIL: ${#designs[@]} designs written, expected $N" >&2; exit 1; }
for d in "${designs[@]}"; do
  grep -qE '^HETATM.* CU ' "$d" || { echo "FAIL: no Cu HETATM in $d — do NOT scale up" >&2; exit 1; }
done
echo "ok: ${#designs[@]} designs in $OUT, Cu retained in all"

# A Cu HETATM can still sit well off the motif, so the smoke run also faces the real verifier —
# the same call that produced the committed BoltzGen baseline. --selectivity is what instantiates
# the motif-selectivity tier; without it the two arms are scored over different tier sets.
# --cutoff 2.9: the Met92 thioether sits at 2.82 A, so the 2.8 default drops it and the site scores
# as N2S1/CN3 instead of the N2S2/CN4 the motif was chosen for. Measured on the input theozyme.
if (( N <= 2 )); then
  command -v uv >/dev/null || {
    echo "BLOCKED: uv not found, so the geometry gate did not run. Before scaling up:" >&2
    echo "  touchstone verify ${designs[0]} --metal Cu2+ --cutoff 2.9 --selectivity Ni2+,Cu2+,Co2+" >&2
    exit 1
  }
  uv run --directory "$TOUCHSTONE" touchstone verify "${designs[0]}" \
    --metal Cu2+ --cutoff 2.9 --selectivity Ni2+,Cu2+,Co2+
fi

echo
echo "Next: LigandMPNN in REFINE mode — a full redesign discards the motif"
echo "  (boltzgen-metal-design.md: 2/24 redesigned vs 8/24 with the coordinators fixed)."
echo "  Pass --fixed_residues for the motif as numbered IN THE RFD2 OUTPUT — guideposted"
echo "  residues are renumbered in the scaffold, so they are not the theozyme's 37/84/87/92."
echo "Then: Chai fold -> touchstone verify --metal Cu2+ --cutoff 2.9 --selectivity Ni2+,Cu2+,Co2+"
echo "  (RFD2's setup.py ships chai.sif and mlfold.sif, so the downstream can stay containerised.)"
