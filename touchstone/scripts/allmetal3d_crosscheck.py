"""Predict a design's metal site with AllMetal3D — a third, specialized co-fold
predictor for the CofoldCrossCheck verifier.

AllMetal3D scans the apo backbone with a 3D-CNN (a different paradigm from
co-folding), so it is independent of both Boltz-2 and Chai-1. Runs on a GPU box
with the allmetal3d package + weights.

Emits a holo PDB (apo protein + the top confident predicted metal as Ni) that
touchstone's cif_provider / coordination_site parses identically to any other
structure — so AllMetal3D drops in as a provider with no verifier change. If no
predicted site clears the probability threshold, emits nothing: the cross-check
then defers (no corroboration available), which is itself an honest signal.

    conda run -n metal3d python scripts/allmetal3d_crosscheck.py design.pdb holo.pdb --pthreshold 0.25
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import typer


def _apo(pdb: Path, out: Path) -> None:
    """Strip metals/waters → apo backbone (AllMetal3D predicts where metals go)."""
    import gemmi

    st = gemmi.read_structure(str(pdb))
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()
    st.write_pdb(str(out))


def _top_metal(metals_pdb: Path) -> tuple[float, np.ndarray | None]:
    """Highest-probability predicted metal (AllMetal3D writes prob in occupancy)."""
    xs = []
    for line in metals_pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            xs.append((float(line[54:60]), xyz))
    return max(xs, key=lambda t: t[0]) if xs else (0.0, None)


def main(pdb: Path, out: Path, pthreshold: float = 0.25) -> None:
    import gemmi

    work = Path(tempfile.mkdtemp())
    apo = work / "apo.pdb"
    _apo(pdb, apo)
    # run AllMetal3D at a low floor so it emits the density peaks; gate on pthreshold below
    subprocess.run(
        [str(Path(sys.executable).with_name("allmetal3d")), "-i", str(apo),
         "-a", "allmetal3d", "-m", "fast", "-p", str(min(pthreshold, 0.1)), "-o", str(work)],
        check=True,
    )
    metals = next(work.glob("*_metals.pdb"), None)
    prob, xyz = _top_metal(metals) if metals else (0.0, None)
    if xyz is None or prob < pthreshold:
        print(f"no confident metal (top prob {prob:.2f} < {pthreshold}) — no prediction emitted")
        return

    # inject the predicted metal as Ni via gemmi (correct PDB columns ⇒ parses downstream)
    merged = gemmi.read_structure(str(apo))
    chain = merged[0].add_chain("Z")
    res = gemmi.Residue()
    res.name, res.seqid = "NI", gemmi.SeqId("1")
    atom = gemmi.Atom()
    atom.name, atom.element, atom.pos = "NI", gemmi.Element("Ni"), gemmi.Position(*map(float, xyz))
    res.add_atom(atom)
    chain.add_residue(res)
    merged.write_pdb(str(out))
    print(f"wrote {out}: apo + predicted Ni (prob {prob:.2f})")


if __name__ == "__main__":
    typer.run(main)
