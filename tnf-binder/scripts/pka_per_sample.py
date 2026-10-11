"""Which diffusion sample gives which pKa — the mapping the robustness scan collapsed away.

`pka_robustness.py` reports min, max and mean per position across samples, which is what the
ranking needs and which throws away the association between a pKa and the structure that
produced it. To ask whether the energetically better structures give the more credible pKa,
that association is exactly what is required.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, "/home/ubuntu/tnf/scripts")
from Bio.PDB import PDBIO  # noqa: E402
from pka_switch import run_propka  # noqa: E402
from score_designs import classify_chains, load  # noqa: E402


def main() -> None:
    directory = Path(sys.argv[1])
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        for cif in sorted(directory.glob("*_model_*.cif")):
            model = load(cif)
            _target, binders = classify_chains(model)
            chain = binders[0]
            positions = [r.id[1] for r in chain if r.get_resname() == "HIS"]
            pdb = work / f"{cif.stem}.pdb"
            io = PDBIO()
            io.set_structure(model)
            io.save(str(pdb))
            pkas = run_propka(pdb, work)
            shown = "  ".join(
                f"His{p}={pkas.get((chain.id, p)):.2f}" if pkas.get((chain.id, p)) is not None
                else f"His{p}=NA" for p in positions)
            print(f"  {cif.stem[-7:]}  chain {chain.id}  {shown}", flush=True)


if __name__ == "__main__":
    main()
