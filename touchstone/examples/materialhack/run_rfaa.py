"""Rank real designed nickel structures with the geometry verifier.

The adapter is source-agnostic: point it at a directory of RFdiffusionAA backbones
(materialhack/run_nickel.sh on pi-a100-80gb) or of LigandMPNN-packed full-atom designs
(run_ligmpnn.sh) — it parses the metal site either way.

    python examples/materialhack/run_rfaa.py <design_pdb_dir>

Judged against the real PDB reference (PDBReference). Bare RFdiffusionAA backbones DEFER —
no sidechains, so no real coordination. After LigandMPNN packs His/Asp/Cys sidechains the
designs gain real donors and the scores climb; full trust still needs clean geometry
(relaxation). That progression is the point — the verifier tracks how real the site is.
"""

import sys

from touchstone import GeometryVerifier, PDBReference, RFdiffusionAdapter, rank


def main(output_dir: str) -> None:
    verifier = GeometryVerifier(PDBReference())
    designs = RFdiffusionAdapter(output_dir, pdb_element="NI", metal_label="Ni2+").design("Ni2+")
    print(f"# RFdiffusionAA → touchstone — {len(designs)} nickel designs\n")
    print(f"{'design':12} {'CN':>3} {'score':>9}  verdict")
    for d, v in rank(designs, verifier):
        flag = "DEFER" if v.ood else ("trust" if v.trust else "weak")
        print(f"{d.sequence:12} {d.site.coordination_number:>3} {v.score:9.2e}  {flag:6} {v.reason}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "output/nickel")
