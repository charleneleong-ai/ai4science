"""Rank real designed nickel structures with the geometry verifier.

The adapter is source-agnostic: point it at a directory of RFdiffusionAA backbones or
of LigandMPNN-packed full-atom designs — it parses the metal site either way.

    python examples/run_rfaa.py <design_pdb_dir>

Judged against the real PDB reference. Bare RFdiffusionAA backbones DEFER — no sidechains,
so no real coordination. After LigandMPNN packs His/Asp/Cys sidechains the designs gain
real donors and the scores climb; full trust still needs clean geometry. That progression
is the point — the verifier tracks how real the site is.
"""

import typer

from touchstone import GeometryVerifier, PDBReference, RFdiffusionAdapter, rank


def main(design_dir: str) -> None:
    verifier = GeometryVerifier(PDBReference())
    designs = RFdiffusionAdapter(design_dir, pdb_element="NI", metal_label="Ni2+").design("Ni2+")
    print(f"# touchstone — {len(designs)} nickel designs\n")
    print(f"{'design':12} {'CN':>3} {'score':>9}  verdict")
    for d, v in rank(designs, verifier):
        flag = "DEFER" if v.ood else ("trust" if v.trust else "weak")
        print(f"{d.sequence:12} {d.site.coordination_number:>3} {v.score:9.2e}  {flag:6} {v.reason}")


if __name__ == "__main__":
    typer.run(main)
