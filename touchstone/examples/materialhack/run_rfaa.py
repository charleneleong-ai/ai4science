"""Rank REAL RFdiffusionAA nickel designs with the geometry verifier.

RFdiffusionAA is run on a GPU box via apptainer (see materialhack/run_nickel.sh on
pi-a100-80gb). Point this at the directory of its output PDBs:

    python examples/materialhack/run_rfaa.py <rfaa_output_dir>

Raw RFdiffusionAA output is backbone-only — no sequence/sidechain design yet — so the
verifier should correctly DEFER on it: proper metal coordination needs the rest of the
pipeline (LigandMPNN sidechains). That deferral is the whole point — the verifier knows
the generator's output isn't a trustworthy metal site until it actually is.
"""

import sys

from touchstone import GeometryVerifier, RFdiffusionAdapter, rank


def main(output_dir: str) -> None:
    designs = RFdiffusionAdapter(output_dir, pdb_element="NI", metal_label="Ni2+").design("Ni2+")
    print(f"# RFdiffusionAA → touchstone — {len(designs)} nickel designs\n")
    print(f"{'design':12} {'CN':>3} {'score':>9}  verdict")
    for d, v in rank(designs, GeometryVerifier()):
        flag = "DEFER" if v.ood else ("trust" if v.trust else "weak")
        print(f"{d.sequence:12} {d.site.coordination_number:>3} {v.score:9.2e}  {flag:6} {v.reason}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "output/nickel")
