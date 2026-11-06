"""Cut the design and verification inputs out of 7OW6 and 7OW5.

Four targets come out of this, because the challenge's three ranking tiers are three separate
assay conditions and generation wants a fifth, cheaper view:

- `pmhc_g12d.pdb`  the intact complex, heavy chain and beta-2 microglobulin and the G12D
  peptide. The generation target. An earlier version handed the generator the groove alone,
  alpha1+alpha2 with the peptide, on the reasoning that alpha3 and beta-2 microglobulin sit
  under the floor and only cost sampling budget. Both designs it returned were beta-2
  microglobulin at 86-91% identity: given the fold with its own partner missing, the generator
  supplied the partner. Every natural interface is therefore kept occupied, leaving the groove
  top as the only vacancy.
- `complex_g12d`, `complex_wt`, `apo_hla`  sequence sets for co-folding, matching the three
  things the lab measures: the G12D complex, the same HLA carrying the wild-type peptide, and
  the peptide-free HLA with beta-2 microglobulin.

Hotspots are peptide chain positions, chosen from measured exposure and reachability rather than
from the mutation's identity: `epitope_geometry.py` finds the position-6 carboxylate unreachable
(a 0.02% escape cone over OD1), so the variants below centre on the exposed backbone bulge at
positions 5-8, where the wild-type complex is about twice as open at position 5 and so where a
tightly packed binder has its best chance of incidental discrimination.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import gemmi
import typer
from common import (
    PEPTIDE_G12D,
    PEPTIDE_WT,
    VARIANTS,
    chain_sequence,
    load,
    one_letter,
)

app = typer.Typer(add_completion=False)

HLA_CHAIN = "A"
B2M_CHAIN = "B"
PEPTIDE_CHAIN = "C"

MUTANT_ID = "7OW6"


# Expected chain lengths in the carved complex, so a renamed or re-deposited chain fails here
# rather than after nine GPU-hours against the wrong surface.
EXPECTED_LENGTHS = {HLA_CHAIN: 276, B2M_CHAIN: 99, PEPTIDE_CHAIN: 10}

# The T-cell receptor in the ternary deposition. It is not part of the target, but it is the
# natural binder of the groove top and so the likeliest thing for a generator to reproduce
# instead of designing; its chains go into the novelty references for that reason. The same
# two chains serve as the one positive control available on this epitope, being a protein
# known to bind it and to discriminate the mutation.
TCR_CHAINS = ("D", "E")



def sequence_of(structure: gemmi.Structure, chain: str) -> str:
    return chain_sequence(structure[0][chain])


def check_chains(structure: gemmi.Structure, expected: dict[str, int]) -> None:
    """Confirm the deposition's chains are the ones the hardcoded IDs assume."""
    found = {ch.name: len(ch) for ch in structure[0]}
    for name, length in expected.items():
        if found.get(name) != length:
            raise typer.BadParameter(
                f"chain {name} has {found.get(name)} residues, expected {length}; "
                f"chains present: {found}"
            )


def check_same_receptor(mutant: str, wildtype: str) -> None:
    """The two depositions must carry the same allele, differing only at the construct ends.

    A difference inside the groove would mean the structures are not the matched pair the
    comparison assumes, which no downstream score could reveal.
    """
    if mutant == wildtype:
        return
    if not (mutant.endswith(wildtype) or wildtype.endswith(mutant)
            or mutant.startswith(wildtype) or wildtype.startswith(mutant)):
        raise typer.BadParameter(
            "the two heavy chains differ internally, not just at the construct ends; "
            "the depositions are not a matched pair"
        )


def copied_chain(source: gemmi.Chain, name: str) -> gemmi.Chain:
    chain = gemmi.Chain(name)
    for residue in source:
        chain.add_residue(residue)
    return chain


def write_target(source: Path, out: Path) -> dict[str, int]:
    """The intact pMHC as a three-chain PDB, with the peptide relabelled P for the hotspot list.

    Built by selection rather than deletion: gemmi 0.6.5 has no Chain.remove_residue, and
    assembling the wanted residues states the intent more directly than removing the rest.
    """
    model = load(source)[0]
    out_structure = gemmi.Structure()
    out_structure.add_model(gemmi.Model("1"))
    out_structure[0].add_chain(copied_chain(model[HLA_CHAIN], HLA_CHAIN))
    out_structure[0].add_chain(copied_chain(model[B2M_CHAIN], B2M_CHAIN))
    out_structure[0].add_chain(copied_chain(model[PEPTIDE_CHAIN], "P"))
    out_structure.setup_entities()
    counts = {ch.name: len(ch) for ch in out_structure[0]}
    out.parent.mkdir(parents=True, exist_ok=True)
    out_structure.write_pdb(str(out))
    return counts


def verify_hotspots(target: Path, chain: str, hotspots: tuple[int, ...], peptide: str) -> str:
    """Read hotspot positions back off the written file and confirm what they spell.

    Hotspots are chain positions, so one unmodelled residue shifts every one of them by a
    place; BoltzGen then completes normally having designed against the wrong surface, and
    nothing downstream can tell. 7OW4's peptide is missing positions 5 and 6 entirely, which
    is exactly how a carved file acquires a gap, and tnf lost a campaign to the silent version
    of this. The check reads the artefact the generator consumes, not what we meant to write.
    """
    residues = list(load(target)[0][chain])
    if len(residues) != len(peptide):
        raise typer.BadParameter(
            f"{target} chain {chain} has {len(residues)} residues, expected {len(peptide)}"
        )
    spelled = "".join(one_letter(residues[h - 1].name) for h in hotspots)
    expected = "".join(peptide[h - 1] for h in hotspots)
    if spelled != expected:
        raise typer.BadParameter(
            f"hotspots {hotspots} in {target} chain {chain} spell {spelled!r}, "
            f"not {expected!r} as {peptide} requires"
        )
    return spelled


def read_exposure(path: Path, chain: str) -> dict[int, float]:
    """Accessible area per peptide position, read from target_sites.py's own output.

    Transcribing these by hand left the spec comments asserting areas that a re-measurement
    would have changed, with nothing to detect the divergence.
    """
    with path.open() as handle:
        return {
            int(row["position"]): float(row["sasa_total"])
            for row in csv.DictReader(handle)
            if row["structure"] == MUTANT_ID and row["chain"] == chain
        }


def spec_for(
    name: str,
    hotspots: tuple[int, ...],
    target: Path,
    length: str,
    exposure: dict[int, float],
    spelled: str,
) -> str:
    binding = ",".join(str(h) for h in hotspots)
    noted = " ".join(f"P{h}={exposure[h]:.0f}A^2" for h in hotspots)
    return f"""# KRAS G12D pMHC binder screen: {name}
# Hotspots are 1-BASED PEPTIDE CHAIN POSITIONS in {target.name}, chain P, which holds
# {PEPTIDE_G12D} at positions 1-10; read back off that file, these spell {spelled}.
# Accessible area per hotspot, receptor removed: {noted}
# The mutated side chain is not a hotspot. Its only water-accessible oxygen, OD1, loses all
# accessible area against a 2.0 A probe and sits at the end of a funnel that admits a few
# percent of approach directions at 5-6 A depth, so a designed side chain cannot engage it.
# These variants target the backbone bulge the mutation reshapes, which is also where the one
# receptor known to tell these two peptides apart makes every one of its position-6 contacts.
entities:
  - protein:
      id: X
      sequence: {length}
  - file:
      path: {target}
      include:
        - chain: {{id: {HLA_CHAIN}}}
        - chain: {{id: {B2M_CHAIN}}}
        - chain: {{id: P}}
      binding_types:
        - chain: {{id: P, binding: "{binding}"}}
      structure_groups: "all"
"""


@app.command()
def main(
    data: Path = typer.Option(Path("data")),
    targets: Path = typer.Option(Path("targets")),
    specs: Path = typer.Option(Path("specs")),
    mutant: str = typer.Option(MUTANT_ID),
    wildtype: str = typer.Option("7OW5"),
    length: str = typer.Option("60..95", help="binder length range handed to the generator"),
    exposure_csv: Path = typer.Option(
        Path("peptide_exposure.csv"), help="target_sites.py output; areas go into the specs"
    ),
) -> None:
    mut, wt = load(data / f"{mutant}_pmhc.cif"), load(data / f"{wildtype}_pmhc.cif")
    check_chains(mut, EXPECTED_LENGTHS)
    hla, b2m = sequence_of(mut, HLA_CHAIN), sequence_of(mut, B2M_CHAIN)
    peptide_g12d, peptide_wt = sequence_of(mut, PEPTIDE_CHAIN), sequence_of(wt, PEPTIDE_CHAIN)
    for label, found, expected in (
        ("mutant", peptide_g12d, PEPTIDE_G12D),
        ("wild type", peptide_wt, PEPTIDE_WT),
    ):
        if found != expected:
            raise typer.BadParameter(f"{label} peptide is {found}, expected {expected}")
    check_same_receptor(hla, sequence_of(wt, HLA_CHAIN))

    counts = write_target(data / f"{mutant}_pmhc.cif", targets / "pmhc_g12d.pdb")
    typer.echo(f"pmhc_g12d.pdb chains {counts}")

    # One HLA and one beta-2 microglobulin across all three conditions. The depositions differ
    # by a leading construct residue on the heavy chain, and carrying that difference into the
    # co-fold would put a sequence change in the very comparison selectivity is read from.
    cofold = [
        {"name": "complex_g12d", "chains": {"A": hla, "B": b2m, "C": peptide_g12d}},
        {"name": "complex_wt", "chains": {"A": hla, "B": b2m, "C": peptide_wt}},
        {"name": "apo_hla", "chains": {"A": hla, "B": b2m}},
    ]

    out = targets / "cofold_targets.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cofold, indent=2) + "\n")
    for target in cofold:
        lengths = {k: len(v) for k, v in target["chains"].items()}
        typer.echo(f"  {target['name']:14s} {lengths}")

    ternary = load(data / f"{mutant}.cif")
    tcr = {c: sequence_of(ternary, c) for c in TCR_CHAINS}
    refs = {"hla_heavy": hla, "b2m": b2m, "peptide_g12d": peptide_g12d, "peptide_wt": peptide_wt}
    refs.update({f"tcr_{c}": s for c, s in tcr.items()})
    (targets / "novelty_references.json").write_text(json.dumps(refs, indent=2) + "\n")
    typer.echo(f"  novelty references: {', '.join(f'{k}({len(v)})' for k, v in refs.items())}")

    (targets / "controls.json").write_text(json.dumps({"tcr": tcr}, indent=2) + "\n")

    specs.mkdir(parents=True, exist_ok=True)
    target_path = (targets / "pmhc_g12d.pdb").resolve()
    exposure = read_exposure(exposure_csv, PEPTIDE_CHAIN)
    for name, hotspots in VARIANTS.items():
        spelled = verify_hotspots(target_path, "P", hotspots, PEPTIDE_G12D)
        if missing := [h for h in hotspots if h not in exposure]:
            raise typer.BadParameter(
                f"{exposure_csv} has no measured area for position(s) {missing}; "
                "re-run target_sites.py"
            )
        (specs / f"{name}.yaml").write_text(
            spec_for(name, hotspots, target_path, length, exposure, spelled)
        )
        typer.echo(f"  spec {name:14s} hotspots {hotspots} spell {spelled}")
    typer.echo(f"wrote {out} and {len(VARIANTS)} specs")


if __name__ == "__main__":
    app()
