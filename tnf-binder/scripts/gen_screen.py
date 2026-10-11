"""Emit the config screen for lifting design_to_target_iptm.

Every architecture so far plateaued at d2t_iptm 0.11-0.30 against a ~0.5 bar. Rather than scale
a config that does not work, screen the untested levers on a small budget and scale the winner.

Levers, in expected order of impact:
  epitope      The original hotspots spanned 8 residues, two protomers and 18 A -- far more than
               the 3-5 tight hotspots that are standard practice, and centred on the LEAST
               designable conserved cluster (concavity 81, 1 apolar). The 113/115/143-147 pocket
               is more concave (104) and more hydrophobic (3 apolar), still 100% conserved, and
               still a genuine TNFR2 contact site.
  helix        All-helical binders are markedly easier to design and to re-predict confidently
               than mixed or beta topologies.
  length       A shallow groove may need more binder to grip it.

Indices in the emitted YAML are 1-based CHAIN POSITIONS. The target files keep 1TNF numbering
(6..157), so position = mature - 5. Asserted against the expected residue identities below --
an unguarded version of this arithmetic already produced a silent five-residue shift once.
"""
import sys
from pathlib import Path

from Bio import PDB
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE2ONE

SD, OUT = Path(sys.argv[1]), Path(sys.argv[2])
TARGET = "/home/ubuntu/tnf/targets/tnf_trimer.pdb"

st = PDB.PDBParser(QUIET=True).get_structure("t", SD / "tnf_trimer.pdb")[0]

pos_of, aa_of = {}, {}
for ch in st:
    for i, r in enumerate((x for x in ch if x.id[0] == " "), start=1):
        pos_of[(ch.id, r.id[1])] = i
        aa_of[(ch.id, r.id[1])] = THREE2ONE.get(r.get_resname(), "X")

# Candidate epitopes, mature numbering, with the identities we expect to find there.
EPITOPES = {
    # concave hydrophobic conserved pocket -- the designability winner
    "pocket": ({"A": [113, 115, 143, 144]}, {"A": "PYLF"}),
    # the original lobe, narrowed from 8 spread residues to 3 tight ones
    "lobe3":  ({"A": [88, 90, 92]}, {"A": "QKN"}),
    # pocket plus two more conserved contacts, still one protomer
    "pocket5": ({"A": [67, 113, 115, 143, 144]}, {"A": "QPYLF"}),
}

def ranges(nums: list[int]) -> str:
    nums, out, i = sorted(nums), [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"{nums[i]}..{nums[j]}" if j > i else f"{nums[i]}")
        i = j + 1
    return ",".join(out)

def binding_block(spec: dict[str, list[int]], expect: dict[str, str]) -> str:
    lines = []
    for ch, nums in spec.items():
        got = "".join(aa_of[(ch, n)] for n in nums)
        want = expect[ch]
        assert got == want, f"numbering guard failed on chain {ch}: {got!r} != {want!r}"
        r = ranges([pos_of[(ch, n)] for n in nums])
        print(f"   {ch}: mature {nums} -> positions {r}  identities {got!r} OK")
        lines.append(f'        - chain: {{id: {ch}, binding: "{r}"}}')
    return "\n".join(lines)

def helix_bundle(n: int) -> str:
    """A three-helix bundle pattern of exactly n residues.

    Uniform 'h' is not a designable topology -- a real bundle needs loops to turn between
    helices. Three helices with two short loops is the classic de novo minibunder fold and is
    much easier to re-predict confidently than a mixed or beta topology.
    """
    loop = 4
    helix_total = n - 2 * loop
    a = helix_total // 3
    b = helix_total - 2 * a
    return "h" * a + "l" * loop + "h" * b + "l" * loop + "h" * a


CONFIGS = [
    # name,            epitope,    length,    ss builder
    ("pocket_free",    "pocket",   "60..95",  None),
    ("pocket_helix",   "pocket",   "60..95",  helix_bundle),
    ("pocket5_helix",  "pocket5",  "70..110", helix_bundle),
    ("lobe3_free",     "lobe3",    "60..95",  None),
]

OUT.mkdir(parents=True, exist_ok=True)
for name, epi, length, ss in CONFIGS:
    spec, expect = EPITOPES[epi]
    print(f"\n=== {name}: epitope={epi} length={length} ss={ss or 'free'}")
    bb = binding_block(spec, expect)
    ss_line = ""
    ss_desc = "unconstrained"
    if ss:
        # Size to the MINIMUM of the length range: BoltzGen rejects a secondary_structure string
        # longer than the sampled sequence (its error says "shorter than the sequence", but the
        # check is len(spec) > n). Trailing positions default to 'u'.
        n_min = int(length.split("..")[0])
        pattern = ss(n_min)
        assert len(pattern) == n_min, f"ss pattern {len(pattern)} != {n_min}"
        ss_line = f'\n      secondary_structure: "{pattern}"'
        ss_desc = f"3-helix bundle over first {n_min} residues"
        print(f"   ss: {pattern}")
    text = f"""# TNF-alpha binder screen config: {name}
# epitope={epi}  length={length}  secondary_structure={ss_desc}
#
# Screening to lift design_to_target_iptm off the 0.11-0.30 plateau. Indices are 1-based CHAIN
# POSITIONS; target files keep 1TNF numbering so position = mature - 5, asserted by gen_screen.py.
entities:
  - protein:
      id: X
      sequence: {length}{ss_line}
  - file:
      path: {TARGET}
      include:
        - chain: {{id: A}}
        - chain: {{id: B}}
        - chain: {{id: C}}
      binding_types:
{bb}
      structure_groups: "all"
"""
    (OUT / f"screen_{name}.yaml").write_text(text)

print(f"\nwrote {len(CONFIGS)} screen specs to {OUT}")
