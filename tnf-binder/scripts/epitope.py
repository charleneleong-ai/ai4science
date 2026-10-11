"""Human/mouse TNF-alpha conservation vs the TNFR2 epitope (from 3ALQ)."""
import sys
from pathlib import Path
from Bio import PDB, Align
from Bio.PDB.SASA import ShrakeRupley

STRUCT_DIR = Path(sys.argv[1])

HUMAN = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTH"
         "TISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")
MOUSE = ("LRSSSQNSSDKPVAHVVANHQVEEQLEWLSQRANALLANGMDLKDNQLVVPADGLYLVYSQVLFKGQGCPDYVLLTH"
         "TVSRFAISYQEKVNLLSAVKSPCPKDTPEGAELKPWYEPIYLGGVFQLEKGDQLSAEVNLPKYLDFAESGQVYFGVIAL")

def align(a: str, b: str) -> tuple[str, str]:
    al = Align.PairwiseAligner(scoring="blastp", open_gap_score=-11, extend_gap_score=-1)
    return tuple(str(x) for x in al.align(a, b)[0])

def conservation(ha: str, ma: str) -> dict[int, tuple[str, str]]:
    """human resnum (1-based mature) -> (human aa, mouse aa or '-')"""
    out, hi = {}, 0
    for hc, mc in zip(ha, ma):
        if hc != "-":
            hi += 1
            out[hi] = (hc, mc)
    return out

def tnfr2_epitope(path: Path, cutoff: float = 4.5) -> dict[str, set[int]]:
    """TNF residues contacting TNFR2, split by which protomer they sit on."""
    st = PDB.PDBParser(QUIET=True).get_structure("x", path)[0]
    tnf = [c for c in st if c.id in "ABCDEF"]
    rec = [c for c in st if c.id in "RSTUVW"]
    ratoms = [a for c in rec for a in c.get_atoms() if a.element != "H"]
    ns = PDB.NeighborSearch(ratoms)
    hits: dict[str, set[int]] = {}
    for c in tnf:
        for res in c:
            if res.id[0] != " ":
                continue
            if any(ns.search(a.coord, cutoff) for a in res if a.element != "H"):
                hits.setdefault(c.id, set()).add(res.id[1])
    return hits

def sasa_rel(path: Path) -> dict[int, float]:
    st = PDB.PDBParser(QUIET=True).get_structure("x", path)[0]
    ShrakeRupley().compute(st, level="R")
    # chain A in trimer context => buried protomer interfaces already accounted
    return {r.id[1]: r.sasa for r in st["A"] if r.id[0] == " "}

ha, ma = align(HUMAN, MOUSE)
cons = conservation(ha, ma)
ident = sum(1 for h, m in cons.values() if h == m)
print(f"Soluble-domain identity: {ident}/{len(cons)} = {100*ident/len(cons):.1f}%\n")

ep = tnfr2_epitope(STRUCT_DIR / "3alq.pdb")
print("TNFR2-contacting TNF residues per chain (3ALQ, 4.5 A):")
for ch, s in sorted(ep.items()):
    print(f"  chain {ch}: n={len(s)}  {sorted(s)}")

# union across chains = the composite two-protomer receptor site
union = sorted(set().union(*ep.values())) if ep else []
sasa = sasa_rel(STRUCT_DIR / "1tnf.pdb")
print(f"\nComposite receptor epitope ({len(union)} unique positions), conservation + SASA in apo trimer:")
print(f"{'pos':>5} {'hum':>4} {'mou':>4} {'same':>5} {'SASA':>7}")
div_basic, cons_list, div_list = [], [], []
for p in union:
    h, m = cons.get(p, ("?", "?"))
    s = sasa.get(p, float("nan"))
    same = "yes" if h == m else "NO"
    print(f"{p:>5} {h:>4} {m:>4} {same:>5} {s:>7.1f}")
    (cons_list if h == m else div_list).append(p)
    if h in "RKH" and s > 20:
        div_basic.append((p, h, m, same))

print(f"\nConserved epitope positions ({len(cons_list)}): {cons_list}")
print(f"DIVERGENT epitope positions ({len(div_list)}): {[(p, cons[p][0], cons[p][1]) for p in div_list]}")
print(f"\nExposed basic residues in epitope (His-repulsion anchors for pH switch):")
for p, h, m, same in div_basic:
    print(f"  {h}{p}  mouse={m}  conserved={same}")
