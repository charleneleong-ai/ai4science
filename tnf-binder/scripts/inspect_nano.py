"""Report CDR diversity across the nanobody designs."""
import glob
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).parent))
from cdr_extract import cdrs, is_vhh
from tnf_common import HUMAN


def is_tnf(s: str) -> bool:
    return any(s[i:i + 25] in HUMAN for i in range(0, max(1, len(s) - 25), 5))


rows = []
for f in sorted(glob.glob(sys.argv[1])):
    st = gemmi.read_structure(f)
    st.setup_entities()
    for ch in st[0]:
        s = gemmi.one_letter_code([r.name for r in ch])
        if len(s) < 160 and not is_tnf(s):
            rows.append((Path(f).stem, len(s), is_vhh(s), cdrs(s)))

print(f"{len(rows)} nanobody designs\n")
hdr = f"{'design':16s} {'len':>4} {'vhh':5s} {'CDR1':<13} {'CDR2':<20} CDR3"
print(hdr)
for n, L, v, c in rows:
    print(f"{n:16s} {L:>4} {str(v):5s} {c['cdr1']:<13} {c['cdr2']:<20} {c['cdr3']}")

for k in ("cdr1", "cdr2", "cdr3"):
    uniq = len({c[k] for _, _, _, c in rows if c[k]})
    print(f"unique {k}: {uniq}/{len(rows)}")
