"""Pull binder sequences out of BoltzGen output complexes into a FASTA.

The binder is identified by low identity to human TNF-alpha, so it does not matter which
chain id BoltzGen assigned after its collision renaming.
"""
import argparse
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).parent))

HUMAN = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTH"
         "TISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")


def is_target(seq: str) -> bool:
    """TNF protomers are ~150 aa and share long exact stretches with the reference."""
    return any(seq[i:i + 25] in HUMAN for i in range(0, max(1, len(seq) - 25), 5))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("roots", type=Path, nargs="+")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min_len", type=int, default=10)
    ap.add_argument("--max_len", type=int, default=250)
    args = ap.parse_args()

    recs, seen = [], set()
    for root in args.roots:
        for f in sorted(root.rglob("*.cif")):
            try:
                st = gemmi.read_structure(str(f))
                st.setup_entities()
            except Exception as e:
                print(f"skip {f.name}: {e}", file=sys.stderr)
                continue
            for ch in st[0]:
                seq = gemmi.one_letter_code([r.name for r in ch]).upper()
                if not (args.min_len <= len(seq) <= args.max_len) or is_target(seq):
                    continue
                if seq in seen:
                    continue
                seen.add(seq)
                recs.append((f"{root.name}__{f.stem}__{ch.name}", seq))

    args.out.write_text("".join(f">{n}\n{s}\n" for n, s in recs))
    print(f"extracted {len(recs)} unique binder sequences -> {args.out}")


if __name__ == "__main__":
    main()
