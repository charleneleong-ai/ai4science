"""Extract VHH CDR loops, so the novelty gate can judge what was actually designed.

A grafted nanobody deliberately keeps a natural framework, so a whole-sequence memorisation
screen rejects every one of them. Only the CDRs are designed, and CDR3 dominates binding, so
that is what must be novel. Loops are located by the conserved framework anchors flanking
them rather than by residue index, because design_insertions make CDR lengths variable.
"""
import re

# Framework anchors immediately flanking each CDR in a VHH.
CDR1 = re.compile(r"C[A-Z]{3}(?P<cdr>[A-Z]{3,14})W[VIFLY][RK]Q")
CDR2 = re.compile(r"W[VIFLY][RK]Q[A-Z]{8,20}?(?P<cdr>[A-Z]{5,19})R[FL][TA][IV]S")
CDR3 = re.compile(r"[YFC][YFAVC]C(?P<cdr>[A-Z]{3,30})WG[QKRP]G")


def cdrs(seq: str) -> dict[str, str]:
    out = {}
    for name, rx in (("cdr1", CDR1), ("cdr2", CDR2), ("cdr3", CDR3)):
        m = rx.search(seq)
        out[name] = m.group("cdr") if m else ""
    return out


def is_vhh(seq: str) -> bool:
    """Needs the two most diagnostic framework motifs, not all three CDRs."""
    return bool(re.search(r"W[VIFLY][RK]Q", seq) and re.search(r"WG[QKRP]G", seq))


def read_fasta(path: str) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name and buf:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name and buf:
        seqs[name] = "".join(buf)
    return seqs


if __name__ == "__main__":
    import sys
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    as_fasta = "--fasta" in sys.argv   # emit CDRs as FASTA for the novelty gate
    # CDR1/CDR2 in a grafted VHH are germline-derived by construction and always match natural
    # sequences (e.g. the YYADSVKG junction), so gating on them rejects every valid design.
    # CDR3 carries the specificity and is where the de novo content actually is.
    only3 = "--cdr3" in sys.argv

    for f in args:
        for name, s in read_fasta(f).items():
            c = cdrs(s)
            if as_fasta:
                # Concatenate the loops: the gate only has to judge designed residues.
                joined = c["cdr3"] if only3 else "".join(c.values())
                if joined:
                    print(f">{name}\n{joined}")
            else:
                print(f"{name}\n  vhh={is_vhh(s)} len={len(s)}")
                for k, v in c.items():
                    print(f"  {k}: {v!r} ({len(v)} aa)")
