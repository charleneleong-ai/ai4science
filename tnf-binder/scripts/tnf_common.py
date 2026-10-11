"""Shared TNF-alpha sequence facts: human/mouse divergence and the pH-switch anchors."""
from Bio import Align

HUMAN = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTH"
         "TISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")
MOUSE = ("LRSSSQNSSDKPVAHVVANHQVEEQLEWLSQRANALLANGMDLKDNQLVVPADGLYLVYSQVLFKGQGCPDYVLLTH"
         "TVSRFAISYQEKVNLLSAVKSPCPKDTPEGAELKPWYEPIYLGGVFQLEKGDQLSAEVNLPKYLDFAESGQVYFGVIAL")

# Conserved, exposed, cationic residues in the TNFR2 epitope. A binder histidine placed near
# one of these is the lever for pH-dependent release: neutral at 7.4 (binds), cationic at 6.0
# (repelled). All three are human/mouse identical, so the switch does not cost cross-reactivity.
PH_ANCHORS = (32, 90, 128)


def divergence() -> dict[int, tuple[str, str]]:
    """Mature human resnum -> (human aa, mouse aa). Equal pair means conserved."""
    al = Align.PairwiseAligner(scoring="blastp", open_gap_score=-11, extend_gap_score=-1)
    ha, ma = (str(x) for x in al.align(HUMAN, MOUSE)[0])
    out, hi = {}, 0
    for hc, mc in zip(ha, ma):
        if hc != "-":
            hi += 1
            out[hi] = (hc, mc)
    return out


def divergent_positions() -> set[int]:
    return {p for p, (h, m) in divergence().items() if h != m}
