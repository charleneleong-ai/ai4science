"""Build a Boltz spec for the mouse TNF-alpha trimer.

Objective 2 wants one binder that recognises human AND mouse TNF-alpha. BindCraft2 supports
multi-target campaigns ("target": [human, mouse]), which optimises both interfaces directly
instead of designing against human and trusting that 79% identity carries the binder across.
That needs a mouse structure, and no mouse TNF-alpha crystal structure is being used here, so
predict the trimer from sequence.

Mouse soluble domain from UniProt P06804, aligned to the human construct boundaries.
"""
import sys
from pathlib import Path

MOUSE = ("LRSSSQNSSDKPVAHVVANHQVEEQLEWLSQRANALLANGMDLKDNQLVVPADGLYLVYSQVLFKGQGCPDYVLLTH"
         "TVSRFAISYQEKVNLLSAVKSPCPKDTPEGAELKPWYEPIYLGGVFQLEKGDQLSAEVNLPKYLDFAESGQVYFGVIAL")

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
chains = []
for cid in "ABC":
    chains.append(f"""  - protein:
      id: {cid}
      sequence: {MOUSE}
      msa: empty""")
(out / "mouse_tnf_trimer.yaml").write_text("version: 1\nsequences:\n" + "\n".join(chains) + "\n")
print(f"mouse soluble domain: {len(MOUSE)} aa")
print(f"wrote {out / 'mouse_tnf_trimer.yaml'} (3 protomers)")
