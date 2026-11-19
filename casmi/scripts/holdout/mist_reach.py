"""Reach check for the MIST msg_aug ranks 2-25 blend, replayed on the natural-product holdout's eval lists.

Kit order (eval_scores.parquet) stands in for the pipeline: rank 1 is kept, ranks 2-25 are reordered by
0.5 z(-log rank) + 0.5 z(MIST cosine). Spectra: train.parquet rows of each eval molecule, timsTOF only when present,
positive mode only (MIST has no negative ions); the formula is the true one (same-formula lists).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import torch
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

sys.path.insert(0, str(Path.home() / "mist_probe"))
from mist.data.featurizers import assign_subformulae
from mist.data.input_signal import collate_spectrum_only
from mist.models.scaled_mist import ScaledMIST
from mist.utils import chem_utils

RDLogger.DisableLog("rdApp.*")
H = Path.home() / "casmi_holdout_np"
DEV = torch.device("cuda")
MAX_PEAKS, PPM, INSTRUMENT, K = 10, 20.0, "Bruker Q-ToF (LCMS)", 25
ION_ALIAS = {"[M+NH4]+": "[M+H3N+H]+"}
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=4096, includeChirality=False)


def featurize(mz, it, prec, formula, ion):
    try:
        return encode(mz, it, prec, formula, ion)
    except ValueError:
        return None


def encode(mz, it, prec, formula, ion):
    mz, it = np.asarray(mz, float), np.asarray(it, float)
    keep = (mz <= prec + 1) & (it > 0)
    if not keep.any():
        return None
    mz, it = mz[keep], np.sqrt(it[keep] / it[keep].max())
    order = np.argsort(it)[::-1][:60]
    mz, it = mz[order], it[order]
    mz, it = mz[it > 0.001], it[it > 0.001]
    tree = assign_subformulae({"cand_form": formula, "cand_ion": ion,
                               "output_tbl": {"mz": mz.tolist(), "ms2_inten": it.tolist()}}, PPM)
    table = tree["output_tbl"] or {"formula": [], "ms2_inten": [], "ions": []}
    top = np.argsort(np.asarray(table["ms2_inten"], float))[::-1][:MAX_PEAKS]
    frags = [table["formula"][i] for i in top]
    width, n = MAX_PEAKS + 1, len(frags) + 1
    form = np.zeros((width, 18), np.int16)
    form[:n] = [chem_utils.formula_to_dense(f) for f in [*frags, formula]]
    ionv = np.zeros(width, np.uint8)
    ionv[:n] = [chem_utils.get_ion_idx(x) for x in [*(table["ions"][i] for i in top), ion]]
    intens = np.zeros(width, np.float16)
    intens[:n] = [*(table["ms2_inten"][i] for i in top), 1.0]
    ptype = np.zeros(width, np.uint8)
    ptype[n - 1] = 3
    return dict(form_vec=form, ion_vec=ionv, intens=intens, peak_type=ptype, lengths=np.uint8(n),
                instruments=np.uint8(chem_utils.get_instr_idx(INSTRUMENT)), row_index=np.int64(0))


def z(x):
    s = x.std()
    return (x - x.mean()) / (s if s > 0 else 1.0)


ckpt = torch.load(Path.home() / "mist_probe/ckpt/mist_msg_aug.ckpt", map_location="cpu", weights_only=False)
model = ScaledMIST(**ckpt["config"]["model"])
model.load_state_dict(ckpt["ema_model"], strict=True)
model = model.to(DEV).eval()

mol = pd.read_parquet(H / "work/molecules.parquet").set_index("key")
lists = pd.read_parquet(H / "runs/lists.parquet", columns=["mol_key", "cand_key", "cand_smiles", "is_true"])
kit = pd.read_parquet(H / "runs/kit/eval_scores.parquet")
lists = lists.merge(kit[["mol_key", "cand_key", "score"]], on=["mol_key", "cand_key"])
keys = lists.mol_key.unique().tolist()
spec = ds.dataset(Path.home() / "casmi_seeds/data/train.parquet").to_table(
    columns=["inchikey14", "ms2_mzs", "ms2_normalized_intensities", "precursor_mz", "adduct", "ionization_mode",
             "instrument_type"], filter=ds.field("inchikey14").isin(keys)).to_pandas()
spec["ion"] = spec.adduct.map(lambda a: ION_ALIAS.get(a, a))
spec = spec[spec.ion.isin(chem_utils.ION_LST)]
tims = spec[spec.instrument_type.eq("timsTOF")]
spec = pd.concat([tims, spec[~spec.inchikey14.isin(set(tims.inchikey14))]])

pred = {}
with torch.inference_mode():
    for key, g in spec.groupby("inchikey14"):
        feats = [f for r in g.head(8).itertuples()
                 if (f := featurize(r.ms2_mzs, r.ms2_normalized_intensities, r.precursor_mz, mol.formula[key], r.ion))]
        if feats:
            batch = {k: v.to(DEV) for k, v in collate_spectrum_only(feats).items() if k != "row_indices"}
            pred[key] = model(batch)[0].float().cpu().numpy().mean(0)

rows, moved = [], []
for key, g in lists.groupby("mol_key"):
    g = g.sort_values("score", ascending=False, kind="stable").head(K).reset_index(drop=True)
    true = np.flatnonzero(g.is_true.to_numpy())
    base_rank = int(true[0]) + 1 if len(true) else None
    new = g
    if key in pred and len(g) > 2:
        fps = np.stack([np.asarray(GEN.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)), np.float32) for s in g.cand_smiles])
        p = pred[key]
        cos = fps @ p / (np.linalg.norm(fps, axis=1) * np.linalg.norm(p) + 1e-9)
        tail = np.arange(1, len(g))
        blend = 0.5 * z(-np.log(tail + 1.0)) + 0.5 * z(cos[tail])
        new = pd.concat([g.iloc[:1], g.iloc[tail[np.argsort(-blend, kind="stable")]]]).reset_index(drop=True)
    t2 = np.flatnonzero(new.is_true.to_numpy())
    new_rank = int(t2[0]) + 1 if len(t2) else None
    rr = lambda r: 0.0 if r is None else 1.0 / r
    rows.append(dict(key=key, lib=mol.library.get(key), scored=key in pred, base=rr(base_rank), new=rr(new_rank),
                     changed=not new.cand_key.equals(g.cand_key)))
    if base_rank and new_rank and base_rank != new_rank:
        moved.append(new_rank - base_rank)

r = pd.DataFrame(rows)
d = r.new - r.base
out = dict(molecules=len(r), scored=int(r.scored.sum()), changed_lists=int(r.changed.sum()),
           true_in_2_25=int(((r.base > 0) & (r.base < 1)).sum()), true_moved=len(moved),
           true_up=int(sum(m < 0 for m in moved)), true_down=int(sum(m > 0 for m in moved)),
           mrr_base=round(r.base.mean(), 4), mrr_new=round(r.new.mean(), 4), delta=round(d.mean(), 4),
           se=round(d.std(ddof=1) / np.sqrt(len(d)), 4),
           by_library={k: dict(n=len(g), delta=round((g.new - g.base).mean(), 4)) for k, g in r.groupby("lib")})
print("REACH", json.dumps(out), flush=True)
(H / "mist_reach.json").write_text(json.dumps(out, indent=2))
