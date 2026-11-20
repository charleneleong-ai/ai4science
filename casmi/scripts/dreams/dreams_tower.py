"""A molecule tower for the frozen official DreaMS embedding model, trained on train.parquet, and its reach check.

Usage (on the A100, from ~/dreams):
  dreams_tower.py embed    mean DreaMS embedding per molecule (positive mode, <= 4 spectra, timsTOF first)
  dreams_tower.py train    MLP from 4,096-bit Morgan r2 fingerprints into the DreaMS space, symmetric InfoNCE
  dreams_tower.py reach    replay ranks 2-25 of the natural-product holdout's kit order with the blend

The DreaMS weights are `DreaMS_embedding_model_torchscript.pt` from HF roman-bushuiev/DreaMS (MIT), never updated.
Molecules come from ~/casmi_holdout_np/work/molecules.parquet: the tower trains on split == train only, early-stops
on split == validation, and the reach replay uses split == eval, whose near neighbours are already excluded.
"""

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import torch
import torch.nn.functional as F
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
HOME = Path.home()
H = HOME / "casmi_holdout_np"
OUT = HOME / "dreams"
MODEL = OUT / "model/DreaMS_embedding_model_torchscript.pt"
TRAIN = HOME / "casmi_seeds/data/train.parquet"
N_PEAKS, PREC_INTENS, PER_MOL, BATCH, FP_BITS, K = 100, 1.1, 4, 384, 4096, 25
DEV = torch.device("cuda" if torch.cuda.is_available() else "cpu")
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=FP_BITS, includeChirality=False)


def log(*a: object) -> None:
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def load_dreams(device: torch.device = DEV) -> torch.jit.ScriptModule:
    """The TorchScript file was exported with its device frozen as CPU; point every `_device` at the target."""
    m = torch.jit.load(str(MODEL), map_location=device).eval()
    for mod in m.modules():
        if hasattr(mod, "_device"):
            mod._device = device
    return m


def preprocess(mzs: np.ndarray, ints: np.ndarray, prec: float) -> np.ndarray:
    """DreaMS's SpectrumPreprocessor: top 100 peaks in m/z order, zero-padded, base-peak relative, precursor first."""
    mzs, ints = np.asarray(mzs, np.float64), np.asarray(ints, np.float64)
    top = np.sort(np.argsort(ints, kind="stable")[-N_PEAKS:])
    spec = np.zeros((N_PEAKS + 1, 2), np.float32)
    spec[0] = (prec, PREC_INTENS)
    if len(top) and ints[top].max() > 0:
        spec[1 : len(top) + 1, 0] = mzs[top]
        spec[1 : len(top) + 1, 1] = ints[top] / ints[top].max()
    return spec


def embed_spectra(model: torch.jit.ScriptModule, specs: np.ndarray, bs: int = BATCH) -> np.ndarray:
    out = np.empty((len(specs), 1024), np.float32)
    with torch.inference_mode():
        for i in range(0, len(specs), bs):
            out[i : i + bs] = F.normalize(model(torch.from_numpy(specs[i : i + bs]).to(DEV)), dim=1).cpu().numpy()
    return out


def fingerprint(smiles: str) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles)
    return np.packbits(GEN.GetFingerprintAsNumPy(mol) if mol else np.zeros(FP_BITS, np.uint8))


class Tower(torch.nn.Module):
    def __init__(self, hidden: int = 2048, dim: int = 1024, dropout: float = 0.1) -> None:
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(FP_BITS, hidden), torch.nn.GELU(), torch.nn.Dropout(dropout),
            torch.nn.Linear(hidden, hidden), torch.nn.GELU(), torch.nn.Dropout(dropout),
            torch.nn.Linear(hidden, dim),
        )
        self.log_temp = torch.nn.Parameter(torch.tensor(np.log(1 / 0.05), dtype=torch.float32))

    def forward(self, fp: torch.Tensor) -> torch.Tensor:
        return F.normalize(self.net(fp), dim=1)


def molecules() -> pd.DataFrame:
    return pd.read_parquet(H / "work/molecules.parquet", columns=["key", "smiles", "split", "library"])


def embed() -> None:
    mol = molecules()
    keep = mol[mol.split.isin(["train", "validation", "eval"])]
    spec = ds.dataset(TRAIN).to_table(
        columns=["inchikey14", "ms2_mzs", "ms2_normalized_intensities", "precursor_mz", "instrument_type"],
        filter=(ds.field("ionization_mode") == "positive") & ds.field("inchikey14").isin(keep.key.tolist()),
    ).to_pandas()
    spec["tims"] = spec.instrument_type.eq("timsTOF")
    spec = spec.sort_values(["inchikey14", "tims"], ascending=[True, False], kind="stable")
    spec = spec.groupby("inchikey14", sort=False).head(PER_MOL).reset_index(drop=True)
    log("SPECTRA", len(spec), "for", spec.inchikey14.nunique(), "molecules")
    arr = np.stack([preprocess(r.ms2_mzs, r.ms2_normalized_intensities, r.precursor_mz) for r in spec.itertuples()])
    t0 = time.time()
    emb = embed_spectra(load_dreams(), arr)
    log("EMBEDDED", len(arr), f"in {time.time() - t0:.0f}s")
    codes, keys = pd.factorize(spec.inchikey14)
    sums = np.zeros((len(keys), 1024), np.float32)
    np.add.at(sums, codes, emb)
    np.save(OUT / "mol_emb.npy", sums / np.linalg.norm(sums, axis=1, keepdims=True))
    pd.Series(keys, name="key").to_frame().to_parquet(OUT / "mol_emb_keys.parquet", index=False)


def training_data() -> tuple[pd.DataFrame, np.ndarray]:
    keys = pd.read_parquet(OUT / "mol_emb_keys.parquet").key
    emb = np.load(OUT / "mol_emb.npy")
    mol = molecules().set_index("key").loc[keys].reset_index()
    fp_path = OUT / "mol_fp.npy"
    if not fp_path.exists():
        with Pool(16) as p:
            np.save(fp_path, np.stack(p.map(fingerprint, mol.smiles, chunksize=2000)))
    return mol.assign(row=np.arange(len(mol))), emb


def to_dense(packed: np.ndarray) -> torch.Tensor:
    return torch.from_numpy(np.unpackbits(packed, axis=1)[:, :FP_BITS].astype(np.float32)).to(DEV)


def retrieval_mrr(tower: Tower, fp: np.ndarray, emb: np.ndarray) -> float:
    """MRR of each validation spectrum embedding against all validation molecules' tower outputs."""
    tower.eval()
    with torch.inference_mode():
        sim = torch.from_numpy(emb).to(DEV) @ tower(to_dense(fp)).T
        rank = (sim > sim.diag()[:, None]).sum(1) + 1
    tower.train()
    return float((1.0 / rank.float()).mean())


def train(epochs: int = 30, bs: int = 4096, lr: float = 3e-4) -> None:
    mol, emb = training_data()
    fp = np.load(OUT / "mol_fp.npy")
    tr, va = mol.row[mol.split.eq("train")].to_numpy(), mol.row[mol.split.eq("validation")].to_numpy()
    torch.manual_seed(20261119)
    tower = Tower().to(DEV)
    opt = torch.optim.AdamW(tower.parameters(), lr=lr, weight_decay=1e-4)
    steps = epochs * (len(tr) // bs)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=lr, total_steps=steps, pct_start=0.05)
    emb_t = torch.from_numpy(emb).to(DEV)
    best, history = -1.0, []
    for epoch in range(1, epochs + 1):
        perm = np.random.default_rng(epoch).permutation(tr)
        losses = []
        for i in range(0, len(perm) - bs + 1, bs):
            idx = perm[i : i + bs]
            logits = tower.log_temp.exp().clamp(max=100) * emb_t[idx] @ tower(to_dense(fp[idx])).T
            target = torch.arange(len(idx), device=DEV)
            loss = (F.cross_entropy(logits, target) + F.cross_entropy(logits.T, target)) / 2
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            sched.step()
            losses.append(loss.item())
        mrr = retrieval_mrr(tower, fp[va], emb[va])
        history.append(dict(epoch=epoch, loss=round(float(np.mean(losses)), 4), val_mrr=round(mrr, 4)))
        log("EPOCH", json.dumps(history[-1]))
        if mrr > best:
            best = mrr
            torch.save(dict(state=tower.state_dict(), epoch=epoch, val_mrr=mrr, fp_bits=FP_BITS), OUT / "tower.pt")
    (OUT / "train_history.json").write_text(json.dumps(dict(best_val_mrr=best, history=history), indent=2))
    log("TRAIN_DONE", json.dumps(dict(best_val_mrr=round(best, 4), train=len(tr), validation=len(va))))


def z(x: np.ndarray) -> np.ndarray:
    s = x.std()
    return (x - x.mean()) / (s if s > 0 else 1.0)


def reach() -> None:
    keys = pd.read_parquet(OUT / "mol_emb_keys.parquet").key
    emb = dict(zip(keys, np.load(OUT / "mol_emb.npy")))
    tower = Tower().to(DEV).eval()
    tower.load_state_dict(torch.load(OUT / "tower.pt", map_location=DEV)["state"])
    lib = molecules().set_index("key").library
    lists = pd.read_parquet(H / "runs/lists.parquet", columns=["mol_key", "cand_key", "cand_smiles", "is_true"])
    kit = pd.read_parquet(H / "runs/kit/eval_scores.parquet", columns=["mol_key", "cand_key", "score"])
    lists = lists.merge(kit, on=["mol_key", "cand_key"])
    rows, moved = [], []
    for key, g in lists.groupby("mol_key"):
        g = g.sort_values("score", ascending=False, kind="stable").head(K).reset_index(drop=True)
        new = g
        if key in emb and len(g) > 2:
            with torch.inference_mode():
                cand = tower(to_dense(np.stack([fingerprint(s) for s in g.cand_smiles]))).cpu().numpy()
            cos = cand @ emb[key]
            tail = np.arange(1, len(g))
            blend = 0.5 * z(-np.log(tail + 1.0)) + 0.5 * z(cos[tail])
            new = pd.concat([g.iloc[:1], g.iloc[tail[np.argsort(-blend, kind="stable")]]]).reset_index(drop=True)
        ranks = [int(np.flatnonzero(x.is_true.to_numpy())[0]) + 1 if x.is_true.any() else None for x in (g, new)]
        rr = [0.0 if r is None else 1.0 / r for r in ranks]
        rows.append(dict(key=key, lib=lib.get(key), scored=key in emb, base=rr[0], new=rr[1],
                         changed=not new.cand_key.equals(g.cand_key)))
        if None not in ranks and ranks[0] != ranks[1]:
            moved.append(ranks[1] - ranks[0])
    r = pd.DataFrame(rows)
    d = r.new - r.base
    out = dict(molecules=len(r), scored=int(r.scored.sum()), changed_lists=int(r.changed.sum()),
               true_in_2_25=int(((r.base > 0) & (r.base < 1)).sum()), true_moved=len(moved),
               true_up=sum(m < 0 for m in moved), true_down=sum(m > 0 for m in moved),
               mrr_base=round(r.base.mean(), 4), mrr_new=round(r.new.mean(), 4), delta=round(d.mean(), 4),
               se=round(d.std(ddof=1) / np.sqrt(len(d)), 4),
               by_library={k: dict(n=len(x), delta=round((x.new - x.base).mean(), 4)) for k, x in r.groupby("lib")})
    log("REACH", json.dumps(out))
    (OUT / "dreams_reach.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    {"embed": embed, "train": train, "reach": reach}[sys.argv[1]]()
