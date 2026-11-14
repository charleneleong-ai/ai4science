"""Prepare the kit cache once and replace its random holdout with a structure-disjoint one.

Usage: build_holdout.py. Splits written to work/molecules.parquet and work/spectra.parquet:
- eval: up to 750 Enveda-180 molecules plus every enveda-np-examples molecule, scored by eval_holdout.py only;
- validation: 1,000 more Enveda-180 molecules, used by the kit for early stopping only;
- excluded: any other molecule with Morgan Tanimoto >= 0.7 to an eval or validation molecule;
- train: everything else.
molecules.parquet also gains a `library` column (the molecule's ingest_lib in train.parquet).
"""

import hashlib
import json

import numpy as np
import pandas as pd

from common import ROOT, TRAIN, WORK, kit_config, load_kit, log

N_EVAL_ENVEDA, N_VAL, TANIMOTO_MAX, CHUNK = 750, 1000, 0.7, 16


def prepare_cache() -> None:
    if (WORK / "molecules.parquet").exists():
        log("CACHE_EXISTS")
        return
    load_kit()
    import prepare

    if "normalized_smiles" not in prepare.ALIASES["smiles"]:
        prepare.ALIASES["smiles"].append("normalized_smiles")
    audit = prepare.prepare(kit_config())
    log("PREP_DONE", audit["spectra"], audit["molecules"])


def ranked(keys: list[str], salt: str) -> list[str]:
    return sorted(keys, key=lambda k: hashlib.sha256((salt + k).encode()).hexdigest())


def choose_split(mol: pd.DataFrame) -> pd.Series:
    enveda = ranked(mol.key[mol.library.eq("enveda-180")].tolist(), "holdout-20261114")
    split = pd.Series("train", index=mol.index)
    split[mol.library.eq("enveda-np-examples") | mol.key.isin(enveda[:N_EVAL_ENVEDA])] = "eval"
    split[mol.key.isin(enveda[N_EVAL_ENVEDA : N_EVAL_ENVEDA + N_VAL])] = "validation"
    return split


def max_tanimoto_to(fps: np.ndarray, held: np.ndarray) -> np.ndarray:
    """Max Morgan Tanimoto of every molecule to any held molecule, on bit-packed fingerprints."""
    words = np.packbits(fps, axis=1).view(np.uint64)
    counts = np.bitwise_count(words).sum(1)
    best = np.zeros(len(words))
    for start in range(0, len(held), CHUNK):
        h = held[start : start + CHUNK]
        inter = np.bitwise_count(words[h][:, None, :] & words[None]).sum(-1)
        sim = inter / np.maximum(counts[h][:, None] + counts[None] - inter, 1)
        best = np.maximum(best, sim.max(0))
        if start % (CHUNK * 16) == 0:
            log("NEIGHBOURS", start, "/", len(held))
    return best


def main() -> None:
    prepare_cache()
    mol = pd.read_parquet(WORK / "molecules.parquet")
    spec = pd.read_parquet(WORK / "spectra.parquet")
    lib = pd.read_parquet(TRAIN, columns=["inchikey14", "ingest_lib"]).drop_duplicates("inchikey14")
    mol["library"] = mol.key.map(dict(zip(lib.inchikey14, lib.ingest_lib)))
    split = choose_split(mol)
    near = max_tanimoto_to(np.load(WORK / "fingerprints.npy"), np.where(split.ne("train"))[0]) >= TANIMOTO_MAX
    split[near & split.eq("train").to_numpy()] = "excluded"
    mol["split"] = split.to_numpy()
    spec["split"] = spec.mol_idx.map(mol.split)
    mol.to_parquet(WORK / "molecules.parquet", index=False)
    spec.to_parquet(WORK / "spectra.parquet", index=False)
    report = dict(
        molecules={k: int(v) for k, v in mol.split.value_counts().items()},
        eval_by_library={k: int(v) for k, v in mol.library[mol.split.eq("eval")].value_counts().items()},
        spectra={k: int(v) for k, v in spec.split.value_counts().items()},
        tanimoto_max=TANIMOTO_MAX,
        unused_kit_keys="val_fraction and seed no longer choose the holdout",
    )
    (ROOT / "holdout_audit.json").write_text(json.dumps(report, indent=2))
    log("SPLIT_DONE", json.dumps(report))


if __name__ == "__main__":
    main()
