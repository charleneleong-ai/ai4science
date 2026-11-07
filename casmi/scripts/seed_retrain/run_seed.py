"""Train the public sanpier v1 spectrum model with the shipped recipe, changing only the seed.

Usage: run_seed.py <seed>. Mirrors bgunshui/casmi26-scorer-train-public (the normalized_smiles alias fix to
prepare.py) but keeps the kit's own config.json schedule (16 x 500 steps), since the longer schedules scored lower on
the public LB (40 x 2000: 0.453 vs 0.472; 80 epochs alone: 0.460).
"""
import hashlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

ROOT = Path.home() / "casmi_seeds"
SEED = int(sys.argv[1])
RUN = ROOT / "runs" / str(SEED)
T0 = time.time()


def log(*a: object) -> None:
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


stage = RUN / "kit"
shutil.copytree(ROOT / "kit", stage, dirs_exist_ok=True)
prep = stage / "prepare.py"
src = prep.read_text()
OLD = "'smiles':['smiles','canonical_smiles']"
assert OLD in src
prep.write_text(src.replace(OLD, "'smiles':['smiles','canonical_smiles','normalized_smiles']"))
sys.path.insert(0, str(stage))
import prepare  # noqa: E402
import train  # noqa: E402

cfg = json.loads((ROOT / "kit" / "config.json").read_text())
cfg.update(seed=SEED, train_path=os.environ.get("CASMI_TRAIN", str(ROOT / "data" / "train.parquet")),
           work=str(RUN / "work"), output=str(RUN / "model"), external_paths=[])
if "CASMI_SMOKE" in os.environ:
    cfg.update(epochs=2, steps_per_epoch=20, max_validation_molecules=50)
log("CFG", json.dumps(cfg))

audit = prepare.prepare(cfg)
log("PREP_DONE", round(time.time() - T0), "s", json.dumps({k: audit[k] for k in ("spectra", "molecules", "rejected")}))
out = Path(cfg["output"])
out.mkdir(parents=True, exist_ok=True)
shutil.copy2(Path(cfg["work"]) / "data_audit.json", out / "data_audit.json")

result = train.run(cfg)
log("TRAIN_RESULT", json.dumps(result))
(out / "config_used.json").write_text(json.dumps(cfg, indent=2))
for name in ("model.py", "features.py", "inference.py", "config.json"):
    shutil.copy2(ROOT / "kit" / name, out / name)
for p in sorted(out.iterdir()):
    log("ARTIFACT", p.name, p.stat().st_size, hashlib.sha256(p.read_bytes()).hexdigest()[:16])
shutil.rmtree(cfg["work"], ignore_errors=True)
log(f"SEED_DONE {SEED} {round(time.time() - T0)}s")
