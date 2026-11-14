"""Train one kit variant on the structure-disjoint split, holdout stage only (no full refit).

Usage: train_variant.py <name> <seed> <epochs> <steps_per_epoch>. Reads the shared cache from build_holdout.py
and writes runs/<name>/model/: best_validation_model.pt plus the kit files `inference.Scorer` needs.
"""

import json
import shutil
import sys

from common import KIT, RUNS, kit_config, load_kit, log


def main(name: str, seed: int, epochs: int, steps: int) -> None:
    load_kit()
    import train

    out = RUNS / name / "model"
    cfg = kit_config(
        seed=seed, epochs=epochs, steps_per_epoch=steps, output=str(out), refit_full=False, holdout_hours=12.0
    )
    log("CFG", json.dumps(cfg))
    result = train.run(cfg)
    (out / "config_used.json").write_text(json.dumps(cfg, indent=2))
    for f in ("model.py", "features.py", "inference.py", "config.json"):
        shutil.copy2(KIT / f, out / f)
    log("VARIANT_DONE", name, json.dumps(result))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
