"""Paths, kit loading and logging shared by the holdout scripts. Shipped next to them on the training host."""

import json
import os
import shutil
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(os.environ.get("CASMI_HOLDOUT_ROOT", Path.home() / "casmi_holdout"))
TRAIN = Path(os.environ.get("CASMI_TRAIN", Path.home() / "casmi_seeds" / "data" / "train.parquet"))
KIT_SOURCE = Path(os.environ.get("CASMI_KIT", Path.home() / "casmi_seeds" / "kit"))
KIT, WORK, RUNS = ROOT / "kit", ROOT / "work", ROOT / "runs"


def log(*a: object) -> None:
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def load_kit() -> None:
    """Stage the public kit under ROOT once and make its modules importable."""
    if not KIT.exists():
        shutil.copytree(KIT_SOURCE, KIT)
    if str(KIT) not in sys.path:
        sys.path.insert(0, str(KIT))


def kit_config(**overrides: Any) -> dict[str, Any]:
    cfg = json.loads((KIT / "config.json").read_text())
    cfg.update(work=str(WORK), train_path=str(TRAIN), external_paths=[], **overrides)
    return cfg
