"""Confirm BindCraft2 resolves the AlphaFold weights we already have on disk.

BINDCRAFT_AF2_PARAMS points at the BindCraft 1 download (same params_2022-12-06 set), so the
5.3 GB is shared rather than fetched twice -- the volume is at 99%.
"""
import os
import sys

os.environ.setdefault("BINDCRAFT_AF2_PARAMS", os.path.expanduser("~/BindCraft/params"))
sys.path.insert(0, os.path.expanduser("~/BindCraft2"))

from bindcraft import model_weights as mw

params = mw.alphafold_parameters(download=False)
print("af2 params resolved to:", params or "(none)")
mpnn = mw.proteinmpnn_weights()
print("proteinmpnn weights:", mpnn)
problems = mw.missing_model_weights(params, mpnn) if params else ("no af2 params",)
print("missing:", list(problems)[:4] if problems else "none - all weights present")
