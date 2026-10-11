"""One-shot status for every running experiment.

Lives on the host and is invoked with no arguments, because inline python inside a
single-quoted ssh command kept tripping the local shell's globbing and quote handling.
"""
import csv
import glob
import os
import subprocess
from pathlib import Path

H = Path.home()
BC = H / "tnf/bindcraft_out/smoke"


def sh(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip()


print("=== GPU ===")
print(sh("nvidia-smi --query-gpu=memory.free,memory.used,utilization.gpu --format=csv,noheader"))
print(sh("nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader"))
print(f"disk free: {sh('df -h ~ | tail -1')}")

print("\n=== BoltzGen screen ===")
for c in ("pocket_free", "pocket_helix", "pocket5_helix", "lobe3_free"):
    d = H / f"tnf/out/screen_{c}"
    nd = len(glob.glob(str(d / "intermediate_designs/*.cif")))
    nf = len(glob.glob(str(d / "**/fold_out_npz/*.npz"), recursive=True))
    print(f"  {c:16s} design={nd:<3} folded={nf:<3}")

print("\n=== BindCraft ===")
# Sort by mtime, not name: "bcsmoke2_..." sorts BEFORE "bcsmoke_..." in ASCII
# ('2' < '_'), so alphabetical ordering kept reporting the older failed run.
logs = sorted(glob.glob(str(H / "tnf/logs/bcsmoke*.log")), key=os.path.getmtime)
if logs:
    txt = Path(logs[-1]).read_text(errors="replace")
    stages = [ln.strip() for ln in txt.splitlines()
              if ln.startswith(("Stage ", "Starting trajectory", "Trajectory ")) or
              "Accepted" in ln or "Rejected" in ln or "BINDCRAFT_RUN_DONE" in ln]
    for s in stages[-8:]:
        print(f"  {s[:140]}")
    for bad in ("RESOURCE_EXHAUSTED", "out of memory", "Traceback"):
        if bad in txt:
            print(f"  !! {bad} present in log")

for name in ("trajectory_stats.csv", "mpnn_design_stats.csv", "final_design_stats.csv"):
    f = BC / name
    if not f.exists():
        continue
    rows = list(csv.DictReader(f.open()))
    print(f"\n  --- {name}: {len(rows)} rows")
    keys = ["Design", "Length", "pLDDT", "pTM", "i_pTM", "i_pAE", "i_pLDDT",
            "Unrelaxed_Clashes", "Binder_Energy_Score"]
    for r in rows[:6]:
        print("   ", {k: r.get(k) for k in keys if k in r})

f = BC / "failure_csv.csv"
if f.exists():
    rows = list(csv.DictReader(f.open()))
    if rows:
        nz = {k: v for k, v in rows[0].items() if v not in ("0", "", "0.0")}
        print(f"\n  failure counters (nonzero): {nz if nz else 'none yet'}")

for sub in ("Accepted", "Rejected", "Trajectory/LowConfidence"):
    d = BC / sub
    if d.exists():
        print(f"  {sub}: {len(list(d.glob('*.pdb')))} pdb")
