"""Compare BindCraft's default vs relaxed filters, and surface hotspot-related knobs.

Hotspot drift is BindCraft's main weakness for this task: its first design reached AF2
i_ptm 0.79 while contacting mature 117-124 instead of the specified 113/115/143-147.
"""
import json
import sys
from pathlib import Path

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home() / "BindCraft"


def fmt(v):
    if isinstance(v, dict):
        t = v.get("threshold")
        hi = v.get("higher")
        if t is None:
            return "none"
        op = ">" if hi else "<"
        return op + str(t)
    return str(v)


d = json.load((root / "settings_filters/default_filters.json").open())
r = json.load((root / "settings_filters/relaxed_filters.json").open())

print("filters that DIFFER between default and relaxed:")
print("  %-36s %-18s %-18s" % ("filter", "default", "relaxed"))
for k in sorted(set(d) | set(r)):
    a, b = fmt(d.get(k)), fmt(r.get(k))
    if a != b:
        print("  %-36s %-18s %-18s" % (k, a, b))

print("\nfilters mentioning hotspot / interface / contact (default set):")
for k in sorted(d):
    if any(t in k.lower() for t in ("hotspot", "interface", "contact", "shape")):
        print("  %-40s %s" % (k, fmt(d[k])))

adv = json.load((root / "settings_advanced/default_4stage_multimer_hardtarget.json").open())
print("\nadvanced knobs relevant to hotspot adherence:")
for k in sorted(adv):
    if any(t in k.lower() for t in ("con", "inter", "hotspot", "weight")):
        print("  %-40s %s" % (k, adv[k]))
