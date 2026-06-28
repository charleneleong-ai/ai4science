"""Convert touchstone-selected BoltzGen winners into BoltzGen training targets (RAFT step 3).

Runs in the `bg` conda env on the A100 (needs `boltzgen` + the `~/.boltz/mols` CCD cache) —
it's the bridge from `rlvr_select`'s kept CIFs to the on-disk layout the trainer's
`DatasetConfig` reads:

    out/structures/<id>.npz   # Structure.dump  (from parse_mmcif → ParsedStructure.data)
    out/records/<id>.json     # Record.dump     (StructureInfo + per-chain ChainInfo)
    out/manifest.json         # Manifest over the kept ids → point DatasetConfig.manifest_path here

De-novo designs carry no MSA, so fine-tune with `use_msa: false` (single-sequence); the trainer
then skips MSA loading and no `msa_dir` files are needed. See docs/specs/2026-06-28-rlvr-boltzgen.md.

    BoltzGen generate → rlvr_select (winners) → winners_to_targets (this) → resume-train → repeat

    python scripts/winners_to_targets.py --cif-dir rlvr_round/dataset --out targets/round1
"""

from __future__ import annotations

from pathlib import Path

import typer
from boltzgen.data import data as D
from boltzgen.data.parse import mmcif

MOLDIR = Path("~/.boltz/mols").expanduser()  # BoltzGen's CCD component cache (NI.pkl etc.)


def _chain_info(chain, design_id: str) -> D.ChainInfo:
    return D.ChainInfo(
        chain_id=int(chain["asym_id"]),
        chain_name=str(chain["name"]),
        mol_type=int(chain["mol_type"]),
        cluster_id=design_id,  # de-novo: every design is its own cluster (no PDB clustering)
        msa_id=-1,             # single-sequence; the trainer skips MSA when use_msa=false
        num_residues=int(chain["res_num"]),
        valid=True,
        entity_id=int(chain["entity_id"]),
    )


def convert_one(cif: Path, struct_dir: Path, record_dir: Path, moldir: Path) -> D.Record:
    """One winner CIF → Structure.npz + Record.json on disk; returns the Record for the manifest."""
    ps = mmcif.parse_mmcif(str(cif), moldir=str(moldir), use_assembly=False)
    design_id = cif.stem
    ps.data.dump(struct_dir / f"{design_id}.npz")
    record = D.Record(
        id=design_id,
        structure=ps.info,
        chains=[_chain_info(c, design_id) for c in ps.data.chains],
        interfaces=[],
        templates=None,
    )
    record.dump(record_dir / f"{design_id}.json")
    return record


def main(
    cif_dir: Path = typer.Option(..., help="dir of selected winner CIFs (rlvr_select's dataset/)"),
    out: Path = typer.Option(..., help="output target_dir (structures/ + records/ + manifest.json)"),
    moldir: Path = typer.Option(MOLDIR, help="CCD mols dir (BoltzGen's ~/.boltz/mols cache)"),
) -> None:
    struct_dir, record_dir = out / "structures", out / "records"
    struct_dir.mkdir(parents=True, exist_ok=True)
    record_dir.mkdir(parents=True, exist_ok=True)

    records, failed = [], []
    for cif in sorted(cif_dir.glob("*.cif")):
        try:
            records.append(convert_one(cif, struct_dir, record_dir, moldir))
        except Exception as e:  # one bad CIF shouldn't sink the batch — report and continue
            failed.append((cif.stem, f"{type(e).__name__}: {e}"))

    D.Manifest(records=records).dump(out / "manifest.json")
    print(f"converted {len(records)} winners → {out} (structures/ records/ manifest.json)")
    for cid, err in failed:
        print(f"  ⚠ skipped {cid}: {err}")
    if not records:
        print("⚠ no targets written — check --cif-dir and --moldir")


if __name__ == "__main__":
    typer.run(main)
