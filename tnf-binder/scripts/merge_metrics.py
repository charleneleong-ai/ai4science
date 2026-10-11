"""Join BoltzGen folding confidence onto the geometric scores, and gate on both.

Two complementary views of a design, and both are needed:

  geometry (score_designs)  measures the pose BoltzGen DESIGNED -- interface area, how much of
                            it lands on the intended epitope, how much on species-divergent
                            residues.
  confidence (folding step) measures whether the model BELIEVES that pose when it re-predicts
                            the complex from scratch -- iptm, design_to_target_iptm,
                            interaction_pae.

A design can have excellent geometry and still be a non-binder: the first nanobody batch buried
plausible interfaces yet returned iptm ~0.15 with interaction_pae ~28 A. Confidence is the
stricter and more predictive filter, so it gates; geometry then orders within the survivors.

Gate on design_to_target_iptm, NOT the global iptm. Global iptm averages over every inter-chain
pair, so a multi-chain target inflates it with its own well-predicted internal interfaces: moving
from the TNF dimer to the trimer lifted global iptm 0.20 -> 0.52 purely by adding three
well-predicted TNF-TNF pairs, while design_to_target_iptm stayed flat at 0.19 and interaction_pae
did not improve. Gating on global iptm would have passed a design whose binder does not dock.

Thresholds follow de novo binder-design practice (BindCraft-style): >= 0.5 for the binder-to-target
interface and interaction pae <= 10. They are explicit and overridable rather than tuned to
whatever this run happened to produce.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from tnf_binder.tables import number, read_rows, write_rows

app = typer.Typer(add_completion=False)
console = Console()

COLUMNS = ["design", "binder_len", "iptm", "design_to_target_iptm", "interaction_pae",
           "design_ptm", "target_ptm", "interface_area", "epitope_focus", "cross_reactivity",
           "divergent_frac", "clashes", "flag", "sequence"]
NO_PAE = 99.0


def merged_row(scored: dict[str, str], folded: dict[str, str], min_d2t: float,
               max_pae: float) -> dict[str, object]:
    d2t = number(folded, "design_to_target_iptm")
    pae = number(folded, "interaction_pae", NO_PAE)
    flags = [f for f in (scored.get("flag") or "").split(",") if f]
    if not folded:
        flags.append("no_fold_metrics")
    else:
        if d2t < min_d2t:
            flags.append(f"low_d2t_iptm({d2t:.2f})")
        if pae > max_pae:
            flags.append(f"high_pae({pae:.1f})")
    return scored | {
        "iptm": number(folded, "iptm"), "design_to_target_iptm": d2t, "interaction_pae": pae,
        "target_ptm": number(folded, "target_ptm"), "design_ptm": number(folded, "design_ptm"),
        "flag": ",".join(flags),
    }


@app.command()
def main(
    scored_csv: Path = typer.Argument(..., help="geometry table from score_designs.py"),
    fold_csv: Path = typer.Argument(..., help="confidence table from read_fold_metrics.py"),
    out: Path = typer.Option(..., "--out", help="the joined, gated table"),
    min_d2t_iptm: float = typer.Option(0.5, "--min-d2t-iptm",
                                       help="threshold on design_to_target_iptm, the "
                                            "binder-specific interface score"),
    max_pae: float = typer.Option(10.0, "--max-pae"),
    relax: bool = typer.Option(False, "--relax",
                               help="keep designs failing the confidence gate, flagged, "
                                    "for inspection"),
) -> None:
    """Join the two tables and keep only what clears both gates."""
    folded = {r["design"]: r for r in read_rows(fold_csv)}
    merged = [merged_row(r, folded.get(r["design"], {}), min_d2t_iptm, max_pae)
              for r in read_rows(scored_csv)]

    passing = [r for r in merged if not r["flag"]]
    kept = merged if relax else passing
    kept.sort(key=lambda r: (-number(r, "design_to_target_iptm"),
                             number(r, "interaction_pae", NO_PAE)))
    write_rows(out, kept, COLUMNS)

    console.print(f"{len(merged)} scored | [bold]{len(passing)}[/] pass both gates "
                  f"(design_to_target_iptm>={min_d2t_iptm}, pae<={max_pae}) | "
                  f"wrote {len(kept)} -> [cyan]{out}[/]")
    if not passing and merged:
        best = max(merged, key=lambda r: number(r, "design_to_target_iptm"))
        console.print(f"  best by design_to_target_iptm: {best['design']} "
                      f"d2t={best['design_to_target_iptm']:.3f} "
                      f"(global iptm={best['iptm']:.3f}) "
                      f"pae={best['interaction_pae']:.1f} -- nothing clears the bar")


if __name__ == "__main__":
    app()
