"""Write the provenance record for a submission assembled from the public parts market.

Every circuit in such a submission was published by someone else under Apache-2.0. This
records who wrote each one, where it came from, and what our own verifier measured for it, so
the authors are credited and the submission can be reproduced byte for byte from the sha256
values. Each part is identified by hashing the file that was actually submitted and looking
that hash up in the parts dataset's manifest, so the record cannot name a circuit the
submission did not contain. For the two single-circuit tasks the verifier's own
`candidate_sha256` is cross-checked against the same hash; the chemistry task's
`submission_sha256` is a digest of the parsed submission rather than of any file, so there it
is the manifest's `factorization_sha256` that the npz must match.

    uv run --extra tracking python scripts/build_parts_record.py \\
        <manifest-dir> <rows.jsonl> <submission-root> \\
        --submission-ref 57834335 --submitted-at 2026-10-31T09:27:57Z
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Annotated, Any

import typer
import yaml

from quantum_compilation.cost import TASKS, dfthcsos_constraint_margins

PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS_DIR = PROJECT_ROOT / "docs" / "experiments"

HEADER = """\
# Provenance for a submission assembled from public Apache-2.0 circuits, fetched from a Kaggle
# parts dataset and re-scored locally with the official verifiers before submitting. Generated
# by scripts/build_parts_record.py from that dataset's manifest.csv and our own rows.jsonl;
# do not hand-edit.
#
# We wrote none of these circuits. `author` published the part and `core_author` records the
# lineage behind it. Only `our_measurement` and `constraint_margins` are ours -- every other
# per-part field is copied from the dataset's manifest and is the part author's own report.
#
# In `constraint_margins`, `fraction_of_limit_used` is closeness-to-rejection in both
# directions, so 1.0 means sitting exactly on the limit whether it is a cap or a floor.
"""

# Per part: who made it, where it came from, and what its author reported for it.
CARRY = (
    "author",
    "core_author",
    "source_notebook_url",
    "notebook_title",
    "source_notebook_version",
    "release_date",
    "licence",
    "qubits",
    "change_summary",
    "notes",
)
# The chemistry factorization's shape, which is the part of the record that drives cost.
CHEM_EXTRA = (
    "R", "B", "C", "J", "brot", "KEEP_BITS",
    "Lambda", "E_gap", "lambda_eff", "T_total", "T_angles", "Q_total",
)
# Manifest fields that are numbers rather than prose, so the record stores them typed and its
# readers need no casts.
NUMERIC = {
    "qubits", "R", "B", "C", "J", "brot", "KEEP_BITS",
    "Lambda", "E_gap", "lambda_eff", "T_total", "T_angles", "Q_total",
}


# Which submitted file identifies each part, and the manifest column its hash appears in.
# dfthcsos is keyed by its npz because the verifier's own digest for that task covers the
# parsed submission rather than any single file.
IDENTIFIED_BY = {
    "invsqrt": ("invsqrt.kmb", "sha256"),
    "dqi": ("dqi.kmb", "sha256"),
    "dfthcsos": ("dfthcsos/submission.npz", "factorization_sha256"),
}


def typed(field: str, raw: str) -> Any:
    if field not in NUMERIC:
        return raw.strip()
    text = raw.strip()
    return int(text) if text.isdigit() else float(text)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def part_record(row: dict[str, str], verdict: dict[str, Any], task: str, file_sha: str) -> dict[str, Any]:
    # For the single-circuit tasks the verifier hashed the same bytes we just did, so a
    # mismatch means the record would describe a file the verifier never saw.
    candidate = verdict["candidate_sha256"]
    if task != "dfthcsos" and candidate != file_sha:
        raise ValueError(
            f"{task}: verifier scored {candidate} but the submitted file hashes to {file_sha}"
        )
    record: dict[str, Any] = {"submitted_file_sha256": file_sha}
    fields = CARRY + (CHEM_EXTRA if task == "dfthcsos" else ())
    record |= {f: typed(f, row[f]) for f in fields if row.get(f, "").strip()}
    record["our_measurement"] = {
        "status": verdict["status"],
        "cost": verdict["score"],
        "seconds": round(verdict["seconds"], 1),
        "verifier_candidate_sha256": candidate,
    }
    if task == "dfthcsos":
        record["constraint_margins"] = dfthcsos_constraint_margins(verdict["detail"])
    return record


def main(
    manifest_dir: Annotated[Path, typer.Argument(help="Directory holding the parts manifest.csv.")],
    rows: Annotated[Path, typer.Argument(help="rows.jsonl written by `vqcc verify`.")],
    submission_root: Annotated[
        Path, typer.Argument(help="The submission directory whose files were uploaded.")
    ],
    dataset: Annotated[str, typer.Option(help="Kaggle dataset the parts came from.")] = (
        "manwithacat/vqcc-apache-public-circuits"
    ),
    submission_ref: Annotated[str, typer.Option(help="Kaggle submission id.")] = "",
    submitted_at: Annotated[str, typer.Option(help="UTC submission time, ISO-8601.")] = "",
    out: Annotated[
        Path | None,
        typer.Option(help="Defaults to docs/experiments/<submission-ref>-parts.yaml."),
    ] = None,
) -> None:
    manifest_rows = list(csv.DictReader((manifest_dir / "manifest.csv").open()))
    measured = {}
    for line in rows.read_text().splitlines():
        if line.strip():
            verdict = json.loads(line)
            measured[verdict["task"]] = verdict

    missing = sorted(set(TASKS) - set(measured))
    if missing:
        raise ValueError(f"{rows} has no verdict for {missing}")

    parts = {}
    for task in TASKS:
        relative, column = IDENTIFIED_BY[task]
        file_sha = sha256(submission_root / relative)
        found = [r for r in manifest_rows if r[column] == file_sha]
        if len(found) != 1:
            raise ValueError(
                f"{task}: {relative} hashes to {file_sha}, which matches {len(found)} "
                f"manifest rows by {column} -- expected exactly one"
            )
        parts[task] = part_record(found[0], measured[task], task, file_sha)

    record = {
        "dataset": dataset,
        "submission_ref": submission_ref,
        "submitted": submitted_at,
        "parts": parts,
    }
    # width keeps the long lineage and change-summary strings on one line each; allow_unicode
    # keeps author names spelled the way their authors spell them.
    body = yaml.safe_dump(
        record, sort_keys=False, allow_unicode=True, default_flow_style=False, width=10**9
    )
    target = out or EXPERIMENTS_DIR / f"{submission_ref}-parts.yaml"
    target.write_text(HEADER + body)
    typer.echo(f"wrote {target}")


if __name__ == "__main__":
    typer.run(main)
