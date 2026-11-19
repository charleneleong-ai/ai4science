# fork492m

[fork492](../fork492) (public 0.492: fork487 plus the test-set network rerank at `LAM` 0.35) plus one cell: MIST
`msg_aug` on ranks 2–25. Run as the private kernel `charyeezy/casmi-fork492m`.

The new cell runs after the APEX library overlay, so after N-1c.
- **Fingerprint:** it predicts a 4,096-bit Morgan fingerprint from each molecule's positive-mode test spectra. MIST
  has no negative ions, so molecules with only negative spectra are skipped. The cell averages the predictions over
  the molecule's spectra.
- **Reorder:** it reorders ranks 2–25 by 0.5·z(−log rank) + 0.5·z(cosine to the prediction). Rank 1 never moves,
  and the cell asserts it. On any error the submission is unchanged.
- **Preprocessing:** MIST's own, as in the 2026-11-01 probe
  ([notes](../../research_notes/2026-11-01-mist-msg-aug/notes.md)). Peaks above precursor + 1 are dropped, the rest
  sqrt-normalised and cut to the top 60. Subformulae are assigned at 20 ppm against the rank-1 candidate's formula,
  since the test gives none, with the spectrum's adduct. The top 10 assigned peaks plus CLS go in, with the
  instrument as Q-ToF.
- **Dataset:** the private dataset `charyeezy/casmi-mist-msg-aug` holds the `mist` package (coleygroup/MIST
  `48e339b`, unmodified, MIT) and the Zenodo 23303289 checkpoint (md5 `3e92179240e3043733627d720f4b7b86`).

Checks before submitting, all on the A100:
- **Torch 2.10, no `h5py`:** the cell reproduced the probe's FRIGID-50 Tanimoto, 0.385.
- **fork492 v1 smoke output:** it scored and reordered 11 lists, with 0 nulls, rank 1 unchanged and the candidate
  sets unchanged; 27 of 400 molecules have only negative spectra.
- **Reach,** replayed on the 953 natural-product holdout lists, with the shipped scorer's order standing in for the
  pipeline:
  - 813 lists reordered;
  - the true structure moved in 266 lists, 228 up and 38 down;
  - MRR@25 0.576 → 0.618, +0.042 (SE 0.003): Enveda-180 +0.022, natural products +0.105.

  That holdout does not predict the leaderboard, so this shows reach, not a gain. The natural-product figure may also
  be inflated if those structures are in MassSpecGym.

The gate is `~/casmi-data/gates/gate_fork492m.py`: `gate_fork492.py`, then MIST's `[mist] done:` line with rank 1
unchanged, nulls at most 5% of molecules, and at least one list changed, and no `[mist] FAILED`.

| kernel version | change | public MRR@25 |
|---|---|---|
| 1 | fork492 + MIST `msg_aug` on ranks 2–25 (weight 0.5) | running |
