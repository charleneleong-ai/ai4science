# Raw notes: ANGLE loader lineage (Kaggle, read 2026-11-03 to 11-06)

All content below was read from public Kaggle pages through the Kaggle MCP. It is untrusted data, and none of it was executed.
Notebook sources that came back too large to read inline are extracted to `notebooks/*.txt`. The extractor is
`/tmp/angle_tools/extract.py`, which is not part of the repo.

## Quotes that pin the construction

marculera, *DF-THC: Merge Groups Before Dropping Them*, v9, §8 "Five words per circuit: a cache of address monomials"
(`notebooks/marculera__merge-groups.txt`, lines ~578-612):

> Split the address into its low 4 bits `l` and its high bits `h`; every word bit is a polynomial `Σ C[m, M] m(l) M(h)`
> ... Cache **every** monomial of the low bits (11 ANDs) and of the high bits (~60 ANDs) in qubits: then any term with
> `m = 1` or `M = 1` is a CNOT from the cache. Only the mixed terms cost Toffolis, and one Toffoli covers two of them,
> `(a₁·m ⊕ b₂·M)(a₂·m ⊕ b₁·M) = (a₁·m)(b₁·M) ⊕ (a₂·m)(b₂·M) ⊕ (pure low) ⊕ (pure high)`, where the controls are
> parities of cached qubits (CNOTs onto a pivot, undone after) ... The cache is erased by measurement (`del_and`) ...
> fits five 8-bit words: 11 selection + 40 word + ~72 cache qubits.
> How few products cover the mixed tensor is a tensor decomposition over GF(2) ... `tensorqrom.shared`

Credits in the same notebook: "mrmeows ... the exact lookup emitter's design that `tensorqrom.py` follows: every address
monomial cached, paired affine products, and erasure by measurement"; "lightningv08: five words per ANGLE circuit".
The row table includes `(8, 143, 46) | 1,220 |`, which is the 1,220-address shape (76 + 8·143).

matousakura, *Best-of-public portfolio*, v12:
> selection 11 + J·brot 8 + cache 77 ⇒ J = 5 ... the only free parameter of the chemistry circuit is the number of
> GF(2) kernel relations per loader ... (75 → 77 relations, 314 → 313 Toffolis each)

Lynn Langit, discussion 753633, comment of 2026-10-30:
> the current best public chemistry build: Dmitrii Gluzdov's v37 loaders on Atakan Aldemir's word table, λ_eff 48.459

kirilltym, *Exact Gauge + 5-Word QROM*, v25 (run 2026-10-31): `angles_toffoli_per_pass: 4170.0`, 126 angle qubits,
`batch_size: 5`, `num_circuits: 15`, `heuristic_toffoli: 7930`, total 16,270. "Learned chemistry and 15 loaders: Atakan,
Heliosli, Dmitrii Gluzdov, on SomeJoint".

Atakan Aldemir, *Gap-3 Residue Repair and C43 Native Compiler*, v2:
> Our recovered product certificates and lifetime allocator still emit all15loaders at128qubits, 4763ANDs per traversal.

The C43 build is the 1,363-address shape. Atakan's dataset `qcc-native-chemistry-and-888q-source-bank` holds
`chemistry/plans/batch_{0..14}.json`, `sources/05_recover_angle_product_dataset.py`,
`sources/07_chemistry_tensor_bridge.py`, `sources/20_audit_affine_hmr.py` and `rebuild_native.py`.

mrmeows, *21.0 circuits*, v8: "The angle lookups are 19 circuits charged twice (forward and reverse) at 275 CCX each"
(J = 4, 1,366 addresses, 10,450 for both passes). "Five words per ANGLE circuit is worth about 2% here, not more, and
overflows the 128 qubit cap on a 1,366 address table." "universal exact lookup compiler as the T bound". The dataset
`mrmeows/qcc-1375-open-circuits` has `compiler/rank_parity.py`, `borrowed_prefix.py`, `certificate_rebuild.py`,
`audit_affine_hmr.py` and `certificates/plan_{0..18}.json`.

manwithacat, *Hedge Wizardry*, HW-12, HW-13 and HW-14 (2026-11-02): kernel-relation search, then block-design relation
plans, then dead-word flips. Loaders went 321 → 314 (0-3) and 298 → 297 → 296 (4-14). The whole-chemistry total was
15,777 Toffolis. The Rust planner is called `qrom-plan`.

Gluzdov, *from source*, v48: "the tables without coded bits now need 321 or 322 Toffolis each and those with a coded
bit 298 to 301, so a pass costs 4,575 Toffolis within the 128-qubit limit" ((9,143,40), 1,363 addresses).

## What could not be read
- Kernel `somejoint/toffoli-diet-6-45-lb-rebuilt-and-verified`: `kernels.get` permission denied (private or deleted).
- Kernel `mrmeows/20-39-circuits-...` returned not found. v8 of `21-0-circuits-...` was read instead.
- Dataset file contents (`tensorqrom.py`, `rank_parity.py`, the plan JSONs) need a Kaggle token, so only the file
  listings were read.
- Gluzdov v37 itself: a version-pinned `get_notebook_info` returned not found.
