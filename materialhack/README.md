# MaterialHack — AI track application

Application notes + the project plan for [MaterialHack](https://www.materialhack.co.uk/). Submitted June 2026.

## The event
- 48-hour sprint, **London, Jun 26–28 2026**. Applications closed **Jun 15**.
- By **Nucleate UK**, powered by **ARIA** (Advanced Research + Invention Agency); "Manufacturing Abundance Opportunity Space."
- Apply as individuals; teams form at the Friday social.

## AI track
> Design protein–metal interactions and hybrid bio-metal systems for **critical-mineral recovery** and advanced metal-based materials — the intersection of biology and materials science, using AI tools and provided datasets.

Three guiding questions:
1. Can we enhance metal-binding properties of proteins by **de novo design**?
2. How can we accurately **predict protein behaviour in extreme environments**?
3. What **metal-relevant industrial processes** would benefit from proteins?

**The angle:** Q2 is the verification question in disguise. Real critical-mineral recovery (bioleaching ores, e-waste, brines, mine tailings) happens in **acidic (pH 1–2), hot, high-salinity leachates** — exactly where a model trained on benign-condition data goes **out-of-distribution and confidently wrong**. The contrarian move (most teams just maximise predicted binding affinity) is to design for **robustness + know when you can't be trusted**. That threads all three questions: de-novo design (Q1) of binders that survive extreme leachate (Q2) for metal recovery (Q3).

## Submitted answer

> *"What is the most unconventional idea you've ever had for solving a scientific or engineering problem, and how would you go about testing it?"* (≤300 chars)

**Submitted (254 chars) — extreme-environment / OOD:**
> De novo design metal-binders for the conditions that break them — hot, acidic, saline leachate. Treat "extreme environment" as an out-of-distribution problem: the model must flag when it can't be trusted. Test vs held-out extremophile metal-binding data.

**Alternative (282 chars) — selectivity-from-failures:**
> Design metal-binding proteins by training on failures, not successes — non-binders map selectivity best, and selectivity is the whole game in critical-mineral recovery. Score by an independent physics oracle (Boltz-2 co-fold) + held-out real metal-site data, not the model's own guess.

## Doing it without a wet lab (48h, in-silico)

The verifier doesn't need to be wet — it just needs to be **independent of the generator**. Two hackathon-doable sources of ground truth:

1. **A physics oracle.** Co-fold each designed protein *with the metal ion* (AlphaFold3 / **Boltz-2** — Boltz-2 also predicts affinity) and read confidence (ipTM / PAE); estimate coordination-site binding energy with a quick QM/MM or MD pass. Independent of the design model.
2. **The experimental record that already exists.** No new chemistry — decades of it sit in **MetalPDB / BioLiP** (curated real metal-binding sites + selectivities). Use a **held-out slice as the judge**: retrospective validation, zero pipetting.

### 48h loop
```
generate            verify (independent)              judge
RFdiffusionAA   →   Boltz-2 co-fold (ipTM/PAE,    →   held-out MetalPDB/BioLiP hits
ProteinMPNN         affinity) + QM/MM site energy     + OOD flag on extreme-condition inputs
```
Negatives (non-binders / decoys) come from existing data + physics scores — so the "train on failures for selectivity" idea stays fully in-silico.

## Pre-sprint prep
- [ ] Skim **Boltz-2** affinity/confidence outputs + API (don't learn it during the sprint).
- [ ] Pull a **MetalPDB / BioLiP** slice for the target metal(s); define the held-out split up front.
- [ ] Decide the OOD signal for "extreme environment" inputs (energy-based score / ensemble disagreement).
- [ ] RFdiffusionAA / ProteinMPNN generation pipeline ready to run.

## Sources
- Event: [materialhack.co.uk](https://www.materialhack.co.uk/) · Nucleate UK · ARIA.
- Tooling (illustrative, ground versions before use): AlphaFold3 / Boltz-2 (co-fold + affinity), RFdiffusionAA & ProteinMPNN (Baker lab, design), MetalPDB / BioLiP (curated metal-binding data).

## Learnings
- **"Verifier" ≠ "wet lab."** My first instinct framed the test around a wet-lab assay — impossible in a 48h in-silico hackathon. The fix: a verifier only has to be *independent of the generator*, and that can be a physics oracle (Boltz-2 co-fold) or the existing experimental record (MetalPDB/BioLiP), no pipetting. Reframing "reality" as *any independent ground truth* is the reusable move.
- **Read the track's own questions before pitching.** The three guiding questions made Q2 ("behaviour in extreme environments") the unlock — extreme leachate = OOD = my wheelhouse — which is sharper than the generic "design a binder" angle most teams take. The brief tells you the angle; use it.
- **Even the verifier needs verifying.** Grounding Boltz-2 ([log entry](../learning-log.md)) surfaced a 2026 reliability study: its affinities aren't trustworthy *quantitatively*. So the plan shifted to using it as a *ranking* signal with held-out data + physics as meta-verifier — a correction recall alone would have missed.
- **Char limits force the idea to its core.** Compressing to 300 chars killed the hedging and exposed the single claim (selectivity-from-failures / OOD-flag). A good idea survives the squeeze; a vague one doesn't.
