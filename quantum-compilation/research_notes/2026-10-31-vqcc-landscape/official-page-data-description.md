# Official competition page: data-description

Captured 2026-10-31 via the Kaggle MCP.

---

Your task is to optimize three quantum circuits relevant for combinatorial optimization, fusion research, and chemistry. The more optimal the circuit, the sooner it will be possible to build working implementations on real quantum hardware. 

## Problem specification

The verifiers in **verifiers/** are the authoritative specification, and the paper explains each contract. A verifier simulates your circuit with kickmix on many inputs and checks every shot: inputs restored, scratch qubits back to \\( |0\rangle \\), and the right phase. Measurement-based uncomputation is allowed. Tests are seeded from your submission's hash plus a secret salt, so any change to a circuit gets fresh tests (paper §II.C).

### invsqrt (the paper's Benchmark I, §III.A)

Map \\( |x\rangle|0\rangle \mapsto e^{2\pi i/\sqrt{x}}|x\rangle|0\rangle \\) for every \\( x \\) from 1 to 3,133,452, held in a 22-qubit register named `x`. Each input is tested 64 times, and every phase must be within \\( 2^{-24}/\sqrt{x} \\) turns. Cost: qubits × Toffoli-equivalents (Toffolis plus T gates / 4, with rotations priced in T gates).

### dqi (the paper's Benchmark II, §III.B)

A reversible decoder for the Reed–Solomon code \\( \mathrm{RS}(4095, 4025) \\) over \\( \mathrm{GF}(2^{12}) \\). Given error syndromes in registers `S0` to `S69` and a classical record \\( M \\) in `E0` to `E4094`, it restores the syndromes and applies the phase \\( (-1)^{\sum_j \mathrm{popcount}(e_j \wedge M_j)} \\) for the error \\( e \\). It is tested on 102,400 shots; failures are allowed but raise the cost. Cost: \\( R(\epsilon_+) \times Q^2 \times T \\), with \\( R \\) a repetition factor that grows with the failure rate (paper Eq. (43)).

### dfthcsos (the paper's Benchmark III, §III.C and App. A)

Circuits that load the Givens-rotation angles of a factorized FeMoCo Hamiltonian, one batch of \\( J \\) angles per circuit, plus the factorization and rotation tree in `submission.npz`. A submission is rejected if any circuit uses more than 128 qubits, if the encoded Hamiltonian is more than 20 mHa from the reference, or if the LCU normalization \\( \Lambda \le 49.0 \\) or \\( E_{\mathrm{gap}} \le 3.914898 \\) Ha. Cost: \\( \lambda_{\mathrm{eff}} \times T_{\mathrm{total}} \times Q_{\mathrm{total}} \\) over the whole block encoding (paper Eq. (44)). The full contract is in **verifiers/dfthcsos/SCORING_README.md**, and the [starter kit](https://www.kaggle.com/datasets/sohier/verifiable-circuit-compilation-starter-kit) builds a complete submission from the provided baseline.

## Files

**verifiable_circuit_compilation.pdf**

The paper behind the competition, also on arXiv as [arXiv:2610.23832](https://arxiv.org/abs/2610.23832).

**sample_submission/**

The reference submission, which scores 1.0: `invsqrt.kmb`, `dqi.kmb` and `dfthcsos/` (`angle_k_0.kmb` to `angle_k_9.kmb` and `submission.npz`). Zip its contents to submit; text `.kmx` circuits are accepted too.

**verifiers/**

Each task's verifier, as the leaderboard runs it, with the script that writes its reference circuit. **verifiers/dfthcsos/** also holds that task's inputs.

**kickmix/**

The kickmix circuit simulator, `kickmix.so`, built for the Kaggle notebook image. Add this folder to `sys.path` to `import kickmix`.

The dfthcsos starter kit is a separate dataset, the [Verifiable Circuit Compilation Starter Kit](https://www.kaggle.com/datasets/sohier/verifiable-circuit-compilation-starter-kit). It builds a complete dfthcsos submission from the provided baseline in three steps, and its page describes its files.

## Pending Kickmix Update
Note that the kickmix library should be available on PyPi in the near future. When it's available there we'll drop it from the competition dataset and migrate the metric to use the PyPi copy.
