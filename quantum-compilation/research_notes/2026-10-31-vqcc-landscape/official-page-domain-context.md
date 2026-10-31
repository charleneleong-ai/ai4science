# Official competition page: Domain Context

Captured 2026-10-31 via the Kaggle MCP.

---

You don't need a background in quantum physics to participate! This section gives you the essentials: what quantum circuits are, what makes them expensive, why the circuits in this competition can be checked automatically, and some more details behind each of the three tasks.

### Quantum circuits in 60 seconds

A classical bit is either 0 or 1. A **[qubit](https://en.wikipedia.org/wiki/Qubit)**, like a classical bit, can be 0 or 1—or in a *[superposition](https://en.wikipedia.org/wiki/Quantum_superposition)* of both. n qubits can be in a superposition of all 2ⁿ bit strings at once, each with its own complex *[amplitude](https://en.wikipedia.org/wiki/Probability_amplitude)*. Quantum algorithms work by steering these amplitudes so that wrong answers cancel out and right answers reinforce each other. This effect is called *[interference](https://en.wikipedia.org/wiki/Wave_interference#Quantum_interference)*.

Quantum programs are specified with **[circuits](https://en.wikipedia.org/wiki/Quantum_circuit)**: a sequence of **[gates](https://en.wikipedia.org/wiki/Quantum_logic_gate)**, usually acting on one to three qubits. Many important quantum gates are reversible versions of ordinary logic, such as NOT, [controlled-NOT](https://en.wikipedia.org/wiki/Controlled_NOT_gate), and the **[Toffoli gate](https://en.wikipedia.org/wiki/Toffoli_gate)**, which flips a target bit when two control bits are both 1. In many quantum algorithms, the subroutines that dominate the cost are essentially classical arithmetic (adders, multipliers, table lookups) applied to a superposition of many inputs at once.

There's a catch: quantum computation must be **[reversible](https://en.wikipedia.org/wiki/Reversible_computing)**. You can't erase or overwrite a value, so temporary results go into extra "scratch" qubits (*[ancillas](https://en.wikipedia.org/wiki/Ancilla_bit)*), and every scratch qubit must be returned to 0 at the end. This clean-up is sometimes called *[uncomputation](https://en.wikipedia.org/wiki/Uncomputation)*. Leftover garbage would destroy the interference that the algorithm relies on.

### What makes a circuit expensive

Physical qubits are noisy, so large quantum algorithms will run on **[error-corrected (logical) qubits](https://en.wikipedia.org/wiki/Physical_and_logical_qubits)**, each built from some number of physical qubits (depending on the [error correcting code](https://en.wikipedia.org/wiki/Quantum_error_correction) used). On such machines most gates are cheap, but **Toffoli gates** (and similar "non-[Clifford](https://en.wikipedia.org/wiki/Clifford_gate)" gates such as T gates) are expensive: each one consumes a special resource state that dedicated "[magic state factories](https://en.wikipedia.org/wiki/Magic_state_distillation)" take a long time to produce. That's why the cost of a fault-tolerant circuit is usually summarized by two numbers:

- the **Toffoli count**, which drives the runtime, and
- the **logical qubit count**, which drives the size of the machine.

Compiling an algorithm into a concrete circuit means navigating a huge design space: which arithmetic decompositions to use, when to allocate and clean up scratch qubits, and how much precision each step really needs. One classic trick is **measurement-based uncomputation**. Instead of spending a second Toffoli to clean up a scratch qubit that holds *a* AND *b*, you measure it and cancel the leftover sign with a cheap gate, halving the cost. Experts have historically discovered tricks like this by hand, one paper at a time.

### Why these circuits can be verified automatically

For general quantum circuits, checking that a candidate does what it should is believed to be intractable: beyond roughly 40 qubits, it's out of reach for classical computers. The subroutines in this competition are different. They're built only from **reversible classical gates, phase gates, and measurement-based uncomputation**, a family called *kickmix circuits*. Each of these instructions maps a definite bit string to another definite bit string, at most changing its phase. So a classical computer can run the circuit one input at a time, tracking just one bit per qubit plus a phase. The open-source [kickmix](https://github.com/quantumlib/Qualtran/tree/main/kickmix) simulator in Google's Qualtran library processes 64 to 512 test inputs with each CPU instruction and handles circuits with more than 100 million gates on a single core.

This is enough to check the quantum behavior, too. Quantum mechanics is linear, so a circuit that maps every bit-string input to the right output with the right phase, and leaves its scratch qubits clean, also behaves correctly on every superposition of those inputs. Nearly correct circuits are nearly right: if a circuit fails on a fraction ε of inputs (weighted the way the real algorithm weights them), its quantum output is within a distance of 2√ε of the ideal output.

So the verifier **fuzz tests** your circuit. Depending on the input domain size, it simulates the circuit across every valid input or across many random inputs drawn the way the full algorithm weights them, checking that the inputs are restored, the scratch qubits are clean, and the phase is right. Random inputs and internal measurement outcomes are seeded with a salted cryptographic hash of your submission (the *Fiat–Shamir heuristic*), so there's no hidden test set and any change to your submission gets fresh random draws. Passing many random tests is strong evidence: for example, a circuit that passes 100,000 of them without a single failure fails on fewer than 0.1% of inputs, with overwhelming confidence.

Because verification is fast and trustworthy, it gives AI coding agents, classical optimization, and automated search an instant reward signal. When the [ecdsa.fail](https://ecdsa.fail) challenge ran a similar verifiable competition for cryptography circuits in 2026, more than 100 participants, most with no quantum computing background, used AI coding agents to cut the Toffoli count roughly in half.

Not every part of a quantum algorithm fits this mold. Steps that create genuine superpositions, such as preparing arbitrary quantum states, can't be checked by simulating one input at a time. Each task therefore isolates the kickmix subroutine that dominates its algorithm's cost, and prices any other parts with analytical formulas, so your score still reflects the full algorithm.

### The three tasks at a glance

| Task | Scientific goal | Quantum algorithm | Bottleneck you'll compress |
|---|---|---|---|
| Fusion Energy | Predict how energy is deposited in fusion fuel | Step-by-step ([Trotter](https://en.wikipedia.org/wiki/Hamiltonian_simulation#Product_formulas)) simulation of electrons on a grid | Phasing by 1/√x for every pair of electrons |
| Combinatorial Optimization | A verifiable quantum advantage for optimal polynomial intersection problem | Decoded Quantum Interferometry (DQI) | A reversible Reed–Solomon syndrome decoder |
| Quantum Chemistry & Clean Energy | Understand how nature [fixes nitrogen](https://en.wikipedia.org/wiki/Nitrogen_fixation) (the FeMoCo cluster) | Phase estimation on an encoded molecular Hamiltonian | Tensor factorizations, Givens trees, and lookup tables that load the Hamiltonian |

### Fusion Energy: stopping power in fusion fuel

In [inertial confinement fusion](https://en.wikipedia.org/wiki/Inertial_confinement_fusion), lasers crush a tiny capsule of fuel to extreme temperatures and densities. Whether the fuel ignites depends on how quickly fast charged particles, such as the [alpha particles](https://en.wikipedia.org/wiki/Alpha_particle) released by fusion reactions, deposit their energy in the surrounding fuel. This rate is called the **[stopping power](https://en.wikipedia.org/wiki/Stopping_power_%28particle_radiation%29)**. On the way to ignition, the fuel passes through *[warm dense matter](https://en.wikipedia.org/wiki/Warm_dense_matter)*, a state too hot and dense for standard theories and hard to probe experimentally, which makes stopping power a key uncertainty in target design ([Rubin et al., PNAS 2024](https://www.pnas.org/doi/10.1073/pnas.2317772121)).

A quantum computer can simulate the electrons directly. It stores every electron's position on a 3D grid and steps time forward, alternating between kinetic-energy updates and potential-energy updates. The potential energy is the [Coulomb repulsion](https://en.wikipedia.org/wiki/Coulomb%27s_law) between every pair of electrons, which is proportional to 1/distance. For each pair, the circuit computes the squared distance x, phases by **1/√x** to about 7 significant digits (a relative error of 2⁻²⁴), and cleans up any scratch work. With 218 electrons, there are 23,653 pairs per time step, so this one subroutine dominates the cost of the whole simulation. The task isolates it: given x, apply the right phase for every possible x.

You may recognize the classical version of this problem. The famous "fast inverse square root" from the video game Quake III combines a clever first guess with one step of Newton's method. The quantum version must be reversible, use fixed-point arithmetic, and be far more accurate. The reference circuit follows the same recipe; you're free to find a better one.

### Combinatorial Optimization: Decoded Quantum Interferometry

Suppose that at every input in a finite field, you are given a set of allowed output values, and your goal is to find a low-degree polynomial whose curve lands inside as many of those allowed sets as possible. This problem is called **Optimal Polynomial Intersection (OPI)** (closely related to *list recovery* in coding theory). Checking how many sets a proposed polynomial hits takes milliseconds on an ordinary computer, yet when the allowed sets are large—about half the field at every point in our benchmark—finding an optimal or near-optimal polynomial appears beyond the reach of known classical algorithms. That combination makes OPI a leading candidate for a verifiable [quantum advantage](https://en.wikipedia.org/wiki/Quantum_supremacy) ([Jordan et al., Nature 2025](https://www.nature.com/articles/s41586-025-09527-5)).

**Decoded Quantum Interferometry (DQI)** tackles OPI using the [quantum Fourier transform](https://en.wikipedia.org/wiki/Quantum_Fourier_transform) to mathematically transform this hard optimization problem into a classic error-correction problem: decoding a **[Reed–Solomon code](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction)**—the same family of codes used to read scratched QR codes, damaged CDs, and noisy deep-space transmissions. Classical computers decode a single Reed–Solomon codeword quickly using tools like the extended Euclidean algorithm, and [Khattar et al. (arXiv:2510.10967)](https://arxiv.org/abs/2510.10967) estimate that DQI could solve classically intractable OPI instances with just a few million Toffoli gates—with the Reed–Solomon decoder accounting for over 90% of the total circuit cost.

The catch is that in DQI, this decoder must run *inside* the quantum computer: reversibly, on a massive superposition of error patterns at once, and leaving zero scratch garbage behind so quantum interference can constructively amplify high-scoring polynomials. In optimized DQI, the error pattern is erased early via measurement-based uncomputation, leaving behind a ±1 phase sign. Your circuit's job is to recover that sign from the error *syndromes* (a compact fingerprint of the errors), cancel it, and clean up all workspace qubits. Unusually, your decoder doesn't even have to be perfect: because every candidate polynomial sampled by DQI is checked classically at the end anyway, an occasional decoding failure only adds a small repetition overhead—a trade-off captured directly in the score.

  
### Quantum Chemistry & Clean Energy: simulating nitrogen fixation

Industrial fertilizer production relies on the [Haber–Bosch process](https://en.wikipedia.org/wiki/Haber_process), which turns nitrogen from the air into ammonia at high temperature and pressure. Bacteria do the same thing at room temperature using the enzyme [nitrogenase](https://en.wikipedia.org/wiki/Nitrogenase), whose active site is a cluster of iron, molybdenum and sulfur atoms known as **[FeMoCo](https://en.wikipedia.org/wiki/FeMoco)**. Its electrons are so strongly correlated that classical computers can't reliably predict its chemistry. A quantum computer could compute the energy of each step of the reaction and reveal how nature does it ([Reiher et al., PNAS 2017](https://www.pnas.org/doi/abs/10.1073/pnas.1619152114)).

The algorithm is *[quantum phase estimation](https://en.wikipedia.org/wiki/Quantum_phase_estimation_algorithm)*, which repeatedly applies a circuit that encodes the molecule's **[Hamiltonian](https://en.wikipedia.org/wiki/Hamiltonian_%28quantum_mechanics%29)** (the operator that gives its energy). The Hamiltonian is specified by an enormous table of numbers: for N [orbitals](https://en.wikipedia.org/wiki/Molecular_orbital), roughly N⁴ of them. The most efficient known algorithms ([Low et al., PRX 2025](https://journals.aps.org/prx/abstract/10.1103/pb2g-j9cw)) compress this table with a *tensor factorization*, stream the compressed data into the quantum computer with **QROM** circuits (quantum read-only memory: lookup tables that work on a superposition of addresses), and apply it as sequences of **Givens rotations** (rotations that mix pairs of orbitals). These lookups and controlled rotations dominate the Toffoli count.

This task couples classical and quantum design. Your choice of factorization affects both how much data the circuits have to load and how many times phase estimation has to repeat them, and it must still reproduce the molecule's energies accurately. The lookup circuits that load the angle data are yours to design, as are the classical Givens tree and bit precision. The state-preparation and Givens-rotation gates themselves, which aren't kickmix circuits, are priced with analytical formulas.

### Key terms

- **[Ancilla (scratch qubit)](https://en.wikipedia.org/wiki/Ancilla_bit):** a temporary qubit that must be returned to 0 by the end of the circuit.
- **Block encoding:** a circuit that embeds a scaled copy of a matrix, such as a Hamiltonian, inside a larger quantum operation, so that algorithms like phase estimation can use it.
- **[Hamiltonian](https://en.wikipedia.org/wiki/Hamiltonian_%28quantum_mechanics%29):** the operator that gives a physical system's energy.
- **[kickmix circuit](https://github.com/quantumlib/Qualtran/tree/main/kickmix):** a circuit made only of reversible classical gates, phase gates, and measurement-based uncomputation, which makes it efficient to verify on a classical computer.
- **[Logical qubit](https://en.wikipedia.org/wiki/Physical_and_logical_qubits):** an error-corrected qubit built from many physical qubits.
- **Measurement-based uncomputation:** cleaning up a scratch qubit by measuring it and cancelling the resulting sign with cheap gates, instead of reversing the computation.
- **QROM:** a quantum lookup table that returns the entry for every address in a superposition at once.
- **Syndrome:** a compact fingerprint of the errors in a codeword, which can be computed without knowing the original message.
- **[Toffoli gate](https://en.wikipedia.org/wiki/Toffoli_gate):** a controlled-controlled-NOT gate; the main expensive operation on error-corrected hardware, and a key cost metric in this competition.
- **[Trotterization](https://en.wikipedia.org/wiki/Hamiltonian_simulation#Product_formulas):** simulating how a system evolves over time by alternating many small steps of simpler evolutions.

### Learn more

- [Khattar et al. (arXiv:2610.23832)](https://arxiv.org/abs/2610.23832), the companion paper for this competition covering the `kickmix` package, verification guarantees, and all three benchmarks.
- [Quantum computing for the very curious](https://quantum.country/qcvc), a free introductory essay by Andy Matuschak and Michael Nielsen, and Wikipedia's [overview of quantum computing](https://en.wikipedia.org/wiki/Quantum_computing).
- [Quirk](https://algassert.com/quirk), a drag-and-drop quantum circuit simulator that runs in your browser.
- The [kickmix documentation](https://github.com/quantumlib/Qualtran/tree/main/docs/kickmix) in Google's open-source Qualtran library: the circuit file format, the instruction set, and the Python API.
- The [Qualtran documentation](https://qualtran.readthedocs.io/en/latest/), for the library kickmix is part of.
- [Babbush et al. (arXiv:2603.28846)](https://arxiv.org/abs/2603.28846), which introduced kickmix circuits and this style of verification to estimate the cost of breaking elliptic-curve cryptography.
- The papers cited in the Description for each task.
