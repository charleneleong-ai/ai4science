# Official competition page: Description

Captured 2026-10-31 via the Kaggle MCP.

---

The gap between theoretical quantum algorithms and practical quantum advantage often comes down to circuit size. Many high-impact algorithms rely on complex subroutines that require far more gates and qubits than near-term hardware can support. Shrinking these bottleneck subroutines directly reduces the hardware needed to run the full algorithm.

Historically, domain experts have optimized these circuits by hand, publishing incremental improvements one paper at a time. While verifying general quantum circuits is intractable on classical computers, many of the costliest fault-tolerant subroutines can be written as **kickmix** circuits—built from reversible classical logic, diagonal phase gates, and measurement-based uncomputation. As detailed in our companion paper, [*Leveraging verifiable structure in quantum circuits to find more efficient compilations with agentic AI*](https://arxiv.org/abs/2610.23832), our new open-source [kickmix](https://github.com/quantumlib/Qualtran/tree/main/kickmix) library in [Qualtran](https://github.com/quantumlib/Qualtran) can build and fuzz-test circuits with tens of millions of gates on classical inputs in seconds.

Rather than training on a traditional dataset, you will use AI coding agents, classical optimization, and algorithmic insights in a closed verification loop to find smaller, faster circuit compilations—including smart approximations within each benchmark’s accuracy bounds. No quantum physics background or quantum hardware is required. Your submissions will target three real-world bottlenecks:

- **Fusion Energy (1/√x Phasing):** Compressing fixed-point inverse-square-root circuits for pairwise Coulomb interactions used to model alpha-particle stopping power in warm dense plasmas during inertial confinement fusion ([PNAS 2024](https://www.pnas.org/doi/10.1073/pnas.2317772121)).
- **Combinatorial Optimization (Decoded Quantum Interferometry):** Optimizing reversible Reed–Solomon syndrome decoders for Decoded Quantum Interferometry (DQI) on the Optimal Polynomial Intersection problem over finite fields ([Nature 2025](https://www.nature.com/articles/s41586-025-09527-5), [arXiv:2510.10967](https://arxiv.org/abs/2510.10967)).
- **Quantum Chemistry & Clean Energy (Hamiltonian Block Encoding):** Co-optimizing classical tensor factorizations, Givens rotation trees, and quantum lookup tables (QROMs) to simulate FeMoCo—the biological catalyst behind nitrogen fixation, where industrial fertilizer production consumes ~2% of global energy ([PNAS 2017](https://www.pnas.org/doi/abs/10.1073/pnas.1619152114), [PRX 2025](https://journals.aps.org/prx/abstract/10.1103/pb2g-j9cw)).

Could we get more out of a quantum computer simply by asking it to do less? We invite you to help us answer that question.

  
