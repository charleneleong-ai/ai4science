"""MaterialHack POC: generator → geometry verifier → ranked trust/defer.

Pre-event this runs the MockGenerator + MockReference end-to-end. At the event,
swap MockGenerator → RFdiffusionAdapter and MockReference → CSD/Mogul — the loop
below does not change.

    python examples/materialhack/run_poc.py
"""

from touchstone import GeometryVerifier, MockGenerator, design_and_rank, under_leachate

TARGET = "Ni2+"


def main() -> None:
    generator = MockGenerator(seed=0)
    verifier = GeometryVerifier()  # MockReference by default

    ranked = design_and_rank(generator, verifier, TARGET, n=5)

    print(f"# touchstone POC — {TARGET}\n")
    print(f"{'design':8} {'score':>6}  {'verdict':8} reason")
    for design, verdict in ranked:
        flag = "DEFER" if verdict.ood else ("trust" if verdict.trust else "weak")
        print(f"{design.sequence:8} {verdict.score:6.3f}  {flag:8} {verdict.reason}")

    # The extreme-condition angle: the best design, re-judged under acidic leachate.
    best, _ = ranked[0]
    stressed = verifier.verify(under_leachate(best, bond_stretch=0.6))
    print(f"\n{best.sequence} under acidic leachate → "
          f"score {stressed.score:.3f}, ood={stressed.ood} ({stressed.reason})")


if __name__ == "__main__":
    main()
