"""In-place lookup synthesis by span completion, the strategy decoded from 6c6a734d's invsqrt lookup.

Qubits: the address and one slot per output, no ancillas. Each register holds a Boolean function
of the address, a truth table packed into an int. CX is free, so a CCX control can carry anything
in V_s, the span of 1 and every register except the target slot s.

1. Every slot of an independent output gets f * g, with f and g random members of the span so far
   (new to it), which makes V + V*V the whole function space.
2. Each such slot in turn takes products f_i * g_i (f_i random in V_s, g_i solved by GF(2)
   elimination) until its value lies in V_s + T (T: span of 1, the address and the targets) with a
   part outside V_s, or differs from an unplaced target by a member of V_s + products. CX then
   removes the V_s part, so the slot holds a new direction of T.
3. A CX network maps the slots onto the targets; dependent outputs (XORs of others and affine
   address functions) are CX'd into their own slots for free.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field

import numpy as np

Gate = tuple  # ("X", t) | ("CX", c, t) | ("CCX", a, b, t)
SAMPLE = 6  # waiting slots planned per step; the cheapest is applied
CANDIDATES = 4  # random multipliers tried per product; the one adding most to the span is kept
VS, TARGET = -1, -2  # generator kinds besides products, which carry their multiplier's index


class TaggedSpan:
    """GF(2) span of ints; each basis vector carries a mask of the generators that make it."""

    def __init__(self, vectors: list[int] = ()) -> None:
        self.pivots: dict[int, tuple[int, int]] = {}
        for v in vectors:
            self.add(v)

    def reduce(self, v: int, tag: int = 0) -> tuple[int, int]:
        while v and (p := self.pivots.get(v.bit_length() - 1)):
            v ^= p[0]
            tag ^= p[1]
        return v, tag

    def add(self, v: int, tag: int = 0) -> bool:
        v, tag = self.reduce(v, tag)
        if v:
            self.pivots[v.bit_length() - 1] = (v, tag)
        return bool(v)

    def __contains__(self, v: int) -> bool:
        while v and (p := self.pivots.get(v.bit_length() - 1)):
            v ^= p[0]
        return not v


class SlotSpans:
    """Generators for completing one slot: V_s and products (`without_t`), and the same plus the
    unplaced targets (`with_t`); a tag bit k means generator k."""

    def __init__(self) -> None:
        self.gens: list[tuple[int, int]] = []  # (VS, V_s index) | (TARGET, j) | (f index, V_s index)
        self.with_t, self.without_t = TaggedSpan(), TaggedSpan()

    def add(self, gen: tuple[int, int], v: int, target_only: bool = False) -> None:
        self.gens.append(gen)
        bit = 1 << (len(self.gens) - 1)
        self.with_t.add(v, bit)
        if not target_only:
            self.without_t.add(v, bit)

    def decode(self, tag: int, vs: list[int], products: int) -> tuple[list[int], int]:
        """Split a tag into each product's second factor and the V_s part (targets need no gates)."""
        g_of, v = [0] * products, 0
        for k, (kind, i) in enumerate(self.gens):
            if tag >> k & 1 and kind != TARGET:
                if kind == VS:
                    v ^= vs[i]
                else:
                    g_of[kind] ^= vs[i]
        return g_of, v


def packed(bits: np.ndarray) -> int:
    return sum(1 << int(i) for i in np.flatnonzero(bits))


def address_functions(entries: int) -> list[int]:
    index = np.arange(entries)
    return [packed((index >> i) & 1) for i in range((entries - 1).bit_length())]


def independent_outputs(table: np.ndarray) -> list[int]:
    """Columns of `table` that need Toffolis: the rest are XORs of earlier ones and affine address
    functions, so CX and X gates produce them for free."""
    span = TaggedSpan([(1 << len(table)) - 1, *address_functions(len(table))])
    return [j for j in range(table.shape[1]) if span.add(packed(table[:, j]))]


def apply(regs: list[int], gate: Gate, one: int) -> None:
    match gate:
        case ("X", t):
            regs[t] ^= one
        case ("CX", c, t):
            regs[t] ^= regs[c]
        case ("CCX", a, b, t):
            regs[t] ^= regs[a] & regs[b]


@dataclass
class InPlaceLookup:
    """Synthesise `table` (entries x outputs, 0/1) onto address + one qubit per output; `gates` is
    the circuit and `slot_of` maps output -> qubit."""

    table: np.ndarray
    seed: int = 0
    gates: list[Gate] = field(default_factory=list, init=False)
    slot_of: dict[int, int] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        entries, self.m = self.table.shape
        self.one = (1 << entries) - 1
        self.address = address_functions(entries)
        self.n = len(self.address)
        self.regs = [*self.address, *([0] * self.m)]
        self.columns = [packed(self.table[:, j]) for j in range(self.m)]
        self.independent = independent_outputs(self.table)
        self.targets = [self.columns[j] for j in self.independent]
        self.rng = random.Random(self.seed)

    @property
    def toffolis(self) -> int:
        return sum(g[0] == "CCX" for g in self.gates)

    def gate(self, *g: object) -> None:
        apply(self.regs, g, self.one)
        self.gates.append(g)

    def express(self, v: int, exclude: int | None = None) -> tuple[bool, list[int]]:
        """v as constant ^ XOR of registers (other than `exclude`)."""
        others = [r for r in range(len(self.regs)) if r != exclude]
        span = TaggedSpan()
        span.add(self.one, 1)
        for k, r in enumerate(others):
            span.add(self.regs[r], 2 << k)
        rest, tag = span.reduce(v)
        assert rest == 0, "not in the span"
        return bool(tag & 1), [r for k, r in enumerate(others) if tag >> (k + 1) & 1]

    def load(self, v: int, target: int, busy: int | None = None) -> tuple[int, list[Gate]]:
        """Make a register (not `target`, not `busy`) hold v; return it and the gates that undo it."""
        flip, rs = self.express(v, exclude=target)
        into = next(r for r in rs if r != busy)
        undo = [("CX", r, into) for r in rs if r != into] + [("X", into)] * flip
        for g in undo:
            self.gate(*g)
        return into, undo

    def product_into(self, f: int, g: int, t: int) -> None:
        a, undo_f = self.load(f, t)
        b, undo_g = self.load(g, t, busy=a)
        self.gate("CCX", a, b, t)
        for undo in reversed(undo_f + undo_g):
            self.gate(*undo)

    def xor_into(self, v: int, t: int) -> None:
        flip, rs = self.express(v, exclude=t)
        for r in rs:
            self.gate("CX", r, t)
        if flip:
            self.gate("X", t)

    def random_in(self, vectors: list[int]) -> int:
        """A random non-constant XOR of `vectors`."""
        while True:
            f = 0
            for v in vectors:
                if self.rng.random() < 0.5:
                    f ^= v
            if f not in (0, self.one):
                return f

    def fill(self, slots: list[int]) -> None:
        """Step 1: each slot gets a product of two random members of the span before it."""
        for s in slots:
            live = self.regs[:s]
            span = TaggedSpan([self.one, *live])
            while (f := self.random_in(live)) & (g := self.random_in(live)) in span:
                pass
            self.product_into(f, g, s)

    def plan(self, s: int, remaining: list[int]) -> tuple[list[tuple[int, int]], int]:
        """Step 2 for slot s: the products, and the V_s part to CX away after them."""
        vs = [self.one, *(v for k, v in enumerate(self.regs) if k != s)]
        spans = SlotSpans()
        for i, v in enumerate(vs):
            spans.add((VS, i), v)
        for j, t in enumerate(remaining):
            spans.add((TARGET, j), t, target_only=True)
        fs: list[int] = []
        while (tag := self.landing(self.regs[s], spans, remaining)) is None:
            if len(fs) > 4 * len(vs):
                raise RuntimeError(f"slot {s} does not converge")
            pool = [self.random_in(vs[1:]) for _ in range(CANDIDATES)]
            fs.append(max(pool, key=lambda f: sum(f & g not in spans.without_t for g in vs)))
            for i, g in enumerate(vs):
                spans.add((len(fs) - 1, i), fs[-1] & g)
        g_of, v = spans.decode(tag, vs, len(fs))
        return [(f, g) for f, g in zip(fs, g_of) if f & g], v

    @staticmethod
    def landing(h: int, spans: SlotSpans, remaining: list[int]) -> int | None:
        """The tag decomposing slot value h, once it lies in V_s + products + T outside V_s +
        products, or differs from an unplaced target by a member of V_s + products."""
        rest, tag = spans.with_t.reduce(h)
        if rest == 0 and h not in spans.without_t:
            return tag
        for t in remaining:
            rest, tag = spans.without_t.reduce(h ^ t)
            if rest == 0:
                return tag
        return None

    def finish(self, slots: list[int]) -> None:
        """Step 3: CX the slots (now spanning T) onto the independent targets, one pivot slot each,
        then CX every dependent output into its own slot."""
        free = list(slots)
        for j in self.independent:
            _, rs = self.express(self.columns[j])
            s = next(r for r in rs if r in free)
            self.xor_into(self.columns[j] ^ self.regs[s], s)
            free.remove(s)
            self.slot_of[j] = s
        for j in set(range(self.m)) - set(self.independent):
            self.slot_of[j] = self.n + j
            self.xor_into(self.columns[j], self.n + j)

    def run(self) -> int:
        """Synthesise; each step plans SAMPLE waiting slots and applies the cheapest. Returns CCX count."""
        slots = [self.n + j for j in self.independent]
        self.fill(slots)
        waiting = list(slots)
        while waiting:
            placed = TaggedSpan([self.one, *self.regs])
            remaining = [t for t in self.targets if t not in placed]
            plans = [(self.plan(s, remaining), s) for s in self.rng.sample(waiting, min(SAMPLE, len(waiting)))]
            (products, v), s = min(plans, key=lambda p: len(p[0][0]))
            for f, g in products:
                self.product_into(f, g, s)
            self.xor_into(v, s)
            waiting.remove(s)
        self.finish(slots)
        return self.toffolis

    def verify(self) -> bool:
        """Replay `gates` from the initial state: address restored, every output in its slot."""
        regs = [*self.address, *([0] * self.m)]
        for g in self.gates:
            apply(regs, g, self.one)
        return regs[: self.n] == self.address and all(regs[s] == self.columns[j] for j, s in self.slot_of.items())
