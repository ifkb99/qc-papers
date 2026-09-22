"""Exact coherent probability queries for one explicit path-predicate Grover family.

Query-level research primitive; no arbitrary oracle or output sampler.
See notes/MC-memory-structure-pilots.md for scope, validation and measurements.
"""
from __future__ import annotations

from fractions import Fraction


def validate(n: int, steps: int, masks: tuple[int, ...]) -> None:
    if not isinstance(n, int) or not 2 <= n <= 22:
        raise ValueError("supported pilot domain is 2 <= n <= 22")
    if not isinstance(steps, int) or not 0 <= steps <= 8:
        raise ValueError("supported pilot domain is 0 <= steps <= 8")
    if any(not isinstance(y, int) or not 0 <= y < (1 << n) for y in masks):
        raise ValueError("mask outside the register")


def oracle(n: int, x: int) -> bool:
    if not 2 <= n <= 22 or not 0 <= x < (1 << n):
        raise ValueError("oracle argument outside pilot domain")
    return bool(x & 1) and not bool(x & (x >> 1))


def query_masks(n: int) -> tuple[int, ...]:
    validate(n, 0, ())
    proposed = (0, 1, 2, 4, (1 << n) - 1,
                sum(1 << i for i in range(0, n, 2)))
    return tuple(dict.fromkeys(y for y in proposed if y < (1 << n)))


def signed_count(n: int, mask: int) -> int:
    """Sum f(x)(-1)^<x,mask>, via the previous-bit automaton."""
    validate(n, 0, (mask,))
    ending_zero, ending_one = 0, (-1 if mask & 1 else 1)
    for i in range(1, n):
        sign = -1 if (mask >> i) & 1 else 1
        ending_zero, ending_one = (ending_zero + ending_one,
                                   sign * ending_zero)
    return ending_zero + ending_one


def coefficients(n: int, marked: int, steps: int) -> tuple[Fraction, Fraction]:
    validate(n, steps, ())
    if not isinstance(marked, int) or not 0 <= marked <= (1 << n):
        raise ValueError("invalid marked count")
    p = Fraction(marked, 1 << n)
    a, b = Fraction(1), Fraction(0)
    for _ in range(steps):
        a, b = (1 - 4 * p) * a - 2 * p * b, 2 * a + b
    return a, b


def probabilities_from_counts(n, steps, masks, counts):
    """Exact rational probabilities, including the coherent y=0 cross term."""
    a, b = coefficients(n, counts[0], steps)
    return tuple((a * int(y == 0) + b * Fraction(counts[y], 1 << n)) ** 2
                 for y in masks)


def compact(n: int, steps: int, masks: tuple[int, ...]):
    validate(n, steps, masks)
    counts = {y: signed_count(n, y) for y in dict.fromkeys((0,) + masks)}
    return probabilities_from_counts(n, steps, masks, counts)


def streaming(n: int, steps: int, masks: tuple[int, ...], block: int = 16384):
    """Stronger exhaustive baseline: blockwise signed counts, no state vector.

    Shares the coefficient identity with compact; is not its independent
    reference. All input enumeration and count construction is charged.
    """
    import numpy as np
    validate(n, steps, masks)
    if not isinstance(block, int) or block < 1:
        raise ValueError("block must be a positive integer")
    keys = tuple(dict.fromkeys((0,) + masks))
    counts = {y: 0 for y in keys}
    for start in range(0, 1 << n, block):
        x = np.arange(start, min(start + block, 1 << n), dtype=np.uint64)
        selected = ((x & 1) != 0) & ((x & (x >> 1)) == 0)
        for y in keys:
            parity = np.zeros(x.size, dtype=np.uint8)
            bits = y
            while bits:
                bit = bits & -bits
                parity ^= ((x & bit) != 0)
                bits ^= bit
            counts[y] += int(np.count_nonzero(selected & (parity == 0)))
            counts[y] -= int(np.count_nonzero(selected & (parity != 0)))
    return probabilities_from_counts(n, steps, masks, counts)


def dense(n: int, steps: int, masks: tuple[int, ...]):
    """Independent vector updates and existing lab Walsh transform.

    No DP, marked count, coefficient recurrence or supplied solution list.
    This reference shares the written oracle specification only.
    """
    import numpy as np
    from walsh import wht
    validate(n, steps, masks)
    x = np.arange(1 << n, dtype=np.uint64)
    selected = (x % 2 == 1)
    for i in range(n - 1):
        selected &= ~((((x >> i) & 1) == 1) & (((x >> (i + 1)) & 1) == 1))
    del x
    # Store sqrt(N)*psi, avoiding irrational normalization until after H.
    state = np.ones(1 << n, dtype=np.float64)
    for _ in range(steps):
        state[selected] *= -1
        mean = state.mean()
        state *= -1
        state += 2 * mean
    transformed = wht(state) / (1 << n)
    return tuple(float(transformed[y] ** 2) for y in masks)
