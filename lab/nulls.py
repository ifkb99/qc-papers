"""Null models with the known sampling traps designed out.

Two of this project's own methodology bugs (NOTES.md SS I) live here as
defaults so they cannot recur:

  * random tables REJECT CONSTANTS -- at small r a meaningful fraction of
    random tables are constant, which contaminated the small-r rows once;
  * there is no "fresh table per sweep point" helper: draw ONE table, sweep
    the other parameter (vary exactly one parameter).

planted_structure() exists to validate structure finders: a finder that
cannot recover a planted (w, eps) has no business reporting absences.
"""
from __future__ import annotations
import numpy as np


def _rng(rng_or_seed) -> np.random.Generator:
    if isinstance(rng_or_seed, np.random.Generator):
        return rng_or_seed
    return np.random.default_rng(rng_or_seed)


def random_boolean(n: int, rng=0) -> np.ndarray:
    """Uniform random truth table over 2^n points (generically structureless,
    density -> 1: the standard negative control)."""
    return _rng(rng).integers(0, 2, size=1 << n).astype(np.int64)


def random_table(r: int, rng=0, reject_constant: bool = True) -> np.ndarray:
    """Random period-r table; constants rejected by default (see module doc)."""
    g = _rng(rng)
    h = g.integers(0, 2, size=r).astype(np.int64)
    while reject_constant and not (0 < int(h.sum()) < r):
        h = g.integers(0, 2, size=r).astype(np.int64)
    return h


def planted_structure(n: int, rng=0, bit: int | None = None,
                      eps: int = 1) -> tuple[np.ndarray, int, int]:
    """Truth table with an exact planted structure at w = 1<<bit.

    g(y) = f(y with `bit` removed) ^ (eps * y_bit), so w = 1<<bit satisfies
    g(y^w) = g(y)^eps by construction: eps=1 plants an affine structure,
    eps=0 a linear one (a dead variable). Returns (g, w, eps).
    """
    if bit is None:
        bit = n - 1
    f = _rng(rng).integers(0, 2, size=1 << (n - 1)).astype(np.int64)
    idx = np.arange(1 << n, dtype=np.int64)
    lo = idx & ((1 << bit) - 1)
    hi = (idx >> (bit + 1)) << bit
    yb = (idx >> bit) & 1
    g = f[hi | lo] ^ (eps * yb)
    return g, 1 << bit, eps
