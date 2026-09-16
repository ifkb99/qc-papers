"""Exact ROBDD node counts of a Boolean function given by its truth table.

This is a counting routine, not a propagator: it never touches gates. The
truth table comes from an existing verified permutation replay
(`walsh.classical_permutation` or `accel.classical_permutation`), and the
count is obtained by bottom-up hash-consing of cofactor pairs, one variable
at a time, in a stated order.

Conventions (both are reported, because published tools differ):

  * `plain`   : two terminals, no complemented edges. A node at variable v is
                a distinct pair (low, high) with low != high.
  * `complement` : complemented edges; a function and its negation share a
                node, and there is one terminal. Equivalently, the number of
                negation classes of distinct non-constant subfunctions that
                depend on their top variable, plus one terminal. The class
                count does not depend on which edge a package keeps regular.

Order: `order[0]` is the top (first-read) variable. Variables are qubit
indices, i.e. bit positions of the truth-table index.
"""
from __future__ import annotations

import numpy as np


def _arrange(table: np.ndarray, n: int, order) -> np.ndarray:
    """Return the table as a flat array whose C-order index reads
    order[0] as the most significant axis and order[-1] as the least."""
    order = [int(q) for q in order]
    if sorted(order) != list(range(n)):
        raise ValueError("order must be a permutation of range(n)")
    t = np.asarray(table).reshape((2,) * n)     # axis k <-> bit n-1-k
    axes = [n - 1 - q for q in order]
    return np.ascontiguousarray(t.transpose(axes)).reshape(-1)


def robdd_size(table, n: int, order, convention: str = "plain",
               per_level: bool = False):
    """Node count of the ROBDD of the 0/1 truth table under `order`.

    Returns an int, or (int, list of per-variable node counts in `order`)
    when per_level is set. Terminals are included in the total only.
    """
    tab = np.asarray(table)
    if tab.size != 1 << n:
        raise ValueError("table size must be 2^n")
    flat = _arrange(tab.astype(np.int8, copy=False), n, order)
    if convention == "plain":
        ids = flat.astype(np.int64)             # terminals 0 and 1
        next_id = 2
        levels = []
        for _ in range(n):                      # bottom variable first
            pairs = ids.reshape(-1, 2)
            lo, hi = pairs[:, 0], pairs[:, 1]
            red = lo == hi
            key = lo * (1 << 32) + hi           # ids stay < 2^32 for n <= 30
            new_ids = lo.copy()
            if (~red).any():
                uniq, inv = np.unique(key[~red], return_inverse=True)
                new_ids[~red] = inv + next_id
                next_id += uniq.size
                levels.append(int(uniq.size))
            else:
                levels.append(0)
            ids = new_ids
        terminals = len(np.unique(flat))
        total = sum(levels) + terminals
    elif convention == "complement":
        # signed ids: terminal TRUE = +1, FALSE = -1; negation flips the sign.
        ids = np.where(flat.astype(bool), 1, -1).astype(np.int64)
        next_id = 2
        levels = []
        for _ in range(n):
            pairs = ids.reshape(-1, 2)
            lo, hi = pairs[:, 0], pairs[:, 1]
            red = lo == hi
            # normalize so the high edge is regular (positive)
            sgn = np.where(hi < 0, -1, 1)
            nlo, nhi = lo * sgn, hi * sgn
            off = np.int64(1 << 31)
            key = (nlo + off) * (1 << 32) + (nhi + off)
            new_ids = lo.copy()
            if (~red).any():
                uniq, inv = np.unique(key[~red], return_inverse=True)
                new_ids[~red] = (inv + next_id) * sgn[~red]
                next_id += uniq.size
                levels.append(int(uniq.size))
            else:
                levels.append(0)
            ids = new_ids
        total = sum(levels) + 1
    else:
        raise ValueError("convention must be 'plain' or 'complement'")
    levels = levels[::-1]                       # report in `order`
    return (total, levels) if per_level else total


def output_bit_table(perm: np.ndarray, j: int) -> np.ndarray:
    """0/1 truth table of output bit j of a basis permutation (the pulled-back
    projector onto |1> at qubit j; C8's object as a Boolean function)."""
    return ((np.asarray(perm) >> j) & 1).astype(np.int8)
