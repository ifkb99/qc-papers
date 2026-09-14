"""Small dense tensor diagnostics; NOT a propagator or a GF(2) rank routine.

Entries use little-endian bit indexing. Ranks are numerical over R/C, with
explicit tolerances. Construction and SVD remain exponential in qubit count;
these diagnostics do not themselves supply a compressed simulation algorithm.
"""
from __future__ import annotations
from operator import index
import numpy as np


def cut_matrix(values, left_qubits):
    values = np.asarray(values)
    if values.ndim != 1 or not values.size or values.size & (values.size - 1):
        raise ValueError("values must be a nonempty power-of-two vector")
    n = values.size.bit_length() - 1
    left = [index(q) for q in left_qubits]
    if len(set(left)) != len(left) or any(q < 0 or q >= n for q in left):
        raise ValueError("left_qubits must be distinct valid indices")
    right = [q for q in range(n) if q not in left]

    def indices(qubits):
        ids = np.arange(1 << len(qubits), dtype=np.int64)
        packed = np.zeros_like(ids)
        for bit, q in enumerate(qubits):
            packed |= ((ids >> bit) & 1) << q
        return packed

    return values[indices(left)[:, None] | indices(right)[None, :]]


def cut_spectrum(values, left_qubits):
    return np.linalg.svd(cut_matrix(values, left_qubits), compute_uv=False)


def numerical_rank(singular_values, rtol=1e-10, atol=1e-12):
    if rtol < 0 or atol < 0:
        raise ValueError("rank tolerances must be nonnegative")
    s = np.asarray(singular_values)
    threshold = max(atol, rtol * float(s[0])) if s.size else atol
    return int(np.count_nonzero(s > threshold))


def rank_profile(values, order=None, rtol=1e-10, atol=1e-12):
    values = np.asarray(values)
    # Also validates the vector when there are no nontrivial cuts.
    cut_matrix(values, [])
    n = values.size.bit_length() - 1
    order = list(range(n)) if order is None else [index(q) for q in order]
    if sorted(order) != list(range(n)):
        raise ValueError("order must be a permutation of all qubits")
    rows = []
    for k in range(1, n):
        s = cut_spectrum(values, order[:k])
        r = numerical_rank(s, rtol, atol)
        rows.append(dict(left=order[:k], rank=r,
                         rank_loose=numerical_rank(s, rtol*100, atol*100),
                         rank_tight=numerical_rank(s, rtol/100, atol/100),
                         smallest_kept=float(s[r-1]) if r else 0.,
                         largest_discarded=float(s[r]) if r < len(s) else 0.))
    return rows


def residue_automaton(table, bits_msb_first):
    """Evaluate a periodic scalar function through a running residue state.

    A concrete finite-state realization, not a minimality claim. The period
    table is supplied explicitly; building it can already cost O(order).
    """
    table = np.asarray(table)
    if table.ndim != 1 or not table.size:
        raise ValueError("table must be a nonempty vector")
    residue = 0
    for b in bits_msb_first:
        if b not in (0, 1):
            raise ValueError("input symbols must be bits")
        residue = (2 * residue + int(b)) % table.size
    return table[residue]
