"""GF(2) structure analysis of Walsh supports.

The affine-aware finder (structures/find_structures) supersedes the bare rank
test from experiment_linstruct.py: a linear structure g(y^w)=g(y) confines the
support to the hyperplane <z,w>=0, which the rank test sees, but an AFFINE
structure g(y^w)=g(y)^1 confines it to the coset <z,w>=1, which can still span
the full space -- invisible to rank. Translating the support by one member and
taking the kernel of the difference set catches both (see NOTES.md SS AF).

Counting lemma worth remembering before running anything: either kind of
structure caps density at 1/2, so density > 1/2 already proves none exists.
"""
from __future__ import annotations
import numpy as np


def rank_kernel(vecs, n: int) -> tuple[int, list[int]]:
    """GF(2) rank of the given bitmask vectors, plus a basis of the orthogonal
    complement {w : <v,w> = 0 for all v}."""
    basis: list[int] = []
    for v in vecs:
        cur = int(v)
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    piv = {b.bit_length() - 1: b for b in basis}
    kernel = []
    for free in range(n):
        if free in piv:
            continue
        w = 1 << free
        for p in sorted(piv):
            if (piv[p] & w).bit_count() % 2:
                w ^= 1 << p
        kernel.append(w)
    return len(basis), kernel


def structures(zs: np.ndarray, n: int) -> list[tuple[int, int]]:
    """All (w, eps) with <z,w> = eps for every support element z.

    eps = 0 -> linear structure of g;  eps = 1 -> affine structure.
    Empty support returns [].
    """
    if zs.size == 0:
        return []
    z0 = int(zs[0])
    _, kernel = rank_kernel((np.asarray(zs, dtype=np.int64) ^ z0).tolist(), n)
    return [(w, (z0 & w).bit_count() & 1) for w in kernel]


def check_structure(g: np.ndarray, w: int, eps: int) -> bool:
    """Pointwise verification g(y ^ w) == g(y) ^ eps over the full table."""
    idx = np.arange(g.size, dtype=np.int64)
    return bool(np.array_equal(g[idx ^ w], g ^ eps))


def find_structures(g: np.ndarray, n: int, tol: float = 1e-12,
                    verify: bool = True):
    """Convenience: spectrum -> support -> candidate structures, each
    pointwise-verified against the table unless verify=False.

    Returns (support, [(w, eps), ...]).
    """
    import lab.measure as measure
    zs = measure.fn_support(g, tol)
    out = structures(zs, n)
    if verify:
        bad = [(w, e) for w, e in out if not check_structure(g, w, e)]
        if bad:
            raise AssertionError(
                f"support-derived structures fail pointwise: {bad} "
                f"(spectrum/table mismatch -- suspect a bug)")
    return zs, out


def coset_split(zs: np.ndarray, w: int) -> tuple[int, int]:
    """(#support with <z,w>=0, #support with <z,w>=1)."""
    par = np.array([(int(z) & w).bit_count() & 1 for z in zs.tolist()])
    return int(np.count_nonzero(par == 0)), int(np.count_nonzero(par == 1))


def quadrant_counts(zs: np.ndarray, bit_i: int, bit_j: int) -> dict:
    """Support counts by the ((z>>bit_i)&1, (z>>bit_j)&1) quadrant."""
    bi = (zs >> bit_i) & 1
    bj = (zs >> bit_j) & 1
    return {(int(a), int(b)): int(np.count_nonzero((bi == a) & (bj == b)))
            for a in (0, 1) for b in (0, 1)}


def slice_fn(g: np.ndarray, n: int, bit: int, val: int) -> np.ndarray:
    """Restrict g to {y : y_bit = val}, as a table over the remaining n-1
    bits (low bits keep their positions, higher bits shift down by one)."""
    idx = np.arange(1 << (n - 1), dtype=np.int64)
    lo = idx & ((1 << bit) - 1)
    hi = (idx >> bit) << (bit + 1)
    return g[hi | (np.int64(val) << bit) | lo]
