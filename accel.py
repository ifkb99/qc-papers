"""Optional CUDA backend for the two hot paths.

Everything expensive in this project is one of two memory-bound array passes:

  * `walsh.classical_permutation` -- one elementwise int64 pass per logical gate
    over 2^n entries. Dominant cost for circuit-level sweeps (a 25-qubit modexp
    is ~6000 gates over 33.5M entries).
  * `walsh.wht` -- n strided add/sub passes over 2^n float64.

Both are bandwidth-limited, so a card with ~10x the memory bandwidth of system
RAM wins by roughly that factor. Measured on an RTX A4500 (20 GiB):

    perm replay q=17 (3373 gates)    0.47 s -> 0.07 s     7.1x
    perm replay q=19 (5047 gates)    3.28 s -> 0.47 s     6.9x
    perm replay q=21 (6721 gates)   29.31 s -> 2.56 s    11.5x
    FWHT 2^24                        0.86 s -> 0.07 s    11.8x
    FWHT 2^26                        3.77 s -> 0.29 s    12.8x
    FWHT 2^28                       18.12 s -> 1.23 s    14.7x

The speedup GROWS with size, which is the useful direction: it is the large
runs that hurt.

DESIGN RULE, deliberately conservative. This module NEVER changes a default.
`walsh.py` is the verified reference and stays untouched; callers opt in
explicitly, or set LAB_GPU=1 for `lab.measure`. `test_accel.py` gates every GPU
routine against its CPU counterpart -- exact equality for permutations,
allclose for transforms. A faster replacement for a verified primitive is
exactly how a headline number gets silently corrupted, so the GPU path is
treated as a claim to be checked rather than an optimisation to be trusted.

CAPACITY, measured. On one 20 GiB card the float64 FWHT fits n = 30 (8 GiB
array + 4 GiB temp = 12 GiB) and fails at n = 31 (24 GiB). Memory doubles per
qubit, so a second card buys exactly ONE more qubit -- see `wht_exact` for a
cheaper lever that buys the same thing without the complexity, by exploiting
the fact that the transform of +/-1 data is exactly integer-valued.
"""
from __future__ import annotations

import os

import numpy as np

try:
    import cupy as _cp
    HAVE_GPU = _cp.cuda.runtime.getDeviceCount() > 0
except Exception:                                    # no cupy, no driver, no card
    _cp = None
    HAVE_GPU = False

# 2^30 int64 = 8 GiB; leave room for a temporary of the same size.
MAX_QUBITS = 30


def enabled() -> bool:
    """GPU available AND opted into (LAB_GPU=1)."""
    return HAVE_GPU and os.environ.get("LAB_GPU", "") == "1"


def device_info() -> str:
    if not HAVE_GPU:
        return "no CUDA device"
    props = _cp.cuda.runtime.getDeviceProperties(0)
    free, total = _cp.cuda.Device(0).mem_info
    return (f"{props['name'].decode()} x{_cp.cuda.runtime.getDeviceCount()}, "
            f"{total / 2**30:.1f} GiB, {free / 2**30:.1f} GiB free")


def _sync():
    _cp.cuda.Stream.null.synchronize()


def free_pool():
    """Release cupy's cached blocks -- call between large independent runs."""
    if HAVE_GPU:
        _cp.get_default_memory_pool().free_all_blocks()


# ---------------------------------------------------------------------------

def classical_permutation(circuit) -> np.ndarray:
    """GPU replay of a classical circuit; returns a host array.

    Mirrors walsh.classical_permutation gate for gate, including its error on
    a non-classical op, so the two cannot silently diverge.
    """
    if not HAVE_GPU:
        raise RuntimeError("no CUDA device available")
    if not circuit.is_classical():
        bad = next(op[0] for op in circuit.logical
                   if op[0] not in ("x", "cnot", "toffoli"))
        raise ValueError(f"circuit is not a permutation (found {bad!r})")
    n = circuit.n
    if n > MAX_QUBITS:
        raise MemoryError(f"n={n} exceeds MAX_QUBITS={MAX_QUBITS}")

    idx = _cp.arange(1 << n, dtype=_cp.int64)
    for op in circuit.logical:
        if op[0] == "x":
            idx ^= (1 << op[1])
        elif op[0] == "cnot":
            c, t = op[1], op[2]
            idx ^= ((idx >> c) & 1) << t
        else:
            a, b, c = op[1], op[2], op[3]
            idx ^= (((idx >> a) & 1) & ((idx >> b) & 1)) << c
    out = _cp.asnumpy(idx)
    del idx
    free_pool()
    return out


def wht(vec) -> np.ndarray:
    """Unnormalised FWHT on the GPU; returns a host array.

    In-place butterfly with a single half-size temporary, rather than the two
    full-size copies the CPU reference allocates -- the memory saving is what
    makes n = 30 fit alongside the input.
    """
    if not HAVE_GPU:
        raise RuntimeError("no CUDA device available")
    a = _cp.asarray(vec, dtype=_cp.float64)
    n = a.size
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        t = a[:, 0, :] - a[:, 1, :]
        a[:, 0, :] += a[:, 1, :]
        a[:, 1, :] = t
        del t
        a = a.reshape(-1)
        h *= 2
    out = _cp.asnumpy(a)
    del a
    free_pool()
    return out


def wht_exact(chi_pm1) -> np.ndarray:
    """EXACT integer FWHT of a +/-1 vector, on the GPU.

    The transform of +/-1 data is a sum/difference of integers at every level,
    so it is exactly integer-valued -- verified bit-for-bit against the float64
    reference. Consequences, both material:

      * **half the memory** of float64, so n = 30 costs 6 GiB rather than 12;
      * **the support test becomes exact** (!= 0) instead of a magnitude
        threshold, which removes the thresholding artifact that otherwise makes
        densities drift below their true value as n grows.

    Intermediate magnitudes are bounded by 2^k after k levels, so int32 is safe
    for n <= 30 (2^30 < 2^31) and this function refuses beyond that.
    """
    if not HAVE_GPU:
        raise RuntimeError("no CUDA device available")
    a = _cp.asarray(chi_pm1, dtype=_cp.int32)
    n = a.size
    if n > (1 << 30):
        raise ValueError("int32 FWHT is only safe to n = 30; use wht() for larger")
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        t = a[:, 0, :] - a[:, 1, :]
        a[:, 0, :] += a[:, 1, :]
        a[:, 1, :] = t
        del t
        a = a.reshape(-1)
        h *= 2
    out = _cp.asnumpy(a)
    del a
    free_pool()
    return out


def pullback_support_exact(circuit, target_qubit: int, perm=None) -> np.ndarray:
    """Walsh support via exact integer arithmetic -- no tolerance parameter.

    Preferred over `pullback_support` whenever n <= 30: same answer, half the
    memory, and no threshold to misclassify a genuinely tiny coefficient.
    """
    if not HAVE_GPU:
        raise RuntimeError("no CUDA device available")
    n = circuit.n
    if n > 30:
        raise ValueError(f"n={n} > 30; use pullback_support (float64)")
    if perm is None:
        if not circuit.is_classical():
            raise ValueError("circuit is not a permutation")
        idx = _cp.arange(1 << n, dtype=_cp.int64)
        for op in circuit.logical:
            if op[0] == "x":
                idx ^= (1 << op[1])
            elif op[0] == "cnot":
                idx ^= ((idx >> op[1]) & 1) << op[2]
            else:
                idx ^= (((idx >> op[1]) & 1) & ((idx >> op[2]) & 1)) << op[3]
    else:
        idx = _cp.asarray(perm)
    g = ((idx >> target_qubit) & 1).astype(_cp.int32)
    del idx
    a = 1 - 2 * g                                   # 0 -> +1, 1 -> -1
    del g
    size = a.size
    h = 1
    while h < size:
        a = a.reshape(-1, 2, h)
        t = a[:, 0, :] - a[:, 1, :]
        a[:, 0, :] += a[:, 1, :]
        a[:, 1, :] = t
        del t
        a = a.reshape(-1)
        h *= 2
    zs = _cp.nonzero(a)[0].astype(_cp.int64)        # EXACT: no tolerance
    del a
    out = _cp.asnumpy(zs)
    del zs
    free_pool()
    return out


def pullback_support(circuit, target_qubit: int, tol: float = 1e-12,
                     perm=None) -> np.ndarray:
    """Walsh support of the pulled-back Z_target, end to end on the GPU.

    The whole point: the intermediate 2^n coefficient array never leaves the
    card, so only the (small) support is transferred back.
    """
    if not HAVE_GPU:
        raise RuntimeError("no CUDA device available")
    n = circuit.n
    if perm is None:
        if not circuit.is_classical():
            raise ValueError("circuit is not a permutation")
        idx = _cp.arange(1 << n, dtype=_cp.int64)
        for op in circuit.logical:
            if op[0] == "x":
                idx ^= (1 << op[1])
            elif op[0] == "cnot":
                idx ^= ((idx >> op[1]) & 1) << op[2]
            else:
                idx ^= (((idx >> op[1]) & 1) & ((idx >> op[2]) & 1)) << op[3]
    else:
        idx = _cp.asarray(perm)

    g = (idx >> target_qubit) & 1
    del idx
    a = _cp.where(g == 1, -1.0, 1.0)
    del g
    size = a.size
    h = 1
    while h < size:
        a = a.reshape(-1, 2, h)
        t = a[:, 0, :] - a[:, 1, :]
        a[:, 0, :] += a[:, 1, :]
        a[:, 1, :] = t
        del t
        a = a.reshape(-1)
        h *= 2
    a /= size
    zs = _cp.nonzero(_cp.abs(a) > tol)[0].astype(_cp.int64)
    del a
    out = _cp.asnumpy(zs)
    del zs
    free_pool()
    return out
