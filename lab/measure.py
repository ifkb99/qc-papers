"""Measurement helpers with a content-addressed cache.

The cache key is a hash of the circuit's gate content plus the observable, so
a changed circuit is a *different key* -- stale-cache bugs are structurally
impossible, which matters in a repo where a silently wrong number costs days.
Cached artifacts live in out/cache/ (gitignored). Set LAB_NO_CACHE=1 to
bypass, or delete out/cache/ to clear.

Function-level helpers (fn_spectrum / fn_support) are uncached -- WHTs of
period-r tables cost milliseconds.
"""
from __future__ import annotations
import hashlib
import os
from pathlib import Path

import numpy as np

import accel
import walsh
import perm_pps

TOL = 1e-12
CACHE_DIR = Path(__file__).resolve().parent.parent / "out" / "cache"


def _cache_enabled() -> bool:
    return os.environ.get("LAB_NO_CACHE", "") != "1"


def _key(circuit, target_qubit: int) -> str:
    h = hashlib.sha256()
    h.update(f"n={circuit.n};t={target_qubit};".encode())
    h.update(repr(circuit.logical).encode())
    h.update(repr(circuit.gates).encode())
    return h.hexdigest()[:32]


# -- circuit level ----------------------------------------------------------

def pullback(circuit, target_qubit: int) -> np.ndarray:
    """Walsh coefficients of the pulled-back Z_target. Uncached passthrough."""
    return walsh.pullback_coefficients(circuit, target_qubit)


def support(circuit, target_qubit: int, tol: float = TOL) -> np.ndarray:
    """Sorted int64 array of z bitmasks in the Walsh support. Cached."""
    if _cache_enabled():
        path = CACHE_DIR / f"supp_{_key(circuit, target_qubit)}.npy"
        if path.exists():
            return np.load(path)
    if accel.enabled() and circuit.n <= accel.MAX_QUBITS and circuit.is_classical():
        # Opt-in only (LAB_GPU=1). Gated against the CPU reference by
        # test_accel.py; never silently substituted.
        zs = accel.pullback_support(circuit, target_qubit, tol)
    else:
        c = walsh.pullback_coefficients(circuit, target_qubit)
        zs = np.nonzero(np.abs(c) > tol)[0].astype(np.int64)
    if _cache_enabled():
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        np.save(path, zs)
    return zs


def sparsity(circuit, target_qubit: int, tol: float = TOL) -> int:
    return int(support(circuit, target_qubit, tol).size)


def density(circuit, target_qubit: int, tol: float = TOL) -> float:
    return support(circuit, target_qubit, tol).size / (1 << circuit.n)


def peak_pps(circuit, zmask: int, **kwargs) -> perm_pps.PermPPSResult:
    """Permutation-native propagation (peak memory, final terms, <O>).
    Uncached: results are small and the caller usually wants the trajectory."""
    return perm_pps.propagate_perm(circuit, zmask, **kwargs)


# -- function level ---------------------------------------------------------

def fn_spectrum(g: np.ndarray) -> np.ndarray:
    """Normalised Walsh coefficients of a 0/1 truth table of length 2^n."""
    return walsh.wht(np.where(g == 1, -1.0, 1.0)) / g.size


def fn_support(g: np.ndarray, tol: float = TOL) -> np.ndarray:
    c = fn_spectrum(g)
    return np.nonzero(np.abs(c) > tol)[0].astype(np.int64)
