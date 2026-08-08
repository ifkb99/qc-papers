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
import json
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


def support(circuit, target_qubit: int, tol: float = TOL,
            exact: bool = False) -> np.ndarray:
    """Sorted int64 array of z bitmasks in the Walsh support. Cached.

    `exact=True` takes the integer FWHT path (`accel.wht_exact`), where the
    support test is `!= 0` rather than `|c| > tol` -- no tolerance to justify,
    and half the memory. It needs a GPU and n <= 30; without one it falls back
    to the thresholded CPU reference, so the flag never changes availability,
    only rigour. Cached under a distinct key so the two can be compared.
    """
    if _cache_enabled():
        prefix = "suppx" if exact else "supp"
        path = CACHE_DIR / f"{prefix}_{_key(circuit, target_qubit)}.npy"
        if path.exists():
            return np.load(path)
    usable = (accel.enabled() and circuit.n <= accel.MAX_QUBITS
              and circuit.is_classical())
    if usable and exact:
        zs = accel.pullback_support_exact(circuit, target_qubit)
    elif usable:
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


def sparsity(circuit, target_qubit: int, tol: float = TOL,
             exact: bool = False) -> int:
    return int(support(circuit, target_qubit, tol, exact).size)


def density(circuit, target_qubit: int, tol: float = TOL,
            exact: bool = False) -> float:
    return support(circuit, target_qubit, tol, exact).size / (1 << circuit.n)


def stats(circuit, target_qubit: int, masks=()) -> dict:
    """|support|, density, and per-mask GF(2) parity counts. Cached (JSON).

    Same numbers `support()` would give, but the support array never reaches
    the host -- at q = 30 it is 4 GiB, which is worth neither the transfer nor
    the cache file when the caller only wants a density. `odd[w] == 0` says w
    is a linear structure of the pulled-back bit function (C30).

    GPU-only, n <= 30. The CPU fallback computes the same values the slow way
    so a cardless machine still runs, just not at these sizes.
    """
    masks = tuple(int(w) for w in masks)
    if _cache_enabled():
        h = hashlib.sha256(f"{_key(circuit, target_qubit)};{masks}".encode())
        path = CACHE_DIR / f"stats_{h.hexdigest()[:32]}.json"
        if path.exists():
            d = json.loads(path.read_text())
            d["odd"] = {int(k): v for k, v in d["odd"].items()}
            return d
    if accel.enabled() and circuit.n <= 30 and circuit.is_classical():
        out = accel.pullback_stats(circuit, target_qubit, masks)
    else:
        # Same definition as the GPU path, not the thresholded one: the FWHT of
        # +/-1 data is exactly integer-valued and float64 is exact on integers
        # to 2^53, so `!= 0` on the UNNORMALISED transform is the exact test.
        perm = (walsh.classical_permutation(circuit) if circuit.is_classical()
                else walsh.permutation_via_statevector(circuit))
        chi = np.where(((perm >> target_qubit) & 1) == 1, -1.0, 1.0)
        zs = np.nonzero(walsh.wht(chi))[0].astype(np.int64)
        out = {"count": int(zs.size),
               "density": zs.size / (1 << circuit.n),
               "odd": {w: int(sum((int(z) & w).bit_count() & 1
                                  for z in zs.tolist())) for w in masks}}
    if _cache_enabled():
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(out))
    return out


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
