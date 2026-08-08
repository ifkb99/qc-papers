"""Walsh-Hadamard analysis of permutation circuits (claim C8).

For a permutation unitary pi with pi|y> = |perm(y)>, back-propagating Z_j gives

    pi^dag Z_j pi |y> = (-1)^{g(y)} |y|,     g(y) = bit j of perm(y)

which is diagonal, hence expands purely in Z-type Pauli strings:

    pi^dag Z_j pi = sum_z c_z Z^z,     c_z = 2^-n sum_y (-1)^{g(y)} (-1)^{y.z}

That last expression is exactly the Walsh-Hadamard transform of (-1)^{g}. So the
number of Pauli terms PPS ends up holding equals the **Walsh sparsity of g** --
a property of the Boolean function being computed, not of the gate set used to
compute it. That is the C8 claim, and this module tests it directly.
"""
from __future__ import annotations
import numpy as np


def classical_permutation(circuit) -> np.ndarray:
    """Replay a classical reversible circuit on all 2^n basis states at once.

    Returns perm with perm[y] = the basis state |y> is mapped to.
    """
    if not circuit.is_classical():
        bad = next(op[0] for op in circuit.logical if op[0] not in ("x", "cnot", "toffoli"))
        raise ValueError(f"circuit is not a permutation (found {bad!r})")

    n = circuit.n
    idx = np.arange(1 << n, dtype=np.int64)
    for op in circuit.logical:
        if op[0] == "x":
            idx ^= np.int64(1 << op[1])
        elif op[0] == "cnot":
            c, t = op[1], op[2]
            idx ^= ((idx >> c) & 1) << t
        else:                                    # toffoli
            a, b, c = op[1], op[2], op[3]
            idx ^= (((idx >> a) & 1) & ((idx >> b) & 1)) << c
    return idx


def permutation_via_statevector(circuit, tol: float = 1e-8) -> np.ndarray:
    """Extract the basis permutation of ANY circuit that happens to implement
    one, without requiring it to be a permutation gate-by-gate.

    Evolves the whole computational basis at once (a 2^n x 2^n identity), so it
    works for Fourier-compiled arithmetic where individual gates are phase
    rotations. Exponential in n -- for cross-checking, not for scale.
    """
    import statevec as sv
    n = circuit.n
    dim = 1 << n
    out = sv.run(circuit, np.eye(dim, dtype=complex))
    mag = np.abs(out)
    img = np.argmax(mag, axis=0)
    peak = mag[img, np.arange(dim)]
    if not np.all(peak > 1 - tol):
        raise ValueError(f"circuit is not a basis permutation "
                         f"(min peak amplitude {peak.min():.4f})")
    return img.astype(np.int64)


def wht(vec: np.ndarray) -> np.ndarray:
    """Unnormalised fast Walsh-Hadamard transform (in-place butterfly)."""
    a = vec.astype(np.float64).copy()
    n = a.size
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        x, y = a[:, 0, :].copy(), a[:, 1, :].copy()
        a[:, 0, :] = x + y
        a[:, 1, :] = x - y
        a = a.reshape(-1)
        h *= 2
    return a


def pullback_coefficients(circuit, target_qubit: int, perm=None) -> np.ndarray:
    """Pauli coefficients c_z of pi^dag Z_target pi, indexed by the z bitmask.

    Uses the fast gate-level replay when the circuit is classical, otherwise
    falls back to extracting the permutation from the state vector.
    """
    if perm is None:
        perm = (classical_permutation(circuit) if circuit.is_classical()
                else permutation_via_statevector(circuit))
    n = circuit.n
    g = (perm >> target_qubit) & 1
    chi = np.where(g == 1, -1.0, 1.0)
    return wht(chi) / (1 << n)


def walsh_sparsity(circuit, target_qubit: int, tol: float = 1e-12) -> int:
    c = pullback_coefficients(circuit, target_qubit)
    return int(np.count_nonzero(np.abs(c) > tol))


def walsh_degree(circuit, target_qubit: int, tol: float = 1e-12) -> int:
    """Max Hamming weight over the Walsh support.

    NOT the algebraic (GF(2)) degree. For an affine function such as
    a0 XOR b0 XOR c0 the Walsh support is the single point {a0,b0,c0}, so this
    returns 3 while the algebraic degree is 1. It measures the maximum *Pauli
    weight* PPS has to carry, which is the quantity that matters here.
    """
    c = pullback_coefficients(circuit, target_qubit)
    nz = np.nonzero(np.abs(c) > tol)[0]
    if nz.size == 0:
        return 0
    return int(max(int(z).bit_count() for z in nz))


def is_affine(circuit, target_qubit: int, tol: float = 1e-12) -> bool:
    """A Boolean function is affine iff its Walsh spectrum has a single point."""
    return walsh_sparsity(circuit, target_qubit, tol) == 1
