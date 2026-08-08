"""Sparse state-vector simulation of Pauli-rotation circuits.

`Circuit.to_unitary()` builds a dense 2^n x 2^n matrix and dies past ~12 qubits.
Modular exponentiation needs ~18. Applying a Pauli rotation directly to the
state vector is O(2^n) per gate instead of O(4^n).

For P = i^{|x&z|} X^x Z^z:
    (P psi)[j] = i^{|x&z|} * (-1)^{popcount((j^x) & z)} * psi[j^x]
and exp(-i t P / 2) psi = cos(t/2) psi - i sin(t/2) (P psi).
"""
from __future__ import annotations
import numpy as np


def apply_pauli(psi: np.ndarray, sigma: tuple[int, int], n: int) -> np.ndarray:
    """Apply P(sigma). `psi` may be a (2^n,) vector or a (2^n, k) batch of
    column states, which lets a whole basis be evolved in one pass."""
    x, z = sigma
    idx = np.arange(psi.shape[0], dtype=np.int64)
    src = idx ^ x
    signs = np.where(np.bitwise_count(src & z) & 1, -1.0, 1.0)
    phase = 1j ** (bin(x & z).count("1") % 4)
    out = psi[src]
    if out.ndim == 2:
        return phase * signs[:, None] * out
    return phase * signs * out


def apply_rotation(psi: np.ndarray, sigma: tuple[int, int], theta: float,
                   n: int) -> np.ndarray:
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    if abs(s) < 1e-15:
        return c * psi
    return c * psi - 1j * s * apply_pauli(psi, sigma, n)


def run(circuit, psi: np.ndarray | None = None) -> np.ndarray:
    n = circuit.n
    if psi is None:
        psi = np.zeros(1 << n, dtype=complex)
        psi[0] = 1.0
    for sigma, theta in circuit.gates:
        psi = apply_rotation(psi, sigma, theta, n)
    return psi


def basis(n: int, index: int) -> np.ndarray:
    psi = np.zeros(1 << n, dtype=complex)
    psi[index] = 1.0
    return psi


def pack(assignments: dict[int, int]) -> int:
    """{qubit: bit} -> basis index."""
    v = 0
    for q, b in assignments.items():
        if b:
            v |= 1 << q
    return v


def read_register(index: int, qubits: list[int]) -> int:
    """Integer value held by `qubits` (qubits[0] = least significant)."""
    return sum(((index >> q) & 1) << i for i, q in enumerate(qubits))


def peak(psi: np.ndarray) -> tuple[int, float]:
    j = int(np.argmax(np.abs(psi)))
    return j, float(abs(psi[j]))
