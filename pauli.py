"""Symplectic Pauli-string algebra for Pauli Path Simulation (PPS).

A Pauli string on n qubits is stored as a pair of bitmasks (x, z):

    P(x, z) = i^{|x & z|} * (X^x Z^z)

The i^{|x&z|} factor makes P Hermitian (it is what turns XZ into Y), so the
coefficients in the Heisenberg expansion stay real.

Every gate is expressed as a Pauli rotation U = exp(-i*theta*sigma/2), which
gives a single conjugation rule (Gharibyan et al. arXiv:2507.10771, Eq. 8):

    U^dag P U = P                              if [P, sigma] = 0
              = cos(t) P + sin(t) (i sigma P)   if {P, sigma} = 0

Clifford gates are just the theta = pi/2 case, where cos(t) = 0 and the term
maps to a single Pauli -- no branching. All branching comes from non-Clifford
angles. That is the whole mechanism PPS cost depends on.
"""
from __future__ import annotations
import numpy as np

_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)


def popcount(v: int) -> int:
    return v.bit_count()


def commutes(p: tuple[int, int], q: tuple[int, int]) -> bool:
    """Symplectic form: P and Q commute iff x1.z2 + z1.x2 == 0 (mod 2)."""
    x1, z1 = p
    x2, z2 = q
    return (popcount(x1 & z2) + popcount(z1 & x2)) % 2 == 0


def pauli_mult(p: tuple[int, int], q: tuple[int, int]) -> tuple[tuple[int, int], int]:
    """Return (R, k) such that P(p) @ P(q) == i^k * P(R).

    From (X^a Z^b)(X^x Z^z) = (-1)^{b.x} X^{a^x} Z^{b^z}, corrected for the
    Hermitian normalisation of all three strings.
    """
    a, b = p
    x, z = q
    rx, rz = a ^ x, b ^ z
    k = popcount(a & b) + popcount(x & z) - popcount(rx & rz) + 2 * popcount(b & x)
    return (rx, rz), k % 4


def i_sigma_p(sigma: tuple[int, int], p: tuple[int, int]) -> tuple[tuple[int, int], float]:
    """Return (R, s) with i * P(sigma) @ P(p) == s * P(R), s = +/-1.

    Only called when sigma and p anticommute, in which case i*sigma*p is
    Hermitian and s is real.
    """
    r, k = pauli_mult(sigma, p)
    k = (k + 1) % 4          # the leading factor of i
    if k == 0:
        return r, 1.0
    if k == 2:
        return r, -1.0
    raise ValueError(f"i*sigma*P not Hermitian (i^{k}); operands must anticommute")


def weight(p: tuple[int, int]) -> int:
    """Number of qubits the string acts on non-trivially."""
    x, z = p
    return popcount(x | z)


def to_matrix(p: tuple[int, int], n: int) -> np.ndarray:
    """Dense matrix of P(x, z). Testing only -- exponential in n."""
    x, z = p
    M = np.array([[1]], dtype=complex)
    for q in range(n):
        xb, zb = (x >> q) & 1, (z >> q) & 1
        m = _I if not xb and not zb else _X if xb and not zb else _Z if zb and not xb else _Y
        M = np.kron(m, M)          # qubit 0 is the least significant bit
    return M


def rotation_matrix(sigma: tuple[int, int], theta: float, n: int) -> np.ndarray:
    """Dense exp(-i theta P(sigma) / 2). Testing only."""
    S = to_matrix(sigma, n)
    return np.cos(theta / 2) * np.eye(2 ** n) - 1j * np.sin(theta / 2) * S
