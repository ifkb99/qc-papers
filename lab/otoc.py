"""Out-of-time-order correlators of basis permutations via boomerang counts.

For U|y> = |pi(y)>, W = U^dag X^b Z^e U and V = X^a Z^d (Z applied first),

    F = 2^-n Tr(W^dag V^dag W V)
      = 2^-n sum_{y in B} (-1)^(d.(y ^ sigma(y)) ^ e.(pi(y) ^ pi(y ^ a))),
    sigma(y) = pi^-1(pi(y) ^ b),
    B = {y : sigma(y ^ a) ^ sigma(y) = a},

and |B| is the boomerang connectivity table entry beta_pi(a, b) of Cid et al.
(Boura & Canteaut, ToSC 2018(3), Definition 3). `otoc_formula` evaluates
this; `otoc_dense_perm` / `otoc_dense_unitary` are independent references
that multiply explicit matrices and take traces.
"""
from __future__ import annotations

import numpy as np


def parity(v):
    return np.bitwise_count(v) & 1


def otoc_formula(perm, a, d, b, e, *, orientation="pi"):
    perm = np.asarray(perm, dtype=np.int64)
    if orientation != "pi":                     # control: use pi^-1 in place of pi
        perm = np.argsort(perm)
    inv = np.argsort(perm)
    ys = np.arange(len(perm), dtype=np.int64)
    sigma = inv[perm ^ b]
    cond = (sigma[ys ^ a] ^ sigma) == a
    sign = parity((d & (ys ^ sigma)) ^ (e & (perm ^ perm[ys ^ a])))
    return float(np.where(cond, np.where(sign == 1, -1.0, 1.0), 0.0).sum()) / len(perm)


def bct(perm, a, b):
    perm = np.asarray(perm, dtype=np.int64)
    inv = np.argsort(perm)
    ys = np.arange(len(perm), dtype=np.int64)
    return int(((inv[perm ^ b] ^ inv[perm[ys ^ a] ^ b]) == a).sum())


def ddt(perm, a, b):
    perm = np.asarray(perm, dtype=np.int64)
    ys = np.arange(len(perm), dtype=np.int64)
    return int(((perm[ys ^ a] ^ perm) == b).sum())


_X = np.array([[0., 1.], [1., 0.]])
_Z = np.diag([1., -1.])


def pauli_xz(x, z, n):
    """Real matrix X^x Z^z by explicit kron (qubit 0 = least significant)."""
    m = np.eye(1)
    for q in reversed(range(n)):
        f = np.eye(2)
        if (x >> q) & 1:
            f = _X @ f
        if (z >> q) & 1:
            f = f @ _Z if (x >> q) & 1 else _Z
        m = np.kron(m, f)
    return m


def otoc_dense_unitary(U, a, d, b, e):
    n = int(U.shape[0]).bit_length() - 1
    P = pauli_xz(b, e, n)
    V = pauli_xz(a, d, n)
    W = U.conj().T @ P @ U
    return float(np.real(np.trace(W.conj().T @ V.T @ W @ V))) / U.shape[0]


def otoc_dense_perm(perm, a, d, b, e):
    D = len(perm)
    U = np.zeros((D, D))
    U[np.asarray(perm), np.arange(D)] = 1.0
    return otoc_dense_unitary(U, a, d, b, e)
