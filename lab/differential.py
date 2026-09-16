"""Differential (DDT) machinery for off-diagonal pullbacks -- C99 (GF(2)^n
Pauli basis) and C100 (Weyl bases of finite abelian groups, lower half).

A basis permutation U|y> = |pi(y)> conjugates X^a Z^b into the signed
permutation M|y> = s(y)|sigma(y)>, s(y) = (-1)^(b.pi(y)),
sigma(y) = pi^-1(pi(y) ^ a). The coefficient of X^c Z^d is
2^-n sum_{y in D_c} s(y)(-1)^(d.y) with D_c = {y : y ^ sigma(y) = c}, and
|D_c| = DDT_{pi^-1}(a, c). `formula_terms` evaluates that identity;
`dense_terms_gpu` is the independent reference, projecting an explicitly
multiplied dense matrix U^T P U onto the Pauli basis without using sigma.

parity, fwht, gf_mul, gf_pow, power_map, ddt_row, delta_of and
formula_terms were first written inline in
experiments/experiment_differential_bridge.py, which stays as their record.
nbits, formula_counts, dense_terms_gpu and the C100 group helpers were
written here directly.
"""
from __future__ import annotations

import numpy as np

POLY = {3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011,
        8: 0b100011011}


def nbits(perm) -> int:
    return int(len(perm)).bit_length() - 1


def parity(v):
    return np.bitwise_count(v) & 1


def fwht(v: np.ndarray) -> np.ndarray:
    """Unnormalised Walsh-Hadamard transform along the last axis."""
    v = np.array(v, dtype=np.float64, copy=True)
    shape = v.shape
    N = shape[-1]
    h = 1
    while h < N:
        w = v.reshape(-1, N // (2 * h), 2, h)
        a = w[:, :, 0, :].copy()
        w[:, :, 0, :] += w[:, :, 1, :]
        w[:, :, 1, :] = a - w[:, :, 1, :]
        h *= 2
    return v.reshape(shape)


# ---------------------------------------------------------------- GF(2^n) S-boxes
def gf_mul(x: int, y: int, n: int) -> int:
    r = 0
    while y:
        if y & 1:
            r ^= x
        y >>= 1
        x <<= 1
        if (x >> n) & 1:
            x ^= POLY[n]
    return r


def gf_pow(x: int, e: int, n: int) -> int:
    r = 1
    while e:
        if e & 1:
            r = gf_mul(r, x, n)
        x = gf_mul(x, x, n)
        e >>= 1
    return r


def power_map(n: int, e: int) -> np.ndarray:
    p = np.array([gf_pow(x, e, n) if x else 0 for x in range(1 << n)], dtype=np.int64)
    if len(set(p.tolist())) != 1 << n:
        raise ValueError(f"x^{e} is not a permutation of GF(2^{n})")
    return p


# ---------------------------------------------------------------- DDT
def ddt_row(perm, a) -> np.ndarray:
    """Row a of the difference distribution table of perm^-1."""
    inv = np.argsort(perm)
    u = np.arange(len(perm))
    return np.bincount(inv[u ^ a] ^ inv[u], minlength=len(perm))


def delta_of(perm) -> int:
    return max(int(ddt_row(perm, a).max()) for a in range(1, len(perm)))


# ---------------------------------------------------------------- identity (I)
def formula_terms(perm, a, b, *, direction="heis", sign=True):
    """({(c, d): |alpha|}, [(c, |D_c|, support size)]) from identity (I)."""
    n = nbits(perm)
    ys = np.arange(1 << n, dtype=np.int64)
    inv = np.argsort(perm)
    sigma = inv[perm ^ a] if direction == "heis" else perm[inv ^ a]
    s = np.where(parity(perm & b) == 1, -1.0, 1.0) if sign else np.ones(1 << n)
    cvec = ys ^ sigma
    out, classes = {}, []
    for c in np.unique(cvec):
        g = np.zeros(1 << n)
        mask = cvec == c
        g[mask] = s[mask]
        W = fwht(g)
        nz = np.nonzero(np.abs(W) > 0.5)[0]
        classes.append((int(c), int(mask.sum()), len(nz)))
        for d in nz:
            out[(int(c), int(d))] = abs(W[d]) / (1 << n)
    return out, classes


def formula_counts(perm, a, b) -> dict[int, np.ndarray]:
    """{c: sorted nonzero d} -- same identity, compact form for large n."""
    terms, _ = formula_terms(perm, a, b)
    byc: dict[int, list[int]] = {}
    for (c, d) in terms:
        byc.setdefault(c, []).append(d)
    return {c: np.sort(np.array(v, dtype=np.int64)) for c, v in byc.items()}


# ---------------------------------------------------------------- dense reference
def dense_terms_gpu(perm, a, b, device=1):
    """{c: sorted nonzero d} for U^T X^a Z^b U by dense GPU matrix algebra.

    U and P are materialised as dense float64 matrices (P by explicit kron of
    2x2 factors), multiplied with cuBLAS, and projected onto X^c Z^d through
    the diagonal V[c, y] = M[y ^ c, y] and a batched Walsh transform. Nothing
    here uses sigma, D_c or the DDT. Real X^a Z^b has the same support and
    magnitudes as the Hermitian Pauli label (a, b).
    """
    import cupy as cp
    n = nbits(perm)
    D = 1 << n
    with cp.cuda.Device(device):
        X = cp.array([[0., 1.], [1., 0.]])
        Z = cp.array([[1., 0.], [0., -1.]])
        I = cp.eye(2)
        P = cp.ones((1, 1))
        for q in range(n):                       # qubit 0 least significant
            f = I
            if (a >> q) & 1 and (b >> q) & 1:
                f = X @ Z
            elif (a >> q) & 1:
                f = X
            elif (b >> q) & 1:
                f = Z
            P = cp.kron(f, P)
        U = cp.zeros((D, D))
        U[cp.asarray(perm), cp.arange(D)] = 1.0
        M = U.T @ (P @ U)
        del P, U
        ys = cp.arange(D)
        cs = cp.arange(D)
        V = M[ys[None, :] ^ cs[:, None], ys[None, :]]      # V[c, y] = M[y^c, y]
        del M
        h = 1
        while h < D:                                        # batched WHT over y
            w = V.reshape(D, -1, 2, h)
            t = w[:, :, 0, :].copy()
            w[:, :, 0, :] += w[:, :, 1, :]
            w[:, :, 1, :] = t - w[:, :, 1, :]
            del t
            h *= 2
        nzc, nzd = cp.nonzero(cp.abs(V) > 0.5)
        del V
        nzc, nzd = cp.asnumpy(nzc), cp.asnumpy(nzd)
        cp.get_default_memory_pool().free_all_blocks()
    out: dict[int, np.ndarray] = {}
    order = np.lexsort((nzd, nzc))
    nzc, nzd = nzc[order], nzd[order]
    for c in np.unique(nzc):
        lo, hi = np.searchsorted(nzc, [c, c + 1])
        out[int(c)] = nzd[lo:hi].astype(np.int64)
    return out


# ---------------------------------------------------------------- general finite abelian 2-groups
# An element of G = Z/n_0 x Z/n_1 x ... is stored as the mixed-radix integer
# sum_j y_j * prod_{i<j} n_i, so GF(2)^n is radix (2,)*n and a register of
# qubits whose top t bits form one Z/2^t digit is radix (2,)*low + (2^t,).
# Weyl operators W(a, chi) = T_a D_chi act as W|y> = chi(y)|y + a>, and
# U^dag W(a, chi) U = chi(pi(y)) |sigma(y)>, sigma = pi^-1 T_a pi. The
# coefficient on W(c, psi) is |G|^-1 sum_{y in D_c} chi(pi(y)) conj(psi(y)),
# D_c = {y : sigma(y) - y = c}; see C99 for G = GF(2)^n.

def _digits(y, radix):
    y = np.asarray(y, dtype=np.int64)
    out, place = [], 1
    for r in radix:
        out.append((y // place) % r)
        place *= r
    return out


def _undigits(ds, radix):
    y, place = 0, 1
    for d, r in zip(ds, radix):
        y = y + (np.asarray(d, dtype=np.int64) % r) * place
        place *= r
    return y


def g_add(x, a, radix):
    return _undigits([u + v for u, v in zip(_digits(x, radix), _digits(a, radix))], radix)


def g_sub(x, a, radix):
    return _undigits([u - v for u, v in zip(_digits(x, radix), _digits(a, radix))], radix)


def g_char(k, y, radix):
    """chi_k(y) = exp(2 pi i sum_j k_j y_j / n_j)."""
    ph = 0.0
    for kj, yj, r in zip(_digits(k, radix), _digits(y, radix), radix):
        ph = ph + (kj * yj) / r
    return np.exp(2j * np.pi * ph)


def g_ft_support(f, radix, tol=1e-9):
    """Nonzero-frequency mask of the group Fourier transform (any sign convention)."""
    F = np.fft.fftn(np.asarray(f).reshape(tuple(reversed(radix))))
    return np.abs(F).reshape(-1) > tol


def group_classes(perm, a, radix):
    inv = np.argsort(perm)
    ys = np.arange(len(perm), dtype=np.int64)
    sigma = inv[g_add(perm, a, radix)]
    return g_sub(sigma, ys, radix)


def group_terms(perm, a, k, radix, tol=1e-9):
    """(total term count, [(c, |D_c|, support)]) of U^dag W(a, chi_k) U in G's Weyl basis."""
    cvec = group_classes(perm, a, radix)
    s = g_char(k, perm, radix)
    total, classes = 0, []
    for c in np.unique(cvec):
        mask = cvec == c
        g = np.where(mask, s, 0)
        cnt = int(g_ft_support(g, radix, tol).sum())
        total += cnt
        classes.append((int(c), int(mask.sum()), cnt))
    return total, classes


def permutation_matrix_terms(sigma, radix, tol=1e-9):
    """Term count of the (unsigned) permutation matrix |sigma(y)><y| in G's Weyl basis."""
    ys = np.arange(len(sigma), dtype=np.int64)
    cvec = g_sub(sigma, ys, radix)
    total = 0
    R = 0
    for c in np.unique(cvec):
        total += int(g_ft_support((cvec == c).astype(float), radix, tol).sum())
        R += 1
    return total, R


def dense_group_terms(perm, radix, tol=1e-9):
    """{(a, k): term count} for all |G|^2 labels by explicit Weyl matrices and traces."""
    D = len(perm)
    ys = np.arange(D)
    U = np.zeros((D, D))
    U[perm, ys] = 1.0
    W = np.zeros((D * D, D, D), dtype=complex)
    for a in range(D):
        tgt = g_add(ys, a, radix)
        for k in range(D):
            W[a * D + k, tgt, ys] = g_char(k, ys, radix)
    Q = np.einsum("ab,jbc,cd->jad", U.T, W, U, optimize=True)
    ptm = (W.conj().reshape(D * D, -1) @ Q.reshape(D * D, -1).T) / D
    counts = (np.abs(ptm) > tol).sum(axis=0)
    return {(j // D, j % D): int(counts[j]) for j in range(D * D)}, ptm
