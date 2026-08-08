"""Modular-arithmetic bookkeeping used across the 2-adic experiments."""
from __future__ import annotations
from math import gcd

import numpy as np


def order(a: int, N: int) -> int:
    """Multiplicative order of a mod N."""
    if gcd(a, N) != 1:
        raise ValueError(f"a={a} not coprime to N={N}")
    r, v = 1, a % N
    while v != 1:
        v = (v * a) % N
        r += 1
    return r


def v2_split(r: int) -> tuple[int, int]:
    """r = beta * 2^alpha with beta odd; returns (alpha, beta)."""
    alpha = 0
    while r % 2 == 0:
        r //= 2
        alpha += 1
    return alpha, r


def ord2(beta: int) -> int:
    """Multiplicative order of 2 mod beta (0 for beta = 1)."""
    if beta == 1:
        return 0
    return order(2, beta)


def carmichael(N: int) -> int:
    """Carmichael function lambda(N)."""
    def lam_pk(p, k):
        if p == 2 and k >= 3:
            return 1 << (k - 2)
        return (p - 1) * p ** (k - 1)

    out, n = 1, N
    p = 2
    while p * p <= n:
        if n % p == 0:
            k = 0
            while n % p == 0:
                n //= p
                k += 1
            v = lam_pk(p, k)
            out = out * v // gcd(out, v)
        p += 1
    if n > 1:
        v = lam_pk(n, 1)
        out = out * v // gcd(out, v)
    return out


def bit_table(N: int, a: int, bit: int) -> np.ndarray:
    """h[c] = bit_j(a^c mod N), length r -- the real modexp bit function's
    period table."""
    r = order(a, N)
    h = np.empty(r, dtype=np.int64)
    x = 1
    for c in range(r):
        h[c] = (x >> bit) & 1
        x = (x * a) % N
    return h


def tile(h: np.ndarray, t: int) -> np.ndarray:
    """g(e) = h[e mod r] over e in [0, 2^t) -- the idealised bit function."""
    r = h.size
    return np.tile(h, (1 << t) // r + 1)[: 1 << t]
