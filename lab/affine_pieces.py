"""Pull back a computational projector as a signed sum of affine-subspace
indicators (projectors of Z-type stabilizer codes).

Object: for a classical reversible circuit U and qubit j, the pulled-back
projector U^dag |1><1|_j U is the diagonal 0/1 function f(y) = bit j of U(y),
the same Boolean function whose Walsh sparsity C8 counts. Here f is kept as

    f = sum_k c_k 1[A_k],     A_k = {y : r . y = v_r for each constraint row r}

with integer c_k. Rules, for the self-inverse gates X, CNOT and Toffoli
(1[A] o g = 1[g(A)]):

  X(q)          flip the right-hand side of every row containing q    1 piece
  CNOT(c,t)     every row containing t also gets c                    1 piece
  Toffoli(a,b;t) with some row containing t:
      if y_a y_b is affine on A, substitute it                        1 piece
      otherwise branch on ONE control u (w the other):
          A and {y_u = 0}                         (product vanishes)
          A with y_t -> y_t + y_w, and {y_u = 1}  (product equals y_w)  2 pieces

`branch` selects u: "first" is op[1], "second" is op[2] of the logical
("toffoli", op[1], op[2], target) record. After every Toffoli an optional
sibling merge replaces c 1[B, row_i = 0] + c 1[B, row_i = 1] by c 1[B].

A key is the fully reduced echelon form as a sorted tuple of ints: bits 0..n-1
hold the row mask, bit n the right-hand side. This is exploratory storage;
its piece count is not a canonical invariant of f (the decomposition is not
unique) and not an allocation measurement.
"""
from __future__ import annotations

import numpy as np


def canon(rows, n):
    """Fully reduced echelon form of the constraint rows, or None if empty."""
    mask = (1 << n) - 1
    basis: dict[int, int] = {}
    for r in rows:
        for p, br in basis.items():
            if (r >> p) & 1:
                r ^= br
        m = r & mask
        if m == 0:
            if r:
                return None
            continue
        p = m.bit_length() - 1
        for q in list(basis):
            if (basis[q] >> p) & 1:
                basis[q] ^= r
        basis[p] = r
    return tuple(sorted(basis.values()))


def _constant(key, form, n):
    """Value of the linear form on A if constant there, else None."""
    r = form
    for row in key:
        p = (row & ((1 << n) - 1)).bit_length() - 1
        if (r >> p) & 1:
            r ^= row
    if r & ((1 << n) - 1):
        return None
    return (r >> n) & 1


def _toffoli(key, a, b, t, n, branch):
    rows = list(key)
    hit = [i for i, r in enumerate(rows) if (r >> t) & 1]
    if not hit:
        return [key]
    r0 = rows[hit[0]]
    rest = [r ^ r0 if (r >> t) & 1 else r
            for i, r in enumerate(rows) if i != hit[0]]
    ca, cb = _constant(key, 1 << a, n), _constant(key, 1 << b, n)
    if ca == 0 or cb == 0:
        return [key]
    if ca == 1 and cb == 1:
        return [canon(rest + [r0 ^ (1 << n)], n)]
    if ca == 1:
        return [canon(rest + [r0 ^ (1 << b)], n)]
    if cb == 1:
        return [canon(rest + [r0 ^ (1 << a)], n)]
    cab = _constant(key, (1 << a) | (1 << b), n)
    if cab == 1:
        return [key]
    if cab == 0:
        return [canon(rest + [r0 ^ (1 << a)], n)]
    u, w = (a, b) if branch == "first" else (b, a)
    out = [canon(rows + [1 << u], n),
           canon(rest + [r0 ^ (1 << w), (1 << u) | (1 << n)], n)]
    return [k for k in out if k is not None]


def _merge_siblings(terms, n):
    flip = 1 << n
    work = list(terms)
    while work:
        nxt = []
        for k in work:
            if k not in terms:
                continue
            c = terms[k]
            for i, r in enumerate(k):
                sib = tuple(sorted(k[:i] + (r ^ flip,) + k[i + 1:]))
                if sib != k and terms.get(sib) == c:
                    del terms[k], terms[sib]
                    u = canon(k[:i] + k[i + 1:], n)
                    v = terms.get(u, 0) + c
                    if v:
                        terms[u] = v
                        nxt.append(u)
                    else:
                        terms.pop(u, None)
                    break
        work = nxt
    return terms


def pullback(circuit, j: int, branch: str = "second", merge: bool = True,
             cap: int = 3_000_000):
    """Return (terms, peak) for the projector onto |1> at qubit j."""
    if branch not in ("first", "second"):
        raise ValueError("branch must be 'first' or 'second'")
    if not circuit.is_classical():
        raise ValueError("circuit is not a gate-by-gate permutation")
    n = circuit.n
    terms = {canon([(1 << j) | (1 << n)], n): 1}
    peak = 1
    for op in reversed(circuit.logical):
        new: dict = {}
        if op[0] == "x":
            for k, c in terms.items():
                k2 = canon([r ^ (1 << n) if (r >> op[1]) & 1 else r for r in k], n)
                new[k2] = new.get(k2, 0) + c
        elif op[0] == "cnot":
            ctl, tgt = op[1], op[2]
            for k, c in terms.items():
                k2 = canon([r ^ (1 << ctl) if (r >> tgt) & 1 else r for r in k], n)
                new[k2] = new.get(k2, 0) + c
        else:
            for k, c in terms.items():
                for k2 in _toffoli(k, op[1], op[2], op[3], n, branch):
                    new[k2] = new.get(k2, 0) + c
        terms = {k: c for k, c in new.items() if c}
        peak = max(peak, len(terms))
        if merge and op[0] == "toffoli":
            terms = _merge_siblings(terms, n)
        if len(terms) > cap:
            raise MemoryError(f"{len(terms)} pieces exceed cap {cap}")
    return terms, peak


def evaluate(terms, n: int) -> np.ndarray:
    """Dense integer value of the signed sum on every basis label (checking only)."""
    ys = np.arange(1 << n, dtype=np.int64)
    out = np.zeros(1 << n, dtype=np.int64)
    for k, c in terms.items():
        ok = np.ones(1 << n, dtype=bool)
        for r in k:
            m, v = r & ((1 << n) - 1), (r >> n) & 1
            par = np.zeros(1 << n, dtype=np.int64)
            while m:
                q = (m & -m).bit_length() - 1
                par ^= (ys >> q) & 1
                m &= m - 1
            ok &= par == v
        out += c * ok
    return out
