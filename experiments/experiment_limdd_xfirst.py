"""How many Pauli-LIMDD nodes does the post-modexp state need after the x register, in x-first orders?

Context: C110 bounds exponent-first orders.  In any order that reads every x qubit before
every exponent qubit (scratch anywhere), the sub-state after x = a^l mod N is the
indicator of P_l = {e < 2^t : e = l mod r}, and Pauli-LIM equivalence of such 0/1 states is
equality of supports up to an XOR shift (C110 Lemma 1's argument).  Write r = beta * 2^alpha
(beta odd), t' = t - alpha, L = bitlen(beta), g = gcd of the nonzero set-bit positions of beta,
Q_k = {e < 2^t' : e = k mod beta}.  Derivation (draft claim C111):

  R   P_l = {l mod 2^alpha} x Q_k (low alpha bits | high t' bits), k = floor(l / 2^alpha), and
      P_l ~ P_l' iff Q_k ~ Q_k'.  So classes(P) = classes(Q).
  U   Complement e -> 2^t' - 1 - e maps Q_k onto Q_{(2^t' - 1 - k) mod beta}: an involution on
      Z_beta with one fixed point, so classes <= (beta + 1)/2.
  X   If Q_k ^ s = Q_k' then (e & s) mod beta is constant on Q_k.  For every window m with an
      e1 in Q_k whose bits vanish on m + supp(beta), e2 = e1 | (beta << m) is in Q_k, so s holds
      all or none of m + supp(beta).  For t' >= 3L - 1 every window m in [0, t' - L] has such an
      e1, and Fine-Wilf makes s g-periodic.  If g = 1, s is 0 or all ones: classes = (beta+1)/2.

PREDICTIONS, WRITTEN BEFORE THIS SCRIPT RAN (a scratch run of the same brute force had
already produced X2's values for the betas below; see note LM):

  P0  (object) for every unit base a of every odd N in [5, 32), t = bitlen(r) + 1: the
      sub-states after all x qubits, read off build()'s output on |e>|0> (x register after
      the gates; scratch checked clean), are exactly the indicators of {e : e = l mod r}
      for y = a^l, and their XOR-shift class count equals the count for the P_l sets built
      from r alone.  (Added before this script's second run, after review
      Va10c9887fa3f4116 of C110 found the analogous step untested there.)
  P1  (R) for r in {6, 12, 20, 24, 28, 40, 44, 56}, t = alpha + t' with t' = 3L - 1 capped so
      2^t <= 2^15: classes of the P_l (t-bit brute force) equal classes of the Q_k (t'-bit).
  P2  (X, g = 1) for every odd beta in [3, 64) with g = 1 and 3L - 1 <= 17: at t' = 3L - 1 the
      only shifts s with Q_k ^ s = Q_k' for some k are 0 and 2^t' - 1, and classes = (beta+1)/2.
  P3  (U) classes <= (beta + 1)/2 at every (beta, t') tested, t' from L + 1 to 3L - 1.
  C1  must fail: "classes = (beta+1)/2 at every t' >= L + 1" (i.e. that the threshold is not
      needed) must be false somewhere (beta = 31 at t' = 6 gave 6 < 16 in note LM's scan).
  C2  must fail: the exact form without the complement merge, "classes = beta", must be false
      at some g = 1 point of P2.
  OPEN (logged, not graded): g > 1 betas at t' = 3L - 1.  The proof gives only g-periodic
      shifts (classes >= beta / 2^g); the scratch run found only trivial shifts there too.

Run:  uv run python -m experiments.experiment_limdd_xfirst     (from research/)
"""
from __future__ import annotations

import math
from functools import reduce

import numpy as np

from lab import Experiment

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)
exp.predict("P0", "x-first sub-states read off build() are the P_l indicators; class counts agree")
exp.predict("P1", "classes(P_l, t bits) == classes(Q_k, t' bits)")
exp.predict("P2", "g = 1, t' = 3L-1: only shifts 0 and all-ones; classes == (beta+1)/2")
exp.predict("P3", "classes <= (beta+1)/2 at every tested (beta, t')")
exp.must_fail("C1", "classes == (beta+1)/2 for all t' >= L+1 is false somewhere")
exp.must_fail("C2", "classes == beta is false at some g = 1 point")


def gpos(beta):
    return reduce(math.gcd, [i for i in range(1, beta.bit_length()) if beta >> i & 1], 0)


def shifts_for(mod, width, k):
    """(s, k2) with {e < 2^width : e = k mod mod} ^ s == the class of k2 (brute force)."""
    M = 1 << width
    Q = np.arange(k, M, mod)
    cand = np.arange(M)
    for chunk in (Q[:16], Q):
        res = (chunk[None, :] ^ cand[:, None]) % mod
        cand = cand[(res == res[:, :1]).all(axis=1)]
    out = []
    for s in cand:
        k2 = int((Q[0] ^ s) % mod)
        if len(range(k2, M, mod)) == len(Q):
            out.append((int(s), k2))
    return out


def classes(mod, width):
    parent = list(range(mod))

    def root(x):
        while parent[x] != x:
            x = parent[x]
        return x
    shifts = set()
    for k in range(mod):
        for s, k2 in shifts_for(mod, width, k):
            shifts.add(s)
            parent[root(k)] = root(k2)
    return len({root(x) for x in range(mod)}), shifts


def split(r):
    alpha = 0
    while r % 2 == 0:
        r //= 2
        alpha += 1
    return r, alpha


rows = []

exp.section("P0  sub-states read off the circuit")
from lab.modarith import order as mult_order
from toffoli_arith import ToffoliModExp
from walsh import classical_images


def support_classes(sets, width):
    canon = set()
    for S in sets:
        canon.add(min(tuple(sorted(x ^ x0 for x in S)) for x0 in S))
    return len(canon)


bad0 = n0 = 0
bases0 = set()
for N in range(5, 32, 2):
    for a in range(2, N):
        if math.gcd(a, N) != 1:
            continue
        r = mult_order(a, N)
        t = r.bit_length() + 1
        me = ToffoliModExp(N, a, n_exp=t)
        ins = np.array([sum(((e >> i) & 1) << q for i, q in enumerate(me.exp)) for e in range(1 << t)],
                       dtype=np.int64)
        outs = classical_images(me.build(), ins)
        xmask = sum(1 << q for q in me.x)
        emask = sum(1 << q for q in me.exp)
        groups = {}
        for e, yv in enumerate(outs.tolist()):
            assert yv & ~xmask & ~emask == 0, "scratch not clean"
            y = sum(((yv >> q) & 1) << i for i, q in enumerate(me.x))
            groups.setdefault(y, set()).add(e)
        state_sets = [tuple(sorted(v)) for v in groups.values()]
        derived = [tuple(range(l, 1 << t, r)) for l in range(r)]
        ok = (sorted(state_sets) == sorted(derived)
              and support_classes(state_sets, t) == classes(r, t)[0])
        bad0 += not ok
        n0 += 1
        bases0.add(a)
exp.check("P0", n0 > 0 and bad0 == 0, f"{n0} (N, a) circuits, {len(bases0)} distinct bases, {bad0} mismatches")

exp.section("P1  reduction to the odd part")
bad = 0
for r in (6, 12, 20, 24, 28, 40, 44, 56):
    beta, alpha = split(r)
    L = beta.bit_length()
    tp = min(3 * L - 1, 15 - alpha)
    cP, _ = classes(r, alpha + tp)
    cQ, _ = classes(beta, tp)
    exp.log(f"r={r} beta={beta} alpha={alpha} t'={tp}: classes P={cP}, Q={cQ}")
    bad += cP != cQ
exp.check("P1", bad == 0, f"{bad} mismatches")

exp.section("P2/P3/C1/C2  odd beta")
bad2 = bad3 = 0
n2 = 0
c1_hit = c2_hit = False
for beta in range(3, 64, 2):
    L = beta.bit_length()
    g = gpos(beta)
    if 3 * L - 1 > 17:
        continue
    for tp in range(L + 1, 3 * L):
        c, shifts = classes(beta, tp)
        rows.append(dict(beta=beta, L=L, g=g, tp=tp, classes=c, nshifts=len(shifts)))
        bad3 += c > (beta + 1) // 2
        if c != (beta + 1) // 2:
            c1_hit = True
        if tp == 3 * L - 1:
            if g == 1:
                n2 += 1
                ok = shifts <= {0, (1 << tp) - 1} and c == (beta + 1) // 2
                bad2 += not ok
                c2_hit |= c != beta
            else:
                exp.log(f"OPEN g={g}: beta={beta} t'={tp} classes={c} (beta+1)/2={(beta + 1) // 2} "
                        f"shifts={sorted(shifts)}")
exp.check("P2", n2 > 0 and bad2 == 0, f"{n2} g = 1 betas at t' = 3L-1, {bad2} violations")
exp.check("P3", bad3 == 0, f"{len(rows)} (beta, t') points, {bad3} above (beta+1)/2")
below = [(x['beta'], x['tp'], x['classes']) for x in rows if x['classes'] != (x['beta'] + 1) // 2]
exp.log(f"points below (beta+1)/2 (beta, t', classes): {below[:12]}{' ...' if len(below) > 12 else ''}"
        f" ({len(below)} total)")
exp.fail_check("C1", c1_hit, "the threshold matters: some t' >= L+1 falls short of (beta+1)/2")
exp.fail_check("C2", c2_hit, "classes == beta fails (the complement merge is real)")

import os
exp.finish(report_path=os.environ.get("LIMDD_REPORT", "out/limdd_xfirst/report.json"), rows=rows)
