"""Close step 6's scope caveat: redo it with the REAL modexp bit function.

Step 6 swept r freely using g(e) = h[e mod r] with a RANDOM table h. That left
open whether the fine structure within beta>1 is a property of the arithmetic
(e mod r sitting inside a binary Walsh basis) or an artifact of h being random.

The real function has exactly the same form, g(e) = h[e mod r], but with the
specific table h[c] = bit_j(a^c mod N). So the question is sharp: does replacing
random h by the real h change the structure?

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  beta = 1: sparsity constant in t and <= 2^alpha. This is proved at
      function level, so it MUST hold -- a failure means the pipeline is wrong.
  P2  beta > 1: the qualitative structure (period-ord2(beta) oscillation plus
      slow drift) should TRANSFER, because it arises from how e mod r sits in
      the GF(2)^t Walsh basis, which is independent of h. Only the specific
      density values should differ.
  P3  If instead the real function shows markedly DIFFERENT structure -- much
      sparser, or clean periodicity without drift -- then modexp has special
      Walsh structure beyond "depends on e mod r", which would be a finding in
      its own right and would need saying in Paper B.
"""
from __future__ import annotations
import numpy as np
from math import gcd
import walsh

TOL = 1e-12


def order(a, N):
    r, v = 1, a % N
    while v != 1:
        v = (v * a) % N
        r += 1
    return r


def split(r):
    a, rr = 0, r
    while rr % 2 == 0:
        rr //= 2
        a += 1
    return a, rr


def ord2(b):
    if b == 1:
        return 0
    k, v = 1, 2 % b
    while v != 1:
        v = (v * 2) % b
        k += 1
    return k


def real_table(N, a, r, bit):
    """h[c] = bit_j(a^c mod N)."""
    h = np.empty(r, dtype=np.int64)
    x = 1
    for c in range(r):
        h[c] = (x >> bit) & 1
        x = (x * a) % N
    return h


def density(h, r, t):
    vals = np.tile(h, (1 << t) // r + 1)[: 1 << t]
    c = walsh.wht(np.where(vals == 1, -1.0, 1.0)) / (1 << t)
    S = int((np.abs(c) > TOL).sum())
    return S / (1 << t), S


TS = list(range(11, 23))

print("=" * 92)
print("REAL modexp bit function g(e) = bit_j(a^e mod N), fixed (N,a,j), sweeping t")
print("=" * 92)

CASES = [(15, 7, 0), (17, 3, 0),          # beta = 1 (r = 4, 16)
         (7, 3, 0), (7, 2, 0),            # beta = 3 (r = 6, 3)
         (13, 2, 0), (11, 2, 0),          # beta = 3 (r=12), beta = 5 (r=10)
         (19, 2, 0), (23, 5, 0)]          # beta = 9 (r=18), beta = 11 (r=22)

for N, a, bit in CASES:
    if gcd(a, N) != 1:
        continue
    r = order(a, N)
    al, beta = split(r)
    o = ord2(beta)
    h = real_table(N, a, r, bit)
    if not (0 < int(h.sum()) < r):
        print(f"\n  N={N} a={a} bit={bit}: r={r}, table is constant — skipped")
        continue
    ds = [density(h, r, t)[0] for t in TS]
    Ss = [density(h, r, t)[1] for t in TS]
    print(f"\n  N={N:3d} a={a:2d} bit={bit}  r={r:3d}  alpha={al} beta={beta:2d} "
          f"ord2(beta)={o}")
    print("    t        : " + " ".join(f"{t:6d}" for t in TS))
    print("    density  : " + " ".join(f"{d:6.3f}" for d in ds))
    if beta == 1:
        print("    sparsity : " + " ".join(f"{s:6d}" for s in Ss))
        const = len(set(Ss)) == 1
        print(f"    -> P1 sparsity constant: {const}  (bound 2^alpha = {1<<al})")
    else:
        groups = {}
        for t, d in zip(TS, ds):
            groups.setdefault(t % o, []).append(d)
        spreads = {k: max(v) - min(v) for k, v in groups.items()}
        exact = all(s < 1e-9 for s in spreads.values())
        drift = max(spreads.values())
        print(f"    -> grouped by t mod {o}: exact periodicity {exact}, "
              f"max within-class spread {drift:.4f}")

print()
print("=" * 92)
print("COMPARISON: same r, real table vs random table")
print("=" * 92)
print("  If P2 holds, both show the same qualitative shape (period + drift) and")
print("  differ only in the numbers.\n")
rng = np.random.default_rng(3)
for N, a in ((7, 3), (13, 2), (11, 2)):
    r = order(a, N)
    al, beta = split(r)
    o = ord2(beta)
    hr = real_table(N, a, r, 0)
    hx = rng.integers(0, 2, size=r)
    while not (0 < int(hx.sum()) < r):
        hx = rng.integers(0, 2, size=r)
    dr = [density(hr, r, t)[0] for t in TS]
    dx = [density(hx, r, t)[0] for t in TS]
    print(f"  r={r:3d} (beta={beta}, ord2={o})")
    print("    real  : " + " ".join(f"{d:6.3f}" for d in dr))
    print("    random: " + " ".join(f"{d:6.3f}" for d in dx))
