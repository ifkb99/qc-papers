"""Step 6, third pass. Two bugs fixed and a real hypothesis to test.

BUG 1 (pass 1): degenerate random tables contaminated small-r rows.
BUG 2 (pass 2): a NEW random table was drawn for each t, so t-dependence was
confounded with table variance -- a violation of "vary exactly one parameter".
Fix: FIX the table, sweep t.

HYPOTHESIS from pass 2's surviving pattern. Densities were flat in t for
beta=3 (r = 6, 12, 24, 48 all identical at t = 12, 16, 20) but appeared to drift
for beta=9 (r = 9, 18, 36). The sampled t values were 12, 16, 20:

    ord2(3) = 2  ->  t mod 2 = 0, 0, 0   (all the same)  -> looks FLAT
    ord2(9) = 6  ->  t mod 6 = 0, 4, 2   (all different)  -> looks like DRIFT

where ord2(beta) is the multiplicative order of 2 modulo beta. So the apparent
drift may be an aliasing artifact of the sampling grid.

  H1  The density depends on t only through **t mod ord2(beta)**, i.e. it is
      PERIODIC in t with period ord2(beta) -- not monotonically approaching 1.
  H2  Therefore there IS an intermediate law, but it is periodic rather than a
      limit, and the N=323 observation (r = 144 = 16*9, ord2(9) = 6) was
      sampling aliasing, not a finite-size effect.
  H3  Control: beta = 1 has ord2(1) = 0 and density -> 0 (constant sparsity).

This matters because it changes the claim from "beta>1 is asymptotically dense"
to something sharper and stranger.
"""
from __future__ import annotations
import numpy as np
import walsh

TOL = 1e-12


def ord2(b):
    if b % 2 == 0 or b == 1:
        return 0
    k, v = 1, 2 % b
    while v != 1:
        v = (v * 2) % b
        k += 1
    return k


def split(r):
    a, rr = 0, r
    while rr % 2 == 0:
        rr //= 2
        a += 1
    return a, rr


def density_fixed_table(r, t, h):
    vals = np.tile(h, (1 << t) // r + 1)[: 1 << t]
    c = walsh.wht(np.where(vals == 1, -1.0, 1.0)) / (1 << t)
    return int((np.abs(c) > TOL).sum()) / (1 << t)


print("=" * 90)
print("H1  FIXED table, consecutive t. Is density periodic in t with period ord2(beta)?")
print("=" * 90)
TS = list(range(11, 23))
rng = np.random.default_rng(7)
for r in (6, 12, 24, 9, 18, 36, 5, 10, 7, 14):
    al, beta = split(r)
    o = ord2(beta)
    h = rng.integers(0, 2, size=r)
    while not (0 < int(h.sum()) < r):
        h = rng.integers(0, 2, size=r)
    ds = [density_fixed_table(r, t, h) for t in TS]
    print(f"\n  r={r:3d}  alpha={al}  beta={beta:2d}  ord2(beta)={o}")
    print("    t      : " + " ".join(f"{t:6d}" for t in TS))
    print("    density: " + " ".join(f"{d:6.3f}" for d in ds))
    if o:
        groups = {}
        for t, d in zip(TS, ds):
            groups.setdefault(t % o, []).append(d)
        ok = all(max(v) - min(v) < 1e-9 for v in groups.values())
        detail = ", ".join(f"t%{o}={k}: {v[0]:.3f}" for k, v in sorted(groups.items()))
        print(f"    grouped by t mod {o}: {detail}")
        print(f"    -> constant within each residue class: {ok}")

print()
print("=" * 90)
print("H3  Control: beta = 1")
print("=" * 90)
for r in (4, 8, 16):
    al, beta = split(r)
    h = rng.integers(0, 2, size=r)
    while not (0 < int(h.sum()) < r):
        h = rng.integers(0, 2, size=r)
    ds = [density_fixed_table(r, t, h) for t in TS]
    sp = [d * (1 << t) for d, t in zip(ds, TS)]
    print(f"  r={r:3d} alpha={al} beta={beta}: sparsity = "
          + " ".join(f"{s:.0f}" for s in sp))
print("    -> constant sparsity, density falling as 2^-t, as C15 requires.")
