"""Step 6, careful redo. The first pass had a methodology bug and a real result.

BUG: for small r the random bit table h has few possibilities (r=3 gives 8, of
which 2 are constant), so the mean density was contaminated by degenerate
tables. That is why beta=3, alpha=0 appeared to have density ~0.000. Fix by
excluding constant tables, using balanced tables, and reporting spread rather
than only the mean.

REAL RESULT that survived and REFUTES the prediction: for beta > 1 the density
does NOT converge to 1. It PLATEAUS at a value that is constant in t and
depends on r -- 0.667 for r=6, 0.833 for r=12, 0.875 for r=24. So sparsity is
Theta(2^t) either way, but the CONSTANT differs. P1 was wrong in its second
half: alpha does not vanish asymptotically, it sets the plateau.

This pass measures the plateau precisely and looks for its closed form.
"""
from __future__ import annotations
import numpy as np
from math import gcd
import walsh

TOL = 1e-12


def density_stats(r, t, trials=12, balanced=True, seed=0):
    """Density of g(e)=h[e mod r], excluding constant tables."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(trials):
        for _try in range(50):
            if balanced and r >= 2:
                h = np.zeros(r, dtype=np.int64)
                h[rng.choice(r, size=r // 2, replace=False)] = 1
            else:
                h = rng.integers(0, 2, size=r)
            if 0 < int(h.sum()) < r:      # reject constants
                break
        vals = np.tile(h, (1 << t) // r + 1)[: 1 << t]
        c = walsh.wht(np.where(vals == 1, -1.0, 1.0)) / (1 << t)
        out.append(int((np.abs(c) > TOL).sum()) / (1 << t))
    a = np.array(out)
    return float(a.mean()), float(a.min()), float(a.max())


print("=" * 88)
print("1. CLEAN density (balanced, non-constant tables), mean [min-max] over 12 trials")
print("=" * 88)
TS = (12, 16, 20)
print(f"  {'r':>5} {'alpha':>6} {'beta':>5} | " +
      "  ".join(f"{'t=' + str(t):^20}" for t in TS))
for r in (2, 3, 4, 5, 6, 8, 9, 12, 16, 24, 48, 5 * 2, 5 * 4, 9 * 2, 9 * 4):
    al = 0
    rr = r
    while rr % 2 == 0:
        rr //= 2
        al += 1
    beta = rr
    cells = []
    for t in TS:
        m, lo, hi = density_stats(r, t)
        cells.append(f"{m:.4f} [{lo:.3f}-{hi:.3f}]")
    print(f"  {r:5d} {al:6d} {beta:5d} | " + "  ".join(f"{c:^20}" for c in cells),
          flush=True)

print("""
  If a row is flat across t, the sparsity is Theta(2^t) with that density as the
  constant. beta=1 rows should instead fall by 16x per column (constant
  sparsity).""")

print()
print("=" * 88)
print("2. THE PLATEAU: does it have a closed form?")
print("=" * 88)
print("  Candidate: density = 1 - (something)/r.  Compare against 1/2^alpha and")
print("  against multiplicative order of 2 mod beta.\n")
print(f"  {'r':>5} {'alpha':>6} {'beta':>5} {'density':>9} {'1-1/beta':>9} "
      f"{'1-2^-al':>9} {'ord2(beta)':>11}")


def ord2(b):
    if b == 1:
        return 0
    k, v = 1, 2 % b
    while v != 1:
        v = (v * 2) % b
        k += 1
    return k


for r in (3, 5, 6, 9, 10, 12, 18, 20, 24, 36, 48, 72):
    al = 0
    rr = r
    while rr % 2 == 0:
        rr //= 2
        al += 1
    beta = rr
    m, _, _ = density_stats(r, 20)
    print(f"  {r:5d} {al:6d} {beta:5d} {m:9.4f} {1-1/beta:9.4f} "
          f"{1-2**-al:9.4f} {ord2(beta):11d}", flush=True)

print("""
  Looking for which column tracks the measured density. ord2(beta) is the
  multiplicative order of 2 modulo the odd part -- the natural candidate, since
  it controls how the doubling map e -> 2e closes up modulo beta, and the Walsh
  basis is built from exactly that doubling structure.""")
