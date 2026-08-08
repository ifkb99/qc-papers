"""TODO step 6 -- is there an INTERMEDIATE 2-adic law, or is the split binary?

Motivation. N=323 has r = 144 = 16*9, i.e. alpha=4 with an odd part beta=9, and
showed partial sparsity at small t (density 0.981 at t=16) washing out to 1.000
by t=24. That looked like it might be a quantitative law in alpha versus t
rather than the binary power-of-two / not split C15 asserts.

PREDICTIONS, WRITTEN BEFORE MEASURING.

Derivation. g depends on e mod r, and by CRT e mod r <-> (e mod 2^alpha,
e mod beta).
  * e mod 2^alpha is exactly the low alpha bits -- a genuine GF(2)-affine
    function of the input bits, and Walsh-friendly.
  * e mod beta for odd beta > 1: writing e = u + 2^alpha * v, this is
    (u + 2^alpha v) mod beta, which depends on ALL bits **including the low
    ones**. Its level sets are arithmetic progressions of odd stride, which are
    not unions of GF(2)-affine subspaces.

So when beta > 1 the mod-beta part RE-ENTANGLES the low bits and the 2^alpha
structure never separates out. Therefore:

  P1  beta > 1  =>  density -> 1 as t grows, for EVERY alpha. Alpha buys no
      asymptotic reduction at all.
  P2  beta = 1  =>  sparsity <= 2^alpha, constant in t. (control, must hold)
  P3  Alpha may affect the finite-size APPROACH RATE to density 1, but not the
      limit. The N=323 observation should be exactly this.
  P4  MUST FAIL: if a genuine intermediate law existed, density would plateau
      at some value strictly between 0 and 1 depending on alpha. Predict no such
      plateau appears.

If P1-P4 land, step 6's premise is answered NEGATIVELY and the dichotomy is
genuinely binary -- which SHARPENS C15 rather than complicating it.

Method. Test the arithmetic in isolation: g(e) = h[e mod r] for a fixed
pseudorandom bit table h of length r. This strips away everything specific to
modular exponentiation and measures only how "mod r" sits inside a binary Walsh
basis, letting r = beta * 2^alpha be swept freely instead of hunting for moduli
with prescribed orders.
"""
from __future__ import annotations
import numpy as np
import walsh

TOL = 1e-12


def density(r, t, seed=0, trials=3):
    """Mean Walsh density of g(e) = h[e mod r] over e in [0, 2^t)."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(trials):
        h = rng.integers(0, 2, size=r)
        vals = np.tile(h, (1 << t) // r + 1)[: 1 << t]
        c = walsh.wht(np.where(vals == 1, -1.0, 1.0)) / (1 << t)
        out.append(int((np.abs(c) > TOL).sum()) / (1 << t))
    return float(np.mean(out))


print("=" * 84)
print("P1/P2/P4  density of g(e) = h[e mod r], r = beta * 2^alpha")
print("=" * 84)
TS = (10, 12, 14, 16, 18, 20)
print(f"  {'beta':>5} {'alpha':>6} {'r':>6} | " + " ".join(f"t={t:<2d}" for t in TS))
for beta in (1, 3, 5, 9):
    for alpha in (0, 1, 2, 3, 4):
        r = beta * (1 << alpha)
        if r < 2:
            continue
        ds = [density(r, t) for t in TS]
        print(f"  {beta:5d} {alpha:6d} {r:6d} | " +
              " ".join(f"{d:5.3f}" for d in ds), flush=True)
    print()

print("""  beta = 1 rows: density must fall as 2^alpha / 2^t  (sparsity CONSTANT).
  beta > 1 rows: if P1 holds, every row climbs to 1.000 whatever alpha is, and
  alpha shows up only in how fast it gets there (P3).""")

print()
print("=" * 84)
print("P2 check: is the beta=1 sparsity exactly 2^alpha?")
print("=" * 84)
print(f"  {'alpha':>6} {'r':>4} {'t':>4} {'sparsity':>10} {'2^alpha':>9} {'equal':>7}")
for alpha in (1, 2, 3, 4, 5):
    r = 1 << alpha
    for t in (12, 18):
        rng = np.random.default_rng(0)
        h = rng.integers(0, 2, size=r)
        vals = np.tile(h, (1 << t) // r + 1)[: 1 << t]
        c = walsh.wht(np.where(vals == 1, -1.0, 1.0)) / (1 << t)
        S = int((np.abs(c) > TOL).sum())
        print(f"  {alpha:6d} {r:4d} {t:4d} {S:10d} {1<<alpha:9d} "
              f"{str(S <= (1 << alpha)):>7}", flush=True)

print()
print("=" * 84)
print("P3  Does alpha affect the approach rate at fixed odd part?")
print("=" * 84)
print("  Same beta=9, varying alpha. If alpha mattered asymptotically the rows")
print("  would separate and stay separated.\n")
print(f"  {'r':>6} {'alpha':>6} | " + " ".join(f"t={t:<2d}" for t in TS))
for alpha in (0, 1, 2, 3, 4):
    r = 9 * (1 << alpha)
    ds = [density(r, t) for t in TS]
    print(f"  {r:6d} {alpha:6d} | " + " ".join(f"{d:5.3f}" for d in ds), flush=True)

print("""
  Reading: a genuine intermediate law would show rows plateauing at distinct
  values set by alpha. Convergence of every row to 1.000, with alpha only
  shifting the curve rightwards, means the dichotomy is binary and the N=323
  observation was a finite-size effect.""")
