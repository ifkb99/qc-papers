"""Pin down the constancy THRESHOLD in n_exp.

experiment_c15c.py showed that for r=4 (alpha=2) the support is NOT constant
from n_exp=2: it grows 15493 -> 32143 between n_exp 2 and 4, then locks. For
r=2 (alpha=1) it is constant from n_exp=2 onward.

Hypothesis: the u_a block for exponent bit i multiplies by a^(2^i) mod N, which
equals 1 exactly when 2^i is a multiple of r, i.e. when i >= alpha. So blocks
i < alpha are non-trivial and blocks i >= alpha are identity-on-the-valid-
subspace. The support locks once at least ONE identity block is present, i.e.
for n_exp > alpha -- strictly greater, not >=.

Predictions:
  alpha=1 (r=2): constant for n_exp >= 2
  alpha=2 (r=4): constant for n_exp >= 3, and n_exp=2 differs
  alpha=3 (r=8): constant for n_exp >= 4, n_exp=2 and 3 differ

RESOLVED 2026-08-08 for alpha = 3 and 4. This file skips q > 22, which is
exactly where the alpha=3 and alpha=4 lock points live, so it could only
report "still growing, as predicted". `experiment_c21_onset.py` runs them on
the GPU (q up to 27) and both lock at exactly n_exp = alpha + 1, at two moduli
for alpha = 3. See NOTES.md SS OS1.
"""
from __future__ import annotations
import numpy as np, math, time
from toffoli_arith import ToffoliModExp
import walsh


def order(a, N):
    r, v = 1, a % N
    while v != 1:
        v = (v * a) % N; r += 1
    return r


def v2(r):
    k = 0
    while r % 2 == 0:
        r //= 2; k += 1
    return k


print("Threshold test: support vs n_exp, stepping by 1.\n")
print(f"  {'N':>4} {'a':>3} {'r':>3} {'alpha':>6} {'n_exp':>6} {'q':>4} "
      f"{'|support|':>10} {'locked?':>8} {'predict lock':>13}")
cases = [(5, 4), (5, 2), (17, 2), (17, 3)]
for N, a in cases:
    if math.gcd(a, N) != 1:
        continue
    r = order(a, N); al = v2(r)
    if (r >> al) != 1:
        note = f"  (beta={r>>al}, not free)"
    else:
        note = ""
    prev = None
    for ne in (1, 2, 3, 4, 5):
        me = ToffoliModExp(N=N, a=a, n_exp=ne)
        qc = me.build()
        if qc.n > 22:
            print(f"  {N:4d} {a:3d} {r:3d} {al:6d} {ne:6d} {qc.n:4d}  "
                  f"skipped (q>22)")
            continue
        c = walsh.pullback_coefficients(qc, me.x[0])
        sp = int((np.abs(c) > 1e-12).sum())
        lock = "-" if prev is None else ("LOCKED" if sp == prev else "grew")
        pred = f"n_exp >= {al+1}" if (r >> al) == 1 else "never"
        print(f"  {N:4d} {a:3d} {r:3d} {al:6d} {ne:6d} {qc.n:4d} {sp:10d} "
              f"{lock:>8} {pred:>13}{note}", flush=True)
        prev = sp
        del c
    print()

print("""  'LOCKED' means identical to the previous row. The prediction is that the
  first LOCKED row appears at n_exp = alpha + 1, because that is the first
  width at which an identity u_a block (a^(2^i) = 1) is present.""")
