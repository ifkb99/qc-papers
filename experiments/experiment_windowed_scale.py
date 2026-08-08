"""TODO 12e -- does C38's 2-periodicity in the tail-window count survive K = 4, 5?

C38 says that in the table-lookup (Gidney-style) modexp the tail applies a fixed
involution W unconditionally, so the circuit is W^K . P_live and |support| is
**2-periodic** in K rather than constant: 28078 / 62680 / 28078 / 62680 at
K = 0..3. Four points is two periods, which is the minimum that can be called a
period at all, and `experiment_windowed.py` stopped there because K = 4 is
q = 27 and K = 5 is q = 29. Both are now minutes on the GPU.

This is worth doing because the alternative reading is not silly: a slow drift
with an even/odd wobble on top would look exactly like 2-periodicity over four
points, and that is precisely the failure mode that produced the retracted
"intermediate 2-adic law" (NOTES §I -- three points of an oscillation read as a
trend). Two more periods either kills that reading or does not.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  |support| at K = 4 equals the K = 0 value exactly, and K = 5 equals the
      K = 1 value exactly, at both beta = 1 moduli. Derived from C38's
      mechanism (W^K with W^2 = id), so it is a prediction about *exact
      equality*, not about a small drift.

  P2  Sets, not counts. If W^2 = id really makes K = 4 the same permutation as
      K = 0 on a wider register, then the K = 4 support restricted to the
      K = 0 coordinates must equal the K = 0 support **as a set**.
      `experiment_windowed.py` checked this for K = 0 vs 2; the mechanism says
      it holds for every even gap, so K = 0 vs 4 and K = 1 vs 5 must too.

  P3  C37's dead tail persists: the live exponent bits stay exactly
      {0 .. alpha-1} at K = 4 and 5, i.e. no tail bit ever enters the support,
      however long the tail gets.

  C1  MUST FAIL -- the beta > 1 lookup control (N = 7, a = 3, r = 6) has no
      tail window, so nothing is dead and |support| must GROW at K = 4 and 5
      rather than repeat. Without this, P1 could be flatness of the
      measurement rather than of the circuit.

OUTCOME (2026-08-08): 4/4. Both beta = 1 moduli repeat exactly over three full
periods -- [28078, 62680, 28078, 62680, 28078, 62680] and [13563, 62570,
13563, 62570, 13563, 62570] -- with the K=0/K=4 and K=1/K=5 supports identical
as SETS (4/4) and C37's live exponent bits still exactly {0..alpha-1}. The
beta = 3 control grows 84,827 -> 266,296,901 over the same range, ~4x per
window. The "slow drift with an even/odd wobble" reading is dead.

Run:  LAB_GPU=1 uv run python -m experiments.experiment_windowed_scale
      K = 5 is q = 29; cold cost is tens of minutes, cached replay is instant.
"""
from __future__ import annotations

import time

import numpy as np

from lab import Experiment, support, stats, order, v2_split
from windowed_arith import WindowedModExp, verify_modexp

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P1", "|S| at K=4 equals K=0 exactly, K=5 equals K=1 exactly")
exp.predict("P2", "the K=0 vs K=4 and K=1 vs K=5 support SETS coincide")
exp.predict("P3", "live exponent bits stay {0..alpha-1} at K=4,5 (C37)")
exp.must_fail("C1", "the beta>1 lookup control must grow at K=4,5")

W = 2
BETA1 = [(5, 4), (7, 6)]          # r = 2, alpha = 1, beta = 1
BETAG = (7, 3)                    # r = 6, beta = 3
KS = (0, 1, 2, 3, 4, 5)


def build(N, a, K):
    me = WindowedModExp(N=N, a=a, n_exp=W * (K + 1), w=W)
    assert verify_modexp(me), f"lookup N={N} a={a} K={K} incorrect"
    return me


def header():
    print(f"    {'N':>4} {'a':>3} {'K':>3} {'n_exp':>6} {'q':>4} "
          f"{'|support|':>12} {'vs K-2':>9} {'live exp bits':>15} {'secs':>8}")


# ---------------------------------------------------------------------------
exp.section("P1/P3  |support| and the live exponent bits, K = 0..5")
header()
runs = {}
for (N, a) in BETA1:
    alpha, beta = v2_split(order(a, N))
    for K in KS:
        me = build(N, a, K)
        qc = me.build()
        t0 = time.time()
        zs = support(qc, me.x[0], exact=True)
        live = [i for i, q_ in enumerate(me.exp)
                if bool(((zs >> q_) & 1).any())]
        prev = runs.get((N, a, K - 2))
        cmp = "-" if prev is None else ("SAME" if prev[0] == zs.size else "DIFF")
        print(f"    {N:4d} {a:3d} {K:3d} {me.n_exp:6d} {qc.n:4d} "
              f"{zs.size:12,d} {cmp:>9} {str(live):>15} "
              f"{time.time() - t0:8.1f}", flush=True)
        runs[(N, a, K)] = (int(zs.size), live, qc.n, alpha)
        # keep only the sets needed by P2 (K = 0,1,4,5); the rest are large
        if K in (0, 1, 4, 5):
            runs[(N, a, K)] = runs[(N, a, K)] + (zs,)
        else:
            del zs
    print()

p1 = []
for (N, a) in BETA1:
    p1.append(("K4==K0", N, a, runs[(N, a, 4)][0] == runs[(N, a, 0)][0]))
    p1.append(("K5==K1", N, a, runs[(N, a, 5)][0] == runs[(N, a, 1)][0]))
for tag, N, a, ok in p1:
    exp.log(f"N={N} a={a} {tag}: {'exact match' if ok else 'MISMATCH'}")
exp.check("P1", all(ok for _, _, _, ok in p1),
          "counts: " + "; ".join(
              f"N={N}: {[runs[(N, a, K)][0] for K in KS]}"
              for (N, a) in BETA1))

exp.section("P3  the tail stays dead however long it gets")
p3 = [(N, a, K, runs[(N, a, K)][1], list(range(runs[(N, a, K)][3])))
      for (N, a) in BETA1 for K in (4, 5)]
for N, a, K, got, want in p3:
    exp.log(f"N={N} a={a} K={K}: live exp bits {got}, expected {want}")
exp.check("P3", all(got == want for _, _, _, got, want in p3),
          f"{sum(got == want for _, _, _, got, want in p3)}/{len(p3)} widths "
          f"confine the support to the sub-2^alpha exponent bits")

exp.section("P2  the supports coincide as SETS across a two-window gap")
p2 = []
for (N, a) in BETA1:
    for lo, hi in ((0, 4), (1, 5)):
        zlo = runs[(N, a, lo)][4]
        zhi = runs[(N, a, hi)][4]
        low = (1 << runs[(N, a, lo)][2]) - 1        # the K=lo coordinates
        same = bool(np.array_equal(np.sort(zhi & low), np.sort(zlo)))
        p2.append(same)
        exp.log(f"N={N} a={a}: K={lo} vs K={hi} sets identical after dropping "
                f"the dead tail: {same}  (|S|={zlo.size:,})")
exp.check("P2", all(p2), f"{sum(p2)}/{len(p2)} pairs identical as sets")

# ---------------------------------------------------------------------------
exp.section("C1  must-fail control -- beta = 3 lookup has no dead tail")
header()
N, a = BETAG
ctrl = []
for K in KS:
    me = build(N, a, K)
    qc = me.build()
    t0 = time.time()
    st = stats(qc, me.x[0])
    prev = ctrl[-2][1] if len(ctrl) >= 2 else None
    cmp = "-" if prev is None else ("SAME" if prev == st["count"] else "DIFF")
    print(f"    {N:4d} {a:3d} {K:3d} {me.n_exp:6d} {qc.n:4d} "
          f"{st['count']:12,d} {cmp:>9} {'-':>15} {time.time() - t0:8.1f}",
          flush=True)
    ctrl.append((K, st["count"]))
grew = all(ctrl[i + 1][1] > ctrl[i][1] for i in range(len(ctrl) - 1))
exp.fail_check("C1", grew,
               f"beta=3 counts {[c for _, c in ctrl]} -- strictly increasing, "
               f"so the 2-periodicity above is a property of the dead tail")

exp.finish()
