"""TODO step 7 -- does the C15 invariance survive TRUNCATION?

Every C15 figure so far is delta = 0. Nobody runs PPS at delta = 0, so the
practical claim ("period-finding precision is free when r is a power of two")
is currently unsupported in the regime where it would be used.

PREDICTIONS, WRITTEN BEFORE MEASURING (derive-then-test).

Derivation. C24 established that for beta = 1 the support splits into the
z_I = 0 and z_I = 1_I halves with exactly one of each pair nonzero, and M2
established that the surviving magnitudes reproduce the smaller circuit's
spectrum BIT-FOR-BIT (max err 0.00e+00). So the *multiset of coefficient
magnitudes* is identical across n_exp for n_exp >= alpha + 1.

  P1  TERMINAL delta-truncation (threshold the finished spectrum) must preserve
      the constancy EXACTLY, at EVERY delta -- a magnitude threshold applied to
      identical magnitude multisets keeps identical counts. This should hold
      with zero exceptions; if it fails, C24 or M2 is wrong.

  P2  INCREMENTAL delta-truncation (what real PPS does) may BREAK it. Larger
      n_exp means more gates, hence more truncation events, hence more
      opportunity for the trajectory to diverge -- and W1 already showed
      incremental and terminal truncation are different operations with
      non-monotonic error.

  P3  CONTROL, must fail: for beta > 1 there should be no constancy at any
      delta, terminal or incremental.

  P4  <O> under aggressive truncation should also be constant across n_exp for
      beta = 1, since the dominant coefficients are magnitude-preserved.

If P1 holds and P2 fails, the honest statement for Paper B is that the
invariance is a property of the exact operator, and that a practitioner running
truncated PPS does not automatically inherit it.
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
from perm_pps import propagate_perm
import walsh

DELTAS = (0.0, 1e-6, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1)


def order(a, N):
    r, v = 1, a % N
    while v != 1:
        v = (v * a) % N; r += 1
    return r


print("=" * 84)
print("P1  TERMINAL truncation: threshold the exact spectrum")
print("=" * 84)
for N, a, label in ((7, 6, "N=7 a=6 (r=2, beta=1)"), (5, 2, "N=5 a=2 (r=4, beta=1)"),
                    (7, 3, "N=7 a=3 (r=6, beta=3) CONTROL")):
    print(f"\n  {label}")
    hdr = "  ".join(f"{d:>9.0e}" for d in DELTAS)
    print(f"    {'n_exp':>6} {hdr}")
    rows = []
    for ne in (3, 4, 5):
        me = ToffoliModExp(N=N, a=a, n_exp=ne)
        qc = me.build()
        c = walsh.pullback_coefficients(qc, me.x[0])
        ac = np.abs(c)
        counts = [int((ac > max(d, 1e-12)).sum()) for d in DELTAS]
        rows.append(counts)
        print(f"    {ne:6d} " + "  ".join(f"{v:9d}" for v in counts), flush=True)
        del c
    const = all(r == rows[0] for r in rows)
    print(f"    -> identical across n_exp at every delta: {const}")

print()
print("=" * 84)
print("P2/P4  INCREMENTAL truncation: what real PPS actually does")
print("=" * 84)
for N, a, label in ((7, 6, "N=7 a=6 (r=2, beta=1)"), (5, 2, "N=5 a=2 (r=4, beta=1)"),
                    (7, 3, "N=7 a=3 (r=6, beta=3) CONTROL")):
    print(f"\n  {label}")
    print(f"    {'n_exp':>6} {'delta':>9} {'N_max':>9} {'N_final':>9} {'<O>':>11}")
    per_delta = {}
    for ne in (3, 4, 5):
        me = ToffoliModExp(N=N, a=a, n_exp=ne)
        qc = me.build()
        for d in (0.0, 1e-4, 1e-2, 1e-1):
            r = propagate_perm(qc, 1 << me.x[0], delta=d, max_terms=8_000_000)
            per_delta.setdefault(d, []).append((r.n_max, len(r.final_terms),
                                                round(r.expectation, 9)))
            print(f"    {ne:6d} {d:9.0e} {r.n_max:9d} {len(r.final_terms):9d} "
                  f"{r.expectation:+11.6f}", flush=True)
    print("    constancy across n_exp, per delta:")
    for d, vals in per_delta.items():
        nm = "SAME" if all(v[0] == vals[0][0] for v in vals) else "differs"
        nf = "SAME" if all(v[1] == vals[0][1] for v in vals) else "differs"
        eo = "SAME" if all(v[2] == vals[0][2] for v in vals) else "differs"
        print(f"      delta={d:8.0e}  N_max {nm:8s} N_final {nf:8s} <O> {eo}")
