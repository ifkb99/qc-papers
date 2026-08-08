"""Follow-up to experiment_weight.py. Two things to settle.

(a) The naive prediction "weight-truncation error = Fourier tail" FAILED.
    Why: PPS truncates INCREMENTALLY, at every gate. Truncating the final
    operator (terminal truncation) is a different operation from truncating
    throughout, because discarded terms would otherwise have kept branching.
    So the Fourier tail describes terminal truncation only. Verify that, and
    measure how far incremental drifts from it.

(b) The weight profile showed modexp is *strongly low-weight concentrated*:
    50% of spectral mass at weight <= 2, versus 2.6% for a random function
    (292x enrichment at weight 1). That is the opposite of what a "nonlinear,
    half-dense" spectrum would suggest, and it should be the same structure
    F13 found (4 of 15493 coefficients reproducing <O> exactly). Check whether
    the dominant coefficients ARE the low-weight ones, and which registers they
    live on -- this is TODO step 8, now directly motivated.
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
from perm_pps import propagate_perm
import walsh

me = ToffoliModExp(N=5, a=2, n_exp=1)
qc = me.build(); tq = me.x[0]
c = walsh.pullback_coefficients(qc, tq)
n = qc.n
z = np.arange(c.size, dtype=np.int64)
w = np.bitwise_count(z)
exact = float(c.sum())

print("=" * 78)
print("(a) TERMINAL vs INCREMENTAL weight truncation")
print("=" * 78)
print(f"  modexp N=5, {n} qubits, exact <O> = {exact:+.6f}\n")
print(f"  {'k':>3} {'terminal <O>':>13} {'term err':>10} {'incremental <O>':>16} "
      f"{'incr err':>10} {'tail mass':>10}")
for k in (1, 2, 3, 4, 6, 8, 10, 12, n):
    keep = w <= k
    term = float(c[keep].sum())
    r = propagate_perm(qc, 1 << tq, delta=0.0, max_weight=k)
    tail = float(np.sum(c[~keep] ** 2))
    print(f"  {k:3d} {term:+13.6f} {abs(term-exact):10.2e} {r.expectation:+16.6f} "
          f"{abs(r.expectation-exact):10.2e} {tail:10.6f}", flush=True)
print("""
  Terminal truncation is exactly described by the Walsh spectrum (it IS the
  spectrum, truncated). Incremental is path-dependent and can land either side
  of it -- discarded terms never branch, so the error is not the tail.""")

print()
print("=" * 78)
print("(b) ARE the dominant coefficients the low-weight ones?")
print("=" * 78)
order = np.argsort(-np.abs(c))
print(f"  top 12 coefficients by |c|, out of {int((np.abs(c)>1e-12).sum())} nonzero:")
print(f"  {'rank':>5} {'|c|':>9} {'weight':>7} {'z (bits set)':>28}")
reg = {}
for q_ in me.b: reg[q_] = "b"
for q_ in me.t: reg[q_] = "t"
for q_ in me.x: reg[q_] = "x"
reg[me.c0] = "c0"; reg[me.anc] = "anc"
for q_ in me.exp: reg[q_] = "e"
for i, idx in enumerate(order[:12]):
    bits = [j for j in range(n) if (int(idx) >> j) & 1]
    tags = ",".join(f"{reg.get(j,'?')}{j}" for j in bits)
    print(f"  {i+1:5d} {abs(c[idx]):9.6f} {int(idx).bit_count():7d} {tags:>28}")

big = np.abs(c) > 0.1
print(f"\n  coefficients with |c| > 0.1 : {int(big.sum())}")
print(f"    their weights            : {sorted(w[big].tolist())}")
print(f"    they sum to              : {float(c[big].sum()):+.8f}   (exact {exact:+.6f})")
print(f"    everything else sums to  : {float(c[~big].sum()):+.2e}")
print(f"    mass in the big ones     : {float(np.sum(c[big]**2)):.6f}")
print(f"    mass in the rest         : {float(np.sum(c[~big]**2)):.6f}")

print()
print("=" * 78)
print("(c) Does the high-weight sea cancel?")
print("=" * 78)
print(f"  {'weight k':>9} {'#terms':>8} {'mass':>10} {'signed sum':>13}")
for k in range(n + 1):
    sel = (w == k) & (np.abs(c) > 1e-12)
    if not sel.any():
        continue
    print(f"  {k:9d} {int(sel.sum()):8d} {float(np.sum(c[sel]**2)):10.6f} "
          f"{float(c[sel].sum()):+13.6f}")
print("""
  A weight level whose signed sum is ~0 contributes nothing to <O> however many
  terms it holds -- that is the cancelling sea. Levels with a large signed sum
  are the ones truncation must not break.""")
