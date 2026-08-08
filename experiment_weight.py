"""TODO step 2 -- weight truncation as Fourier tail mass.

PPS has two standard truncation knobs: coefficient threshold (delta) and Pauli
WEIGHT (drop terms acting on more than k qubits). All work so far covers only
delta. This closes the gap.

The observation. For a permutation circuit with a computational-basis
observable, the surviving Paulis are Z-strings Z^z, and the Pauli weight of Z^z
is popcount(z) -- which is exactly the Fourier DEGREE of the Walsh coefficient
c_z. Therefore weight-k truncation is precisely low-degree Fourier truncation of
the pulled-back Boolean function, and its error is the Fourier tail:

    operator L2 error^2 = sum_{|z| > k} c_z^2         (Parseval)
    <O> error           = |sum_{|z| > k} c_z|          (<0|Z^z|0> = 1 for all z)

Both are computable directly from the Walsh spectrum, with no simulation. If
that holds, the exact cost model extends from delta-truncation to
weight-truncation.

The practical question this decides: is weight truncation VIABLE for reversible
arithmetic? It is viable only if spectral mass sits at low degree. A uniformly
random Boolean function has mass distributed binomially, C(n,k)/2^n, peaked at
k = n/2 -- i.e. hopeless. Where does modular exponentiation sit?
"""
from __future__ import annotations
import numpy as np
from math import comb
from circuits import ripple_adder
from toffoli_arith import ToffoliModExp
from perm_pps import propagate_perm
import walsh


def mass_by_weight(c: np.ndarray, n: int):
    """Spectral mass sum(c_z^2) and signed sum(c_z), bucketed by popcount(z)."""
    z = np.arange(c.size, dtype=np.int64)
    w = np.bitwise_count(z)
    mass = np.zeros(n + 1)
    ssum = np.zeros(n + 1)
    for k in range(n + 1):
        sel = w == k
        mass[k] = float(np.sum(c[sel] ** 2))
        ssum[k] = float(np.sum(c[sel]))
    return mass, ssum


print("=" * 80)
print("1. VERIFY: is weight-truncation error exactly the Fourier tail?")
print("=" * 80)
print("  Predicted <O> error at cutoff k = |sum of c_z over |z| > k|,")
print("  computed from the Walsh spectrum alone. Measured = weight-truncated PPS.\n")
me = ToffoliModExp(N=5, a=2, n_exp=1)
qc = me.build(); tq = me.x[0]
c = walsh.pullback_coefficients(qc, tq)
mass, ssum = mass_by_weight(c, qc.n)
exact = float(c.sum())
print(f"  modexp N=5, {qc.n} qubits, exact <O> = {exact:+.6f}")
print(f"  {'k':>4} {'terms kept':>11} {'measured <O>':>13} {'predicted <O>':>14} "
      f"{'agree':>7} {'tail mass':>11}")
for k in (2, 4, 6, 8, 10, 12, qc.n):
    r = propagate_perm(qc, 1 << tq, delta=0.0, max_weight=k)
    pred = float(np.sum(ssum[: k + 1]))
    tail = float(np.sum(mass[k + 1:]))
    ok = abs(r.expectation - pred) < 1e-9
    print(f"  {k:4d} {len(r.final_terms):11d} {r.expectation:+13.6f} {pred:+14.6f} "
          f"{'YES' if ok else 'NO':>7} {tail:11.6f}", flush=True)

print()
print("=" * 80)
print("2. MASS-BY-WEIGHT PROFILES -- is low-degree truncation viable?")
print("=" * 80)


def profile(label, c, n, note=""):
    mass, _ = mass_by_weight(c, n)
    cum = np.cumsum(mass)
    peak = int(np.argmax(mass))
    # smallest k capturing 99% of the mass
    k99 = int(np.searchsorted(cum, 0.99)) if cum[-1] >= 0.99 else n
    frac_low = cum[min(3, n)]
    print(f"  {label:>26} n={n:3d}  peak weight={peak:3d}  k(99% mass)={k99:3d}  "
          f"mass at weight<=3: {frac_low:.6f}   {note}")
    return mass


print("  A uniformly random function has mass ~ C(n,k)/2^n, peaking at k=n/2.")
print("  Concentration at LOW weight is what makes weight truncation work.\n")

add, lay = ripple_adder(4)
profile("adder b0 (affine)", walsh.pullback_coefficients(add, lay["b"][0]), add.n)
profile("adder b2", walsh.pullback_coefficients(add, lay["b"][2]), add.n)
for N, a, ne in ((5, 2, 1), (7, 3, 1), (15, 7, 1)):
    m = ToffoliModExp(N=N, a=a, n_exp=ne); q = m.build()
    profile(f"modexp N={N}", walsh.pullback_coefficients(q, m.x[0]), q.n)

rng = np.random.default_rng(0)
for n in (14,):
    g = rng.integers(0, 2, size=1 << n)
    cr = walsh.wht(np.where(g == 1, -1.0, 1.0)) / (1 << n)
    profile(f"random f", cr, n, "(baseline)")

print()
print("=" * 80)
print("3. FULL PROFILE, modexp vs binomial baseline")
print("=" * 80)
me = ToffoliModExp(N=5, a=2, n_exp=1); qc = me.build()
c = walsh.pullback_coefficients(qc, me.x[0])
mass, _ = mass_by_weight(c, qc.n)
binom = np.array([comb(qc.n, k) for k in range(qc.n + 1)], dtype=float)
binom /= binom.sum()
print(f"  {'weight k':>9} {'modexp mass':>13} {'binomial':>11} {'ratio':>8}")
for k in range(qc.n + 1):
    if mass[k] < 1e-12 and binom[k] < 1e-6:
        continue
    rt = mass[k] / binom[k] if binom[k] > 0 else float('nan')
    print(f"  {k:9d} {mass[k]:13.6f} {binom[k]:11.6f} {rt:8.3f}")
print("""
  ratio ~1 everywhere would mean modexp is indistinguishable from a random
  function under weight truncation -- i.e. weight truncation is useless for it.""")
