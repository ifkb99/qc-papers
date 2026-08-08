"""C12: is PPS-hardness of a reversible circuit the same thing as
cryptographic nonlinearity?

The bridge, exactly. With normalised coefficients c_z (so the unnormalised
Walsh value is W(z) = 2^n c_z), the standard cryptanalytic quantities are

    linearity     L(f)  = max_z |W(z)|            = 2^n max_z |c_z|
    nonlinearity  NL(f) = 2^(n-1) - L(f)/2        = 2^(n-1) (1 - max_z |c_z|)
    bent bound    NL_max = 2^(n-1) - 2^(n/2 - 1)

and Parseval gives sum_z c_z^2 = 1. So for a spectrum spread flat over S
nonzero coefficients, each has magnitude 1/sqrt(S), whence

    NL ~ 2^(n-1) (1 - 1/sqrt(S))          <-- ties PPS cost S to nonlinearity

Two extremes pin the correspondence exactly:
  * affine f  <=> S = 1  <=> NL = 0            <=> PPS holds ONE term
  * bent f    <=> S = 2^n, all |c| = 2^-(n/2)  <=> PPS holds EVERY Z-string,
                  NL maximal

So "reversible circuits that are expensive for Pauli propagation" and "Boolean
functions that resist linear approximation" are, for computational-basis
observables, the same class. This script measures where real arithmetic sits.
"""
from __future__ import annotations
import numpy as np, math
from circuits import ripple_adder
from toffoli_arith import ToffoliModExp
import walsh


def profile(c: np.ndarray, n: int, label: str):
    """Cryptanalytic profile from normalised Walsh coefficients."""
    absc = np.abs(c)
    S = int((absc > 1e-12).sum())
    maxc = float(absc.max())
    L = (1 << n) * maxc
    NL = 2 ** (n - 1) * (1 - maxc)
    NL_bent = 2 ** (n - 1) - 2 ** (n / 2 - 1)
    parseval = float((c ** 2).sum())
    flat_pred = 1.0 / math.sqrt(S) if S else float("nan")
    return dict(label=label, n=n, S=S, dens=S / (1 << n), maxc=maxc, L=L,
                NL=NL, ratio=NL / NL_bent if NL_bent > 0 else float("nan"),
                parseval=parseval, flat_pred=flat_pred)


def show(rows):
    print(f"  {'function':>26} {'n':>3} {'sparsity':>9} {'dens':>6} "
          f"{'max|c|':>8} {'1/sqrt(S)':>9} {'NL/NLbent':>10} {'Parseval':>9}")
    for r in rows:
        print(f"  {r['label']:>26} {r['n']:3d} {r['S']:9d} {r['dens']:6.3f} "
              f"{r['maxc']:8.5f} {r['flat_pred']:9.5f} {r['ratio']:10.5f} "
              f"{r['parseval']:9.6f}", flush=True)


print("=" * 92)
print("C12.1  Where does real arithmetic sit on the affine <-> bent axis?")
print("=" * 92)
rows = []

# affine end: ripple adder low bit
add, lay = ripple_adder(5)
c = walsh.pullback_coefficients(add, lay["b"][0])
rows.append(profile(c, add.n, "adder b0 (affine)"))
for i in (1, 2, 4):
    c = walsh.pullback_coefficients(add, lay["b"][i])
    rows.append(profile(c, add.n, f"adder b{i}"))

# nonlinear end: modular exponentiation
for N, a in ((5, 2), (7, 3), (15, 7), (21, 2)):
    if math.gcd(a, N) != 1:
        continue
    me = ToffoliModExp(N=N, a=a, n_exp=2)
    qc = me.build()
    c = walsh.pullback_coefficients(qc, me.x[0])
    rows.append(profile(c, qc.n, f"modexp N={N} a={a}"))

# baseline: uniformly random Boolean function
rng = np.random.default_rng(0)
for n in (14, 18):
    g = rng.integers(0, 2, size=1 << n)
    c = walsh.wht(np.where(g == 1, -1.0, 1.0)) / (1 << n)
    rows.append(profile(c, n, f"random f, n={n}"))

show(rows)

print("""
  Reading the columns:
    dens      = sparsity / 2^n.  1.000 means every Z-string is present.
    max|c|    vs 1/sqrt(S): equal => the spectrum is FLAT (Parseval-saturating).
    NL/NLbent = 1.000 would be a bent function, the maximally nonlinear case
                and simultaneously the PPS worst case.
    Parseval  = sum c_z^2, must be 1.000000; a check on the transform.""")

print()
print("=" * 92)
print("C12.2  Does the Parseval prediction NL ~ 2^(n-1)(1 - 1/sqrt(S)) hold?")
print("=" * 92)
print(f"  {'function':>26} {'NL (measured)':>14} {'NL (predicted)':>15} {'rel.err':>9}")
for r in rows:
    pred = 2 ** (r['n'] - 1) * (1 - r['flat_pred'])
    rel = abs(r['NL'] - pred) / max(abs(r['NL']), 1e-12)
    print(f"  {r['label']:>26} {r['NL']:14.1f} {pred:15.1f} {rel:9.4f}", flush=True)
print("""
  The prediction assumes a perfectly flat spectrum, so it is a bound rather
  than an identity: it is exact for affine and bent, and over-estimates NL
  wherever the spectrum has a heavy tail.""")

print()
print("=" * 92)
print("C12.3  The correspondence, stated as a table")
print("=" * 92)
print("""
    Boolean-function property   Walsh spectrum        PPS behaviour
    -------------------------   ------------------    --------------------------
    affine (NL = 0)             one spike             1 term, exact at any delta
    low nonlinearity            few large spikes      few terms, cheap
    high nonlinearity           spread out            many terms, expensive
    bent (NL maximal)           2^n flat coefficients every Z-string, worst case

  If this holds, decades of results bounding the nonlinearity of specific
  function families transfer directly into PPS resource estimation -- and,
  read the other way, a reversible circuit that is cheap for PPS is exactly one
  whose output bits admit a good linear approximation.""")
