"""C15: does the 2-adic structure of the order r actually set PPS cost for a
REAL circuit, or only for the idealised bit function?

The tension to resolve. At function level the dichotomy is absolute:
  r | 2^k         -> Walsh sparsity <= 2^k, CONSTANT in t
  r has odd factor -> density exactly 1.000
But at circuit level every instance measured so far sits at density ~0.5
regardless of r (N=5 r=4 -> 0.473; N=7 r=6 -> 0.474). Those two cannot both be
the whole story.

Controlled design: fix N (so circuit size, ancilla layout and gate structure are
identical) and vary a, which changes r and nothing else. N=7 gives
r in {3, 6, 2}: pure odd, mixed, pure power of two. Then sweep n_exp, because
with a 2-bit exponent register the periodicity in e has almost no room to show.

Hypothesis H_A: r-dependence is real but masked by ancilla/garbage structure,
                and emerges as the exponent register grows relative to it.
Hypothesis H_B: r-dependence is a property of the idealised function only, and
                circuit-level PPS cost is ~2^q/2 regardless -- in which case
                C15 is FALSE as stated for real circuits.
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


def live_variables(c, n, tol=1e-12):
    """How many input bits does the function actually depend on?

    Bit i is 'dead' if every Walsh coefficient with z_i set vanishes; a dead
    variable halves the achievable density.
    """
    nz = np.nonzero(np.abs(c) > tol)[0]
    live = 0
    for i in range(n):
        if np.any((nz >> i) & 1):
            live += 1
    return live


def ideal_density(N, a, t, bit=0):
    r = order(a, N)
    table = np.empty(r, dtype=np.int64); x = 1
    for i in range(r):
        table[i] = x; x = (x * a) % N
    vals = np.tile(table, (1 << t) // r + 1)[: 1 << t]
    chi = np.where(((vals >> bit) & 1) == 1, -1.0, 1.0)
    c = walsh.wht(chi) / (1 << t)
    return int((np.abs(c) > 1e-12).sum()) / (1 << t)


print("=" * 90)
print("C15.1  CONTROLLED: fix N=7 (identical circuit), vary a to change r")
print("=" * 90)
print(f"  {'a':>3} {'r':>3} {'v2(r)':>6} {'odd':>4} {'n_exp':>6} {'qubits':>7} "
      f"{'2^q':>9} {'sparsity':>10} {'density':>8} {'live':>5} {'ideal dens':>11}")
N = 7
for a in (2, 3, 6):                      # r = 3, 6, 2
    r = order(a, N)
    for ne in (2, 4, 6, 8):
        me = ToffoliModExp(N=N, a=a, n_exp=ne)
        qc = me.build()
        t0 = time.time()
        c = walsh.pullback_coefficients(qc, me.x[0])
        sp = int((np.abs(c) > 1e-12).sum())
        dim = 1 << me.n_qubits
        lv = live_variables(c, me.n_qubits)
        idl = ideal_density(N, a, max(ne, 4))
        print(f"  {a:3d} {r:3d} {v2(r):6d} {r>>v2(r):4d} {ne:6d} {me.n_qubits:7d} "
              f"{dim:9d} {sp:10d} {sp/dim:8.5f} {lv:5d} {idl:11.5f}"
              f"   ({time.time()-t0:.0f}s)", flush=True)
        del c
    print()

print("""  live = number of input bits the pulled-back function actually depends on.
  If live < qubits, the missing variables force density <= 2^-(qubits-live).
  ideal dens = density of the idealised g(e) = bit0(a^e mod N) at t = n_exp.""")

print()
print("=" * 90)
print("C15.2  CONFIRMATION at a second modulus: N=21, r in {6, 3, 2}")
print("=" * 90)
print(f"  {'a':>3} {'r':>3} {'v2(r)':>6} {'n_exp':>6} {'qubits':>7} {'sparsity':>10} "
      f"{'density':>8} {'live':>5}")
for a in (2, 4, 8):                      # r = 6, 3, 2
    r = order(a, 21)
    for ne in (2, 4):
        me = ToffoliModExp(N=21, a=a, n_exp=ne)
        qc = me.build()
        c = walsh.pullback_coefficients(qc, me.x[0])
        sp = int((np.abs(c) > 1e-12).sum()); dim = 1 << me.n_qubits
        print(f"  {a:3d} {r:3d} {v2(r):6d} {ne:6d} {me.n_qubits:7d} {sp:10d} "
              f"{sp/dim:8.5f} {live_variables(c, me.n_qubits):5d}", flush=True)
        del c
    print()

print("=" * 90)
print("VERDICT")
print("=" * 90)
print("""  H_A (r-dependence emerges at circuit level) is supported if density
  separates by v2(r) and the gap widens with n_exp.
  H_B (idealised-function artifact) is supported if density stays ~0.5 for all
  r -- in which case C15 must be restricted to the idealised function and the
  claim that it governs PPS cost for Shor circuits is withdrawn.""")
