"""C7: are the small-n results pre-asymptotic?

SUPERSEDED IN RANGE, NOT IN CONTENT (2026-08-08). This file's circuit series
stops at 24 qubits because that was the CPU limit. `experiment_c7_scale.py`
re-runs it on the exact integer path -- all six rows reproduce bit-for-bit --
and extends it to **30 qubits** in two series, slope 1.006 bits/qubit. Cite
that one for the scaling numbers; this one is still the record of how C7 was
answered.

This was the standing threat to the whole project: every instance so far is
<= 15 qubits, against the 127 qubits of the circuits the PPS resource framework
was calibrated on.

F9 makes the threat tractable without a faster simulator. Since the Pauli term
count for a permutation circuit *equals* the Walsh sparsity of the pulled-back
bit function, the exact PPS cost can be computed in O(2^n n) without running
PPS at all. That reaches 24 qubits where PPS stalls near 15 -- a 512x larger
Hilbert space.

Two probes:
  A. circuit-level  -- the real Toffoli modexp, ancillas and all
  B. function-level -- the idealised map e -> bit j of (a^e mod N), which
                       isolates the algorithm from ancilla-layout artifacts
                       (the confound that invalidated the old C3/C4).

The quantity that decides C7 is the *density* sparsity / 2^n. If density stays
near a constant as n grows, PPS cost is Theta(2^n) for this circuit family --
asymptotically no better than a state-vector simulation -- and the small-n
results are representative rather than pre-asymptotic.
"""
from __future__ import annotations
import numpy as np, time, math
from toffoli_arith import ToffoliModExp
from pps import propagate
import walsh

# ---------------------------------------------------------------- validation
print("=" * 78)
print("A0. VALIDATE the tool one size beyond where it was established")
print("=" * 78)
print("  (n_exp=1 keeps the PPS side tractable; the identity is already")
print("   established at n_exp=2 in experiment_c8.py, 6/6 instances)\n")
print(f"  {'N':>4} {'qubits':>7} {'walsh':>9} {'pps':>9} {'match':>6} {'t_walsh':>9} {'t_pps':>8}")
for N, a in ((5, 2), (7, 3), (15, 7)):
    me = ToffoliModExp(N=N, a=a, n_exp=1)
    qc = me.build(); q = me.x[0]
    t0 = time.time(); c = walsh.pullback_coefficients(qc, q); tw = time.time() - t0
    t0 = time.time()
    r = propagate(qc, {(0, 1 << q): 1.0}, delta=0.0, max_terms=4_000_000)
    tp = time.time() - t0
    ws = int((np.abs(c) > 1e-9).sum())
    print(f"  {N:4d} {me.n_qubits:7d} {ws:9d} {r.n_terms[-1]:9d} "
          f"{'YES' if ws == r.n_terms[-1] else 'NO':>6} {tw:8.2f}s {tp:7.1f}s", flush=True)

# ------------------------------------------------------------ A. circuit level
print()
print("=" * 78)
print("A. CIRCUIT-LEVEL SCALING  (real Toffoli modexp, exact PPS cost via Walsh)")
print("=" * 78)
print(f"  {'N':>4} {'n':>3} {'r':>4} {'qubits':>7} {'2^q':>11} {'sparsity':>10} "
      f"{'density':>8} {'time':>8}")
rows = []
for N, a in ((5, 2), (7, 3), (15, 7), (21, 2), (33, 5), (35, 3)):
    if math.gcd(a, N) != 1:
        continue
    me = ToffoliModExp(N=N, a=a, n_exp=2)
    qc = me.build(); q = me.x[0]
    t0 = time.time()
    c = walsh.pullback_coefficients(qc, q)
    el = time.time() - t0
    sp = int((np.abs(c) > 1e-12).sum())
    dim = 1 << me.n_qubits
    rows.append((me.n_qubits, sp, dim))
    print(f"  {N:4d} {me.n:3d} {me.order():4d} {me.n_qubits:7d} {dim:11d} "
          f"{sp:10d} {sp/dim:8.4f} {el:7.1f}s", flush=True)
    del c

print("\n  density = sparsity / 2^qubits.  0.5 would mean 'half of all possible")
print("  Z-strings are present' -- i.e. PPS carries Theta(2^n) terms.")

# ----------------------------------------------------------- B. function level
print()
print("=" * 78)
print("B. FUNCTION-LEVEL SCALING  g(e) = bit j of (a^e mod N), e in [0, 2^t)")
print("=" * 78)
print("  Isolates the algorithm from ancilla layout. a^e mod N has period r, so")
print("  the table is built once and tiled -- t can go well past any circuit.\n")
print(f"  {'N':>5} {'a':>3} {'r':>4} {'t':>3} {'2^t':>10} {'sparsity':>10} "
      f"{'density':>8} {'r pow2?':>8}")


def modexp_bit_spectrum(N, a, t, bit=0):
    r = 1; v = a % N
    while v != 1:
        v = (v * a) % N; r += 1
    table = np.empty(r, dtype=np.int64)
    x = 1
    for i in range(r):
        table[i] = x
        x = (x * a) % N
    reps = (1 << t) // r + 1
    vals = np.tile(table, reps)[: 1 << t]
    g = (vals >> bit) & 1
    chi = np.where(g == 1, -1.0, 1.0)
    c = walsh.wht(chi) / (1 << t)
    return r, int((np.abs(c) > 1e-12).sum())


for N, a in ((15, 7), (21, 2), (35, 3), (143, 5), (323, 5)):
    if math.gcd(a, N) != 1:
        continue
    for t in (12, 16, 20, 24):
        r, sp = modexp_bit_spectrum(N, a, t)
        ispow2 = (r & (r - 1)) == 0
        print(f"  {N:5d} {a:3d} {r:4d} {t:3d} {1<<t:10d} {sp:10d} "
              f"{sp/(1<<t):8.4f} {'yes' if ispow2 else 'no':>8}", flush=True)
    print()

# ------------------------------------------------------------------- verdict
print("=" * 78)
print("VERDICT ON C7")
print("=" * 78)
if len(rows) >= 3:
    q = np.array([r[0] for r in rows], float)
    s = np.array([r[1] for r in rows], float)
    slope = np.polyfit(q, np.log2(s), 1)[0]
    dens = s / np.array([r[2] for r in rows], float)
    print(f"  log2(sparsity) grows {slope:.3f} bits per qubit "
          f"(1.000 would be exactly Theta(2^n))")
    print(f"  density across {len(rows)} sizes: "
          + ", ".join(f"{d:.3f}" for d in dens))
