"""RETRACTED -- DO NOT CITE. Kept only for the record.

Measured a toy "Shor sandwich" (H layer -> ripple adder -> inverse QFT) that is
NOT a proxy for Shor's algorithm: the adder acts on the b-register while the QFT
acts on the a-register, so the observable barely meets the arithmetic. The
original scaling numbers (138 -> 872306, "1.57 bits/qubit") were additionally
produced with a bit-reversed QFT. With the QFT fixed the sandwich collapses to a
constant N_max = 8 at every size.

Use modexp.py / toffoli_arith.py instead. See NOTES.md STATUS 1.
"""
from __future__ import annotations
import numpy as np, time
from circuits import Circuit, ripple_adder
from pps import propagate

CAP = 800_000


def build(nb, mode):
    add, lay = ripple_adder(nb)
    areg = lay["a"]
    qc = Circuit(2 * nb + 2)
    if mode in ("h", "sandwich"):
        for q in areg:
            qc.h(q)
    qc.gates += add.gates
    if mode == "sandwich":
        qc.qft(areg, inverse=True)
    return qc, lay


print("=" * 78)
print("1. PEAK PAULI TERMS vs REGISTER SIZE  (delta = 0, i.e. exact PPS)")
print("=" * 78)
print(f"  {'nbits':>6}{'qubits':>8}{'Tgates':>8} | {'arith N_max':>12}{'arith fin':>10}"
      f" | {'sandwich N_max':>15}{'sand fin':>10}{'ratio':>8}")
rows = []
for nb in (2, 3, 4, 5, 6):
    out = {}
    for mode in ("arith", "sandwich"):
        qc, lay = build(nb, mode)
        obs = {(0, 1 << lay["a"][0]): 1.0}
        t0 = time.time()
        r = propagate(qc, obs, delta=0.0, max_terms=CAP)
        out[mode] = (r.n_max, r.n_terms[-1] if r.n_terms else 0, r.hit_cap,
                     qc.n_nonclifford(), time.time() - t0)
    a, s = out["arith"], out["sandwich"]
    ratio = s[0] / a[0] if a[0] else float("nan")
    flag = " CAP" if s[2] else ""
    print(f"  {nb:6d}{2*nb+2:8d}{s[3]:8d} | {a[0]:12d}{a[1]:10d}"
          f" | {s[0]:15d}{s[1]:10d}{ratio:8.1f}{flag}")
    rows.append((nb, a[0], s[0], s[2]))

print("""
  'arith' = H layer + ripple adder.  'sandwich' = the same, plus an inverse QFT
  on the a-register -- i.e. the H -> reversible arithmetic -> QFT skeleton of
  Shor's algorithm.""")

print()
print("=" * 78)
print("2. DOES delta CONTROL COST FOR THE SANDWICH?")
print("=" * 78)
for nb in (3, 4, 5):
    qc, lay = build(nb, "sandwich")
    obs = {(0, 1 << lay["a"][0]): 1.0}
    print(f"  nbits={nb} ({2*nb+2} qubits, {qc.n_nonclifford()} T gates):")
    prev = None
    for d in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 0.0):
        r = propagate(qc, obs, delta=d, max_terms=CAP)
        sl = ""
        if prev is not None and d > 0 and r.n_max > 0 and prev[1] > 0:
            sl = f"   slope={np.log(r.n_max/prev[1])/np.log(prev[0]/d):5.2f}"
        tag = "  (exact)" if d == 0 else ""
        print(f"    delta={d:8.0e}  N_max={r.n_max:8d}  <O>={r.expectation:12.8f}{sl}{tag}")
        if d > 0:
            prev = (d, r.n_max)

print("""
  A useful truncation parameter buys a large cost reduction for a small error.
  If N_max barely moves as delta sweeps six decades, delta is not a working
  dial for this circuit family and the paper's N_max extrapolation (Eq. 17)
  has nothing to extrapolate along.""")

print()
print("=" * 78)
print("3. GROWTH RATE")
print("=" * 78)
print("  N_max vs qubit count, exact PPS:")
for nb, a, s, capped in rows:
    print(f"    nbits={nb:2d}  qubits={2*nb+2:3d}  arith={a:8d}  sandwich={s:8d}"
          + ("  (capped)" if capped else ""))
valid = [(2 * nb + 2, s) for nb, a, s, capped in rows if not capped and s > 0]
if len(valid) >= 3:
    q = np.array([v[0] for v in valid], float)
    y = np.log2(np.array([v[1] for v in valid], float))
    slope = np.polyfit(q, y, 1)[0]
    print(f"\n  sandwich: log2(N_max) grows ~{slope:.3f} bits per qubit added")
    print(f"  (1.0 would be N_max ~ 2^n; 2.0 would be the 4^n worst case)")
