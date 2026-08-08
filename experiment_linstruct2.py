"""Follow-up: is the linear structure ALGORITHMIC or a CONSTRUCTION ARTIFACT?

Found: modexp pullbacks have GF(2) rank n-1 with linear structure
w = {b_msb, anc} -- the sign bit XOR the modular-reduction ancilla. That caps
density at exactly 1/2 and explains C7's convergence to 0.498 from below.

The pairing is visible in the source: cc_add_mod does
    cnot(msb, anc)              # set the flag
    ... ; x(msb); cnot(msb, anc); x(msb)   # uncompute it
so flipping msb and anc together is a symmetry of the reduction step.

That raises a question with real consequences for Paper B. If the structure is
an artifact of THIS ancilla discipline, then C7's density constant 1/2 is
construction-specific (the Theta(2^n) conclusion survives regardless, since
1/2 * 2^n is still Theta(2^n)) -- and that must be said.

  Q1  Does the Fourier (Beauregard) compilation, with a different ancilla
      layout, have the same structure?
  Q2  Does w depend on the observable, or is it a property of the circuit?
  Q3  Does the rank deficiency stay exactly 1 across sizes, or grow?
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
from modexp import ModExp
import walsh

TOL = 1e-12


def kernel_of_support(c, n):
    zs = np.nonzero(np.abs(c) > TOL)[0]
    basis = []
    for v in zs.tolist():
        cur = int(v)
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    piv = {b.bit_length() - 1: b for b in basis}
    kern = []
    for free in range(n):
        if free in piv:
            continue
        w = 1 << free
        for p in sorted(piv):
            if (piv[p] & w).bit_count() % 2:
                w ^= (1 << p)
        kern.append(w)
    return zs.size, len(basis), kern


print("=" * 82)
print("Q2/Q3  Toffoli compilation: does w depend on the observable?")
print("=" * 82)
for N, a, ne in ((7, 6, 2), (5, 2, 2)):
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    qc = me.build()
    reg = {}
    for q_ in me.b: reg[q_] = "b"
    for q_ in me.t: reg[q_] = "t"
    for q_ in me.x: reg[q_] = "x"
    reg[me.c0] = "c0"; reg[me.anc] = "anc"
    for q_ in me.exp: reg[q_] = "e"
    print(f"\n  N={N} a={a} n_exp={ne}  (q={qc.n}, b_msb=b{me.b[-1]}, anc={me.anc})")
    print(f"    {'observable':>12} {'|supp|':>8} {'rank':>5} {'defect':>7} {'w (qubits)':>22}")
    for lbl, tq in ([(f"Z_x{i}", me.x[i]) for i in range(min(3, len(me.x)))]
                    + [("Z_b0", me.b[0]), ("Z_t0", me.t[0])]):
        c = walsh.pullback_coefficients(qc, tq)
        S, rank, kern = kernel_of_support(c, qc.n)
        wtxt = "-"
        if kern:
            bits = [i for i in range(qc.n) if (kern[0] >> i) & 1]
            wtxt = ",".join(f"{reg.get(i,'?')}{i}" for i in bits)
        print(f"    {lbl:>12} {S:8d} {rank:5d} {qc.n-rank:7d} {wtxt:>22}", flush=True)
        del c

print()
print("=" * 82)
print("Q1  Fourier (Beauregard) compilation -- different ancilla layout")
print("=" * 82)
for N, a, ne in ((5, 2, 1), (5, 2, 2), (7, 3, 1)):
    me = ModExp(N=N, a=a, n_exp=ne)
    qc = me.build()
    reg = {}
    for q_ in me.b: reg[q_] = "b"
    for q_ in me.x: reg[q_] = "x"
    reg[me.anc] = "anc"
    for q_ in me.exp: reg[q_] = "e"
    c = walsh.pullback_coefficients(qc, me.x[0])
    S, rank, kern = kernel_of_support(c, qc.n)
    wtxt = "none (full rank)"
    if kern:
        bits = [i for i in range(qc.n) if (kern[0] >> i) & 1]
        wtxt = ",".join(f"{reg.get(i,'?')}{i}" for i in bits)
    print(f"  N={N} a={a} n_exp={ne}  q={qc.n} b_msb=b{me.b[-1]} anc={me.anc}: "
          f"|supp|={S} rank={rank} defect={qc.n-rank}")
    print(f"    density={S/(1<<qc.n):.6f}   w = {wtxt}", flush=True)
    del c

print("""
  If the Fourier compilation shows the SAME defect-1 structure with w spanning
  its own sign bit and ancilla, the effect follows from the modular-reduction
  discipline that both share, not from one implementation. If it shows full
  rank, the density-1/2 ceiling is specific to the Toffoli construction and
  C7's constant must be reported as construction-dependent.""")
