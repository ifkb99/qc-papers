"""Step 10, second pass. First hypothesis refuted; better explanation to test.

PASS 1 RESULT: conjugating the ancilla with extra CNOT/Toffoli couplings (from
t0, from t0&t1, from b0) left the structure completely unchanged -- defect 1
with the SAME w = b_msb XOR anc in every variant. So the linear structure is NOT
an artifact of the particular CNOT pairing in cc_add_mod.

WHY IT IS ROBUST. Working it out:

  (a) Flipping the msb of an m-bit register is the same operation as adding
      2^(m-1), and mod-2^m addition is commutative, so
          (b XOR 2^(m-1)) + c  ==  (b + c) XOR 2^(m-1).
      The msb flip therefore COMMUTES with every _add_const in the reduction and
      propagates untouched through the whole arithmetic.

  (b) anc is coupled to msb only through XOR (cnot(msb,anc), and the
      x/cnot/x uncomputation). Flipping both together is exactly the invariance
      of XOR, so anc's value is restored.

  (c) The msb is b[m-1] = b[n], and the cswaps in u_a use zip(x, b[:n]) -- they
      EXCLUDE the msb. So the msb never reaches the x register, and the observed
      bit cannot see the flip.

That explains why pass 1's variants failed: adding more XOR couplings cannot
break an XOR symmetry. Breaking it requires coupling anc to msb NONLINEARLY.

  P1  v4: toffoli(msb, t0, anc) conjugation -- anc ^= (msb AND t0), nonlinear in
      msb. Predict the symmetry BREAKS: defect 0, density above 1/2.
  P2  v5: toffoli(msb, t0, t1) conjugation -- touches msb nonlinearly but not
      anc. Predict it also breaks (a) or (c) rather than (b).
  P3  Correctness must survive: the wraps are controlled on scratch qubits that
      are 0 on the valid subspace.
  P4  Direct check of explanation (a): verify that msb-flip commutes with the
      adder, i.e. that b XOR 2^(m-1) then +c equals +c then XOR 2^(m-1).
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
from circuits import Circuit
import walsh
import statevec as sv

TOL = 1e-12


class Variant(ToffoliModExp):
    def __init__(self, *a, mode="v0", **k):
        self.mode = mode
        super().__init__(*a, **k)

    def _wrap(self, qc):
        msb = self.b[self.m - 1]
        if self.mode == "v0":
            return
        if self.mode == "v4":
            qc.toffoli(msb, self.t[0], self.anc)      # anc ^= msb AND t0
        elif self.mode == "v5":
            qc.toffoli(msb, self.t[0], self.t[1])     # msb used nonlinearly
        elif self.mode == "v6":
            qc.cswap(self.t[0], msb, self.anc)        # swap msb<->anc if t0
        else:
            raise ValueError(self.mode)

    def cc_add_mod(self, qc, c1, c2, c):
        self._wrap(qc)
        super().cc_add_mod(qc, c1, c2, c)
        self._wrap(qc)


def analyse(qc, target, n):
    c = walsh.pullback_coefficients(qc, target)
    zs = np.nonzero(np.abs(c) > TOL)[0]
    basis = []
    for v in zs.tolist():
        cur = int(v)
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur); basis.sort(reverse=True)
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


print("=" * 86)
print("P4  Check explanation (a): does an msb flip commute with mod-2^m addition?")
print("=" * 86)
m = 4
M = 1 << m
bad = 0
for b in range(M):
    for c in range(M):
        if ((b ^ (1 << (m - 1))) + c) % M != (((b + c) % M) ^ (1 << (m - 1))):
            bad += 1
print(f"  m={m}: violations over all {M*M} (b,c) pairs = {bad}  -> "
      f"{'COMMUTES' if bad == 0 else 'does not commute'}")
print("  (so an msb flip propagates through every _add_const untouched)")

print()
print("=" * 86)
print("P3  Correctness of the nonlinear variants")
print("=" * 86)
N, a, ne = 5, 2, 2
ok = []
for mode in ("v0", "v4", "v5", "v6"):
    me = Variant(N=N, a=a, n_exp=ne, mode=mode)
    good = True
    for e in range(1 << ne):
        qc = Circuit(me.n_qubits)
        for i, eq in enumerate(me.exp):
            if (e >> i) & 1:
                qc.x(eq)
        qc.extend(me.build())
        j, amp = sv.peak(sv.run(qc))
        if sv.read_register(j, me.x) != pow(a, e, N) or abs(amp - 1) > 1e-8:
            good = False; break
    print(f"  {mode}: computes a^e mod N correctly = {good}")
    if good:
        ok.append(mode)

print()
print("=" * 86)
print("P1/P2  Does a NONLINEAR coupling break the structure?")
print("=" * 86)
print(f"  {'mode':>5} {'N':>3} {'q':>4} {'|supp|':>9} {'density':>9} {'defect':>7} "
      f"{'w':>22}")
for mode in ok:
    for N_, a_ in ((5, 2), (7, 6)):
        me = Variant(N=N_, a=a_, n_exp=2, mode=mode)
        qc = me.build()
        reg = {}
        for q_ in me.b: reg[q_] = "b"
        for q_ in me.t: reg[q_] = "t"
        for q_ in me.x: reg[q_] = "x"
        reg[me.c0] = "c0"; reg[me.anc] = "anc"
        for q_ in me.exp: reg[q_] = "e"
        S, rank, kern = analyse(qc, me.x[0], qc.n)
        wtxt = "none (FULL RANK)"
        if kern:
            bits = [i for i in range(qc.n) if (kern[0] >> i) & 1]
            wtxt = ",".join(f"{reg.get(i,'?')}{i}" for i in bits)
        print(f"  {mode:>5} {N_:3d} {qc.n:4d} {S:9d} {S/(1<<qc.n):9.6f} "
              f"{qc.n-rank:7d} {wtxt:>22}", flush=True)
