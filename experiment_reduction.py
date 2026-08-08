"""TODO step 10 -- does a different modular reduction lift the density ceiling?

Established (NOTES.md §L): the support carries a linear structure
w = b_msb XOR anc, from the CNOT pairing in cc_add_mod that sets and uncomputes
the comparison flag. It confines the support to a hyperplane, capping density at
exactly 1/2, and it is present in BOTH compilations because they share the
add / subtract-N / conditional-restore reduction.

Question: is the ancilla discipline really the cause? If a reduction that
entangles the flag with something OTHER than the sign bit removes the structure
and lifts the ceiling, the diagnosis is confirmed -- and it would be the FIRST
genuine cost increase from a compilation choice found in this project, since
everything else has turned out compilation-invariant.

DESIGN. The variant must not change what the circuit computes, or the comparison
is meaningless. The lever: conjugate the ancilla with an operation controlled on
SCRATCH qubits, which are 0 on the valid subspace (so a^e mod N is untouched)
but not on the full space (so the permutation, and hence the pullback, changes).

  V0  baseline                      -- expect defect 1, w = msb XOR anc
  V1  conjugate anc with CNOT from t[0]
  V2  conjugate anc with Toffoli from t[0], t[1]
  V3  conjugate anc with CNOT from a NON-msb b qubit (b[0])

PREDICTIONS.
  P1  V0 reproduces defect 1 with w = msb XOR anc.
  P2  At least one variant reaches full rank (defect 0), lifting the ceiling.
  P3  CORRECTNESS: every variant must still compute a^e mod N exactly. Checked
      first -- a variant that breaks correctness proves nothing.
  P4  If a variant has defect 1 but with a DIFFERENT w, that is also
      informative: the ceiling would be structural to the reduction pattern
      rather than to the specific qubit pair.
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
from circuits import Circuit
import walsh
import statevec as sv

TOL = 1e-12


class Variant(ToffoliModExp):
    """ToffoliModExp with the ancilla conjugated by a scratch-controlled op."""

    def __init__(self, *a, mode="v0", **k):
        self.mode = mode
        super().__init__(*a, **k)

    def _wrap(self, qc):
        if self.mode == "v0":
            return
        if self.mode == "v1":
            qc.cnot(self.t[0], self.anc)
        elif self.mode == "v2":
            qc.toffoli(self.t[0], self.t[1], self.anc)
        elif self.mode == "v3":
            qc.cnot(self.b[0], self.anc)
        else:
            raise ValueError(self.mode)

    def cc_add_mod(self, qc, c1, c2, c):
        self._wrap(qc)
        super().cc_add_mod(qc, c1, c2, c)
        self._wrap(qc)


def rank_defect(qc, target, n):
    c = walsh.pullback_coefficients(qc, target)
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


print("=" * 84)
print("P3  CORRECTNESS FIRST -- every variant must still compute a^e mod N")
print("=" * 84)
N, a, ne = 5, 2, 2
ok_modes = []
for mode in ("v0", "v1", "v2", "v3"):
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
            good = False
            break
    print(f"  {mode}: computes a^e mod N correctly = {good}")
    if good:
        ok_modes.append(mode)

print()
print("=" * 84)
print("P1/P2/P4  Rank defect and linear structure per variant")
print("=" * 84)
print(f"  {'mode':>5} {'q':>4} {'|supp|':>9} {'density':>9} {'rank':>5} "
      f"{'defect':>7} {'w (qubits)':>26}")
for mode in ok_modes:
    for N_, a_ in ((5, 2), (7, 6)):
        me = Variant(N=N_, a=a_, n_exp=2, mode=mode)
        qc = me.build()
        reg = {}
        for q_ in me.b: reg[q_] = "b"
        for q_ in me.t: reg[q_] = "t"
        for q_ in me.x: reg[q_] = "x"
        reg[me.c0] = "c0"; reg[me.anc] = "anc"
        for q_ in me.exp: reg[q_] = "e"
        S, rank, kern = rank_defect(qc, me.x[0], qc.n)
        wtxt = "none (FULL RANK)"
        if kern:
            bits = [i for i in range(qc.n) if (kern[0] >> i) & 1]
            wtxt = ",".join(f"{reg.get(i,'?')}{i}" for i in bits)
        print(f"  {mode:>5} {qc.n:4d} {S:9d} {S/(1<<qc.n):9.6f} {rank:5d} "
              f"{qc.n-rank:7d} {wtxt:>26}   (N={N_})", flush=True)

print("""
  P2 wants at least one variant at defect 0 with density climbing above 1/2.
  Defect 1 with a DIFFERENT w (P4) would mean the ceiling is structural to the
  reduction pattern rather than to the specific msb/anc pair.""")
