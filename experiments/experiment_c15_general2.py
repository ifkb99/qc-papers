"""Step 9, second pass. The first pass was VACUOUS -- fixed here.

BUG in pass 1: the synthetic blocks acted only on b-register qubits while the
observable was Z_x0. A permutation that never touches bit x0 leaves
bit_x0(perm(y)) = bit_x0(y), so the pullback is unchanged for a trivial reason.
Both the involution and the 3-cycle came out constant, and the control failing
to fail is what exposed it.

Fix: the synthetic block must ACT ON the observed bit, and a vacuity check is
added -- adding one block must change the support relative to none.

  P1  involution touching x0 (controlled swap x0 <-> b0, order 2)
        -> support CONSTANT in the number of blocks
  P2  non-involution touching x0 (controlled 3-cycle x0 -> b0 -> b1, order 3)
        -> support NOT constant.  This is the control that must fail.
  P0  vacuity: support(k=1) != support(k=0) for both, or the test says nothing.
"""
from __future__ import annotations
import numpy as np
from circuits import Circuit
from toffoli_arith import ToffoliModExp
import walsh

TOL = 1e-12


def perm_order(perm, cap=13):
    idn = np.arange(perm.size)
    cur = perm.copy()
    for k in range(1, cap):
        if np.array_equal(cur, idn):
            return k
        cur = perm[cur]
    return None


def add_block(qc, ctrl, me, kind):
    """Synthetic block CONTROLLED on ctrl, acting on the observed bit x0."""
    if kind == "involution":
        qc.cswap(ctrl, me.x[0], me.b[0])                 # order 2
    elif kind == "cycle3":
        qc.cswap(ctrl, me.x[0], me.b[0])
        qc.cswap(ctrl, me.b[0], me.b[1])                 # order 3
    else:
        raise ValueError(kind)


def build(N, a, k, kind):
    me = ToffoliModExp(N=N, a=a, n_exp=1)
    base = me.build()
    qc = Circuit(me.n_qubits + k)
    qc.gates = list(base.gates)
    qc.logical = list(base.logical)
    for i in range(k):
        add_block(qc, me.n_qubits + i, me, kind)
    return me, qc


print("=" * 80)
print("P3  Order of the synthetic blocks (in isolation)")
print("=" * 80)
for kind in ("involution", "cycle3"):
    me = ToffoliModExp(N=7, a=6, n_exp=1)
    qc = Circuit(me.n_qubits + 1)
    add_block(qc, me.n_qubits, me, kind)
    o = perm_order(walsh.classical_permutation(qc))
    print(f"  {kind:12s}: order = {o}  "
          f"{'(involution)' if o == 2 else '(NOT an involution)'}")

print()
print("=" * 80)
print("P0/P1/P2  Support vs number of blocks (blocks now act on x0)")
print("=" * 80)
for N, a in ((7, 6), (5, 4)):
    for kind in ("involution", "cycle3"):
        print(f"\n  N={N} a={a}  block = {kind}")
        print(f"    {'k':>4} {'qubits':>7} {'|support|':>10} {'vs k=1':>9}")
        base = None
        vals = []
        for k in (0, 1, 2, 3, 4):
            me, qc = build(N, a, k, kind)
            c = walsh.pullback_coefficients(qc, me.x[0])
            S = int((np.abs(c) > TOL).sum())
            vals.append(S)
            if k == 1:
                base = S
            tag = "-" if k == 0 else ("SAME" if S == base else f"{S-base:+d}")
            print(f"    {k:4d} {qc.n:7d} {S:10d} {tag:>9}", flush=True)
            del c
        vacuous = vals[0] == vals[1]
        const = len(set(vals[1:])) == 1
        print(f"    -> vacuity check (k=0 vs k=1 differ): {not vacuous}")
        print(f"    -> constant for k>=1: {const}")

print("""
  Verdict rule. If involution=constant and cycle3=NOT constant, the conjecture
  holds: V^2 = id is the operative condition and the theorem covers any
  construction with an involutive repeated block. If both are constant even
  after the vacuity check passes, the real condition is weaker than
  involutivity and needs identifying.""")
