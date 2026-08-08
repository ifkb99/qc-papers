"""TODO step 9 -- how far does the C15 proof actually reach?

The proof was stated for the multiply-swap-unmultiply construction, where the
a=1 block is V = A^-1 S A and S is a product of disjoint transpositions, giving
V^2 = id. But re-reading the proof, steps (i)-(iv) never use V's internals
beyond V^2 = id:

  (i)   V^2 = id
  (ii)  the controlled blocks commute (same V, controls V does not touch)
  (iii) they compose to V^p with p the parity of the controls
  (iv)  character averaging confines the support to z_I in {0, 1_I}

"Identity on the valid subspace" was CONTEXT for why a^(2^i)=1 produces such a
block -- it is not a proof ingredient. So the conjecture is:

    CONJECTURE. The invariance holds for ANY construction whose repeated block
    is an involution, regardless of how that block is built.

That is a much larger and cleanly checkable class than "Vedral/Beauregard-style
multiply-swap-unmultiply".

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Replace the identity block with a SYNTHETIC INVOLUTION (a controlled swap
      of two ancilla qubits, order 2). Support must be CONSTANT in the number of
      such blocks.
  P2  Replace it with a SYNTHETIC NON-INVOLUTION (a controlled 3-cycle on three
      ancilla qubits, order 3). Support must NOT be constant -- otherwise V^2=id
      is not the operative condition and the proof is resting on something else.
  P3  Verify directly that block^2 = id for the first and block^3 = id (but
      block^2 != id) for the second.

P2 is the control that must fail. If BOTH stay constant, the conjecture is
wrong about the mechanism and the real condition is weaker than involutivity.
"""
from __future__ import annotations
import numpy as np
from circuits import Circuit
from toffoli_arith import ToffoliModExp
import walsh

TOL = 1e-12


def order_of_perm(perm):
    """Multiplicative order of a permutation given as an image array."""
    idn = np.arange(perm.size)
    cur = perm.copy()
    for k in range(1, 13):
        if np.array_equal(cur, idn):
            return k
        cur = perm[cur]
    return None


def build_with_blocks(N, a, n_base, k, kind):
    """modexp prefix on n_base exponent qubits, then k synthetic blocks each
    controlled on its own fresh qubit."""
    me = ToffoliModExp(N=N, a=a, n_exp=n_base)
    base = me.build()
    nq = me.n_qubits + k
    qc = Circuit(nq)
    qc.gates = list(base.gates)
    qc.logical = list(base.logical)
    ctrls = [me.n_qubits + i for i in range(k)]
    for c in ctrls:
        if kind == "involution":            # controlled swap: order 2
            qc.cswap(c, me.b[0], me.b[1])
        elif kind == "cycle3":              # controlled 3-cycle: order 3
            qc.cswap(c, me.b[0], me.b[1])
            qc.cswap(c, me.b[1], me.b[2])
        else:
            raise ValueError(kind)
    return me, qc, ctrls


def single_block_order(N, a, kind):
    me = ToffoliModExp(N=N, a=a, n_exp=1)
    nq = me.n_qubits + 1
    qc = Circuit(nq)
    c = me.n_qubits
    if kind == "involution":
        qc.cswap(c, me.b[0], me.b[1])
    else:
        qc.cswap(c, me.b[0], me.b[1])
        qc.cswap(c, me.b[1], me.b[2])
    return order_of_perm(walsh.classical_permutation(qc))


print("=" * 84)
print("P3  Confirm the synthetic blocks have the intended order")
print("=" * 84)
for kind in ("involution", "cycle3"):
    o = single_block_order(7, 6, kind)
    print(f"  {kind:12s}: permutation order = {o}   "
          f"{'(involution)' if o == 2 else '(NOT an involution)'}")

print()
print("=" * 84)
print("P1/P2  Is the support constant in the number of blocks?")
print("=" * 84)
for N, a in ((7, 6), (5, 4)):
    for kind in ("involution", "cycle3"):
        print(f"\n  N={N} a={a}  block = {kind}")
        print(f"    {'k blocks':>9} {'qubits':>7} {'|support|':>10} {'vs k=1':>9}")
        base = None
        for k in (1, 2, 3, 4):
            me, qc, _ = build_with_blocks(N, a, 1, k, kind)
            c = walsh.pullback_coefficients(qc, me.x[0])
            S = int((np.abs(c) > TOL).sum())
            if base is None:
                base = S
            tag = "SAME" if S == base else f"{S - base:+d}"
            print(f"    {k:9d} {qc.n:7d} {S:10d} {tag:>9}", flush=True)
            del c

print("""
  P1 expects SAME down the involution columns; P2 expects growth down the
  cycle3 columns. If cycle3 also stays constant, involutivity is not the
  operative condition and the proof rests on something weaker.""")
