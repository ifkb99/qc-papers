"""TODO step 3, final piece -- the parity mechanism.

Attempt 3 found: h(y) = g(y XOR v) with v a SINGLE bit, always the *previous*
exponent qubit. That is the fingerprint of a parity dependence.

Proposed mechanism. For beta = 1 and i >= alpha, a^(2^i) mod N = 1, so every
such block applies the SAME fixed permutation V, each controlled on its own
exponent qubit. If V is an INVOLUTION (V^2 = id), the composite depends only on
how many of those controls are on, modulo 2 -- i.e. on the single parity bit

    p = XOR_{i >= alpha} e_i

Everything then follows:
  * the pulled-back function has ONE effective variable for the whole identity
    tail, however many qubits it spans        -> support size constant;
  * that variable involves every one of those qubits
                                              -> the added qubits stay LIVE;
  * toggling the new control has the same effect as toggling any other one
                                              -> h(y) = g(y XOR e_prev), the
                                                 character found in attempt 3;
  * with zero identity blocks there is no such variable, and the first one
    appears at i = alpha                      -> locks at n_exp = alpha + 1;
  * for beta > 1 no block is the identity, V differs per block, and none of it
    applies                                   -> support grows.

Tests:
  V1  Is V = u_a(ctrl, 1) an involution on the full space?
  V2  Does the pulled-back function depend on the identity-block controls ONLY
      through their parity? Equivalent test: flipping any TWO of them together
      leaves the function unchanged.
  V3  Control: for beta > 1, is the corresponding block an involution? (Should
      generally not be.)
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
from circuits import Circuit
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


def block_perm(me, ctrl, a):
    qc = Circuit(me.n_qubits)
    qc.extend(me.u_a(ctrl, a))
    return walsh.classical_permutation(qc)


print("=" * 78)
print("V1  Is the identity block V = u_a(ctrl, 1) an involution?")
print("=" * 78)
print(f"  {'N':>4} {'a arg':>6} {'q':>4} {'V^2 == id?':>11} {'#fixed by V^2':>15}")
for N, a in ((5, 2), (5, 4), (7, 6), (7, 3), (21, 8)):
    me = ToffoliModExp(N=N, a=a, n_exp=2)
    for aa in (1, a % N):
        p = block_perm(me, me.exp[0], aa)
        p2 = p[p]
        idn = np.arange(p.size)
        ok = bool(np.array_equal(p2, idn))
        print(f"  {N:4d} {aa:6d} {me.n_qubits:4d} {str(ok):>11} "
              f"{int(np.count_nonzero(p2 == idn)):15d}", flush=True)
    print()

print("=" * 78)
print("V2  Does the function depend on the identity-block controls ONLY")
print("    through their PARITY?")
print("=" * 78)
print("  Test: flip TWO identity-block controls together (parity preserved).")
print("  If the pulled-back bit function is unchanged, the dependence is parity.\n")
print(f"  {'N':>4} {'a':>3} {'r':>3} {'alpha':>6} {'n_exp':>6} {'id blocks':>10} "
      f"{'parity-only?':>13}")
for N, a, ne in ((7, 6, 4), (5, 4, 4), (5, 2, 5), (7, 3, 4)):
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    qc = me.build()
    r = order(a, N); al = v2(r)
    idblocks = [me.exp[i] for i in range(ne) if i >= al and pow(a, 1 << i, N) == 1]
    if len(idblocks) < 2:
        print(f"  {N:4d} {a:3d} {r:3d} {al:6d} {ne:6d} {len(idblocks):10d} "
              f"{'n/a (<2 blocks)':>13}")
        continue
    perm = walsh.classical_permutation(qc)
    g = (perm >> me.x[0]) & 1
    q0, q1 = idblocks[0], idblocks[1]
    flip = (1 << q0) | (1 << q1)
    idx = np.arange(perm.size, dtype=np.int64)
    same = bool(np.array_equal(g, g[idx ^ flip]))
    print(f"  {N:4d} {a:3d} {r:3d} {al:6d} {ne:6d} {len(idblocks):10d} "
          f"{str(same):>13}", flush=True)

print("""
  A True here means: flipping two identity-block controls simultaneously leaves
  the pulled-back bit function pointwise unchanged, so those controls enter only
  via their XOR. That is one effective variable regardless of how many qubits
  the identity tail spans -- which is exactly the invariance.""")
