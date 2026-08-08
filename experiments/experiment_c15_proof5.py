"""TODO step 3 -- the formal proof, and its sharpest testable consequence.

THEOREM (sketch, now with both previously-owed links closed).

Let ord_N(a) = r = 2^alpha (beta = 1). In the Toffoli modexp circuit the block
attached to exponent qubit e_i is u_a(e_i, a^(2^i)), and for i >= alpha we have
a^(2^i) = 1, so every such block is the SAME controlled operation C_{e_i}(V)
with

    V = u_a(., 1) = A^-1 . S . A

where A is the controlled multiply-accumulate and S the cswap layer.

(i)   S is a product of transpositions on DISJOINT qubit pairs (x_i, b_i), so
      S^2 = id. Hence V^2 = A^-1 S A A^-1 S A = A^-1 S^2 A = id.
      **V is an involution because it is a conjugate of one.**  [was numerical]

(ii)  V never modifies any exponent qubit (they appear only as controls), and
      all identity blocks apply powers of the SAME V, so the C_{e_i}(V) commute.

(iii) Therefore the identity tail composes to
          prod_{i >= alpha} C_{e_i}(V) = V^{sum e_i} = V^{p},  p = XOR_{i>=alpha} e_i
      using V^2 = id. The circuit depends on that whole tail only through the
      single parity bit p.                                      [was numerical]

(iv)  Let I = {i : i >= alpha} be the tail and write f(y, e_I) = F(y, p).
      Averaging the Walsh character over e_I:
          z_I = 0        ->  ( (-1)^F(y,0) + (-1)^F(y,1) ) / 2
          z_I = all-ones ->  ( (-1)^F(y,0) - (-1)^F(y,1) ) / 2
          otherwise      ->  0
      because any other z_I gives a character independent of the parity
      character. So the support is confined to z_I in {0, all-ones} -- TWO
      values, regardless of |I|. Hence

          |support| = #{z_y : c(z_y, 0) != 0} + #{z_y : c(z_y, all-ones) != 0}

      which does not depend on |I| at all.  QED.

This simultaneously gives: support size constant in n_exp; the added qubits
LIVE (they all appear in the all-ones mask); exactly one of each branch pair
surviving; and the onset at n_exp = alpha + 1 (the first identity block).

THE TEST BELOW is (iv)'s sharpest consequence, not previously checked:
    every support element z must have z restricted to the identity tail equal
    to either all-zeros or all-ones.
"""
from __future__ import annotations
import numpy as np
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


print("=" * 78)
print("Test of (iv): is the support confined to z_I in {0, all-ones}?")
print("=" * 78)
print(f"  {'N':>4} {'a':>3} {'r':>3} {'al':>3} {'n_exp':>6} {'|I|':>4} "
      f"{'|support|':>10} {'z_I=0':>9} {'z_I=all1':>9} {'other':>7} {'holds':>6}")
for N, a, ne in ((7, 6, 3), (7, 6, 4), (7, 6, 5), (5, 4, 4), (5, 2, 5),
                 (5, 2, 4), (21, 8, 3)):
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    qc = me.build()
    r = order(a, N); al = v2(r)
    if (r >> al) != 1:
        continue
    tail = [me.exp[i] for i in range(ne) if i >= al]
    if not tail:
        continue
    mask = 0
    for q_ in tail:
        mask |= 1 << q_
    c = walsh.pullback_coefficients(qc, me.x[0])
    zs = np.nonzero(np.abs(c) > 1e-12)[0]
    sub = zs & mask
    n0 = int(np.count_nonzero(sub == 0))
    n1 = int(np.count_nonzero(sub == mask))
    other = zs.size - n0 - n1
    print(f"  {N:4d} {a:3d} {r:3d} {al:3d} {ne:6d} {len(tail):4d} {zs.size:10d} "
          f"{n0:9d} {n1:9d} {other:7d} {str(other == 0):>6}", flush=True)
    del c

print()
print("=" * 78)
print("CONTROL: beta > 1 has no identity tail, so the confinement must FAIL")
print("=" * 78)
print(f"  {'N':>4} {'a':>3} {'r':>3} {'n_exp':>6} {'|support|':>10} "
      f"{'z_I=0':>9} {'z_I=all1':>9} {'other':>8} {'holds':>6}")
for N, a, ne in ((7, 3, 4), (21, 4, 3)):
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    qc = me.build()
    r = order(a, N); al = v2(r)
    # use the same nominal tail indices for comparison
    tail = [me.exp[i] for i in range(ne) if i >= al]
    mask = 0
    for q_ in tail:
        mask |= 1 << q_
    c = walsh.pullback_coefficients(qc, me.x[0])
    zs = np.nonzero(np.abs(c) > 1e-12)[0]
    sub = zs & mask
    n0 = int(np.count_nonzero(sub == 0))
    n1 = int(np.count_nonzero(sub == mask))
    other = zs.size - n0 - n1
    print(f"  {N:4d} {a:3d} {r:3d} {ne:6d} {zs.size:10d} {n0:9d} {n1:9d} "
          f"{other:8d} {str(other == 0):>6}", flush=True)
    del c

print("""
  If 'other' is 0 for every beta=1 row and large for the controls, then (iv) is
  confirmed and the support-size formula is proved: the count splits into the
  z_I=0 and z_I=all-ones halves, neither of which depends on how many qubits
  the tail spans.""")
