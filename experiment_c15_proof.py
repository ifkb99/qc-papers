"""TODO step 3 -- attempt a proof of the circuit-level C15 invariance.

Target: for beta = 1 and n_exp >= v2(r)+1, adding another exponent qubit leaves
the Walsh support size unchanged, even though the added qubits are LIVE (so it
is not a trivial embedding).

Proof strategy. Adding an exponent qubit appends the block U = u_a(ctrl, 1),
which is the identity on the valid subspace but not on the whole space. The
pullback becomes

    pi'^dag Z pi'  =  pi^dag (U^dag Z U) pi

For a permutation p, conjugation acts on a diagonal operator by permuting its
diagonal: D -> D o p. Permuting a diagonal does NOT generally preserve Walsh
sparsity -- EXCEPT when p is AFFINE over GF(2). If p(y) = My + c with M
invertible, then the Walsh coefficients are merely relabelled,

    c'_z = (-1)^{<z, ...>} c_{M^T z}      =>  |support| is preserved EXACTLY.

So: **if u_a(ctrl, 1) is affine over GF(2), the invariance is proved.** X and
CNOT are affine; Toffoli is not, so this is a real question about whether the
nonlinearities cancel in the composition.

H1  Is u_a(ctrl, 1) affine over GF(2)?
H2  If so, is the whole identity-block-suffix affine?
H3  Sanity: the non-identity blocks u_a(ctrl, a) with a != 1 should NOT be
    affine -- otherwise the argument would prove constancy for every r, which
    is false.
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
from circuits import Circuit
import walsh


def is_affine_perm(perm: np.ndarray, n: int):
    """Is y -> perm(y) affine over GF(2)?  perm(y) = perm(0) XOR (+) L(e_i).

    Returns (bool, n_violations). Builds L from the basis images, then checks
    the prediction on every input.
    """
    c0 = int(perm[0])
    L = [int(perm[1 << i]) ^ c0 for i in range(n)]
    idx = np.arange(perm.size, dtype=np.int64)
    pred = np.full(perm.size, c0, dtype=np.int64)
    for i in range(n):
        bit = (idx >> i) & 1
        pred ^= np.where(bit == 1, np.int64(L[i]), np.int64(0))
    bad = int(np.count_nonzero(pred != perm))
    return bad == 0, bad


def block_circuit(me: ToffoliModExp, ctrl: int, a: int) -> Circuit:
    qc = Circuit(me.n_qubits)
    qc.extend(me.u_a(ctrl, a))
    return qc


print("=" * 78)
print("H1  Is the identity block u_a(ctrl, 1) affine over GF(2)?")
print("=" * 78)
print(f"  {'N':>4} {'q':>4} {'block':>18} {'affine?':>9} {'violations':>12}")
for N, a in ((5, 2), (5, 4), (7, 6), (7, 3)):
    me = ToffoliModExp(N=N, a=a, n_exp=2)
    ctrl = me.exp[0]
    for aa, label in ((1, "u_a(ctrl, 1)"), (a % N, f"u_a(ctrl, {a%N})")):
        qc = block_circuit(me, ctrl, aa)
        perm = walsh.classical_permutation(qc)
        ok, bad = is_affine_perm(perm, qc.n)
        print(f"  {N:4d} {qc.n:4d} {label:>18} {str(ok):>9} {bad:12d}", flush=True)
    print()

print("=" * 78)
print("H2  Direct consequence check: does conjugating by the identity block")
print("    preserve Walsh sparsity for a diagonal observable?")
print("=" * 78)
print("  Compare support of Z_j vs support of U^dag Z_j U, for U = u_a(ctrl,1).\n")
print(f"  {'N':>4} {'q':>4} {'|supp Z_j|':>11} {'|supp U+ Z_j U|':>16} {'equal?':>8}")
for N, a in ((5, 2), (7, 6)):
    me = ToffoliModExp(N=N, a=a, n_exp=2)
    qc = block_circuit(me, me.exp[0], 1)
    tq = me.x[0]
    c = walsh.pullback_coefficients(qc, tq)
    s_after = int((np.abs(c) > 1e-12).sum())
    print(f"  {N:4d} {qc.n:4d} {1:11d} {s_after:16d} "
          f"{str(s_after == 1):>8}", flush=True)

print("""
  |supp Z_j| = 1 trivially (it IS a single Z-string). If conjugation by the
  identity block leaves it a single Z-string, the block acts affinely on the
  observable and sparsity is preserved term-for-term.""")

print()
print("=" * 78)
print("H3  Does affineness track the free/expensive branch?")
print("=" * 78)
print("  If u_a(ctrl, a) were affine for every a, the argument would prove")
print("  constancy for every r -- which is false. So non-identity blocks must")
print("  be non-affine.\n")
print(f"  {'N':>4} {'a arg':>6} {'affine?':>9} {'violations':>12}  note")
me = ToffoliModExp(N=7, a=3, n_exp=2)
for aa in (1, 2, 3, 4, 5, 6):
    qc = block_circuit(me, me.exp[0], aa)
    perm = walsh.classical_permutation(qc)
    ok, bad = is_affine_perm(perm, qc.n)
    note = "identity on valid subspace" if aa == 1 else ""
    print(f"  {7:4d} {aa:6d} {str(ok):>9} {bad:12d}  {note}", flush=True)
