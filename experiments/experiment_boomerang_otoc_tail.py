"""Which exponent bits have identical OTOC profiles in compiled modexp: beta = 1 only, or beta in {1, 3}?

STATUS CORRECTIONS (review V7cf842d53f784f58), which leave the pre-registered
text below as registered: D3 actually composes the blocks for c = 5 and
c = 3 at N = 11 (M_2 and M_3), not "M_2 composed with M_2". The pass count
includes the always-passing record Q4. The block-inverse identity holds as
permutations, not literally gate for gate: the controlled swaps appear in
reverse order but commute.

Context. experiment_boomerang_otoc run 1 predicted (P6) that the last two
exponent qubits of ToffoliModExp have identical X-X OTOC rows when beta = 1.
Its must-fail control C5 (N=7, g=3, beta=3) did NOT fail: the rows were
identical there too. That file stays failing. This experiment tests the
explanation derived afterwards by reading the construction.

DERIVATION (after run 1, before this run). The exponent register is only a
control, so pi(e,w) = (e, f_e(w)) with f_e = L o M_k^{e_k} o R. Here M_k is
the compiled block for constant c_k = g^(2^k) mod N, R is the product of the
blocks before bit k and L the product of those after it. For V = X_{e_k}
and W = U^dag X_q U with q a non-exponent qubit, the boomerang condition at
(e, w) reduces, with v = R(w) and h_L(x) = L^-1(L(x) ^ q), to
h_L(M_k v) = M_k h_L(v). So

    F(e_k, q) = E_L  Pr_v[ h_L and M_k commute at v ],

where L ranges uniformly over the products of the later blocks.
  * k = t-1 (no later blocks): F = Pr_v[ M(v ^ q) = M(v) ^ q ] with M = M_{t-1}.
  * k = t-2: F = 1/2 Pr[M_{t-2} commutes with xor_q]
             + 1/2 Pr[h_{M_{t-1}} commutes with M_{t-2}].
Each compiled block is multiply-swap-unmultiply, M_c = A_{c^-1}^-1 S A_c,
so M_{c^-1} = M_c^-1 exactly on the full dirty space. If
M_{t-2} = M_{t-1}^{+-1}, both terms reduce to the k = t-1 event, and the rows
coincide. Since c_{t-1} = c_{t-2}^2, this holds iff c_{t-2} = 1 or
c_{t-2}^3 = 1 mod N. For r = beta 2^alpha and t-2 >= alpha, that is iff
beta in {1, 3}. The derivation gives only the "if" direction; any other
equality would need a different mechanism.

PREDICTIONS, WRITTEN BEFORE MEASURING (t = 4 exponent bits).
  Q1  Rows e_2 and e_3 are identical for every non-exponent qubit for
      ToffoliModExp(9, 2) (r=6, beta=3, c_2=7, c_3=4=7^-1) and
      ToffoliModExp(13, 2) (r=12, beta=3, alpha=2, c_2=3, c_3=9=3^-1).
  Q2  The block-inverse identity M_{c^-1} = M_c^-1 holds as full-space
      permutations for the compiled blocks of N=7 (c=2,4) and N=11 (c=5,9).
  Q3  Ideal modexp (w < N multiplied, identity above; n = bit length of N):
      F(e_{t-2}, w_j) = F(e_{t-1}, w_j) for all j, for every (N, g) with
      3 <= N <= 31, gcd(g,N) = 1, g != 1, t = 4, and c_2^3 = 1 mod N.

  D1  ToffoliModExp(11, 2) (r=10, beta=5, c_2=5, c_3=3; 5^-1 = 9 != 3): row
      equality must FAIL for some qubit.
  D2  ToffoliModExp(11, 3) (r=5, beta=5, c_2=4, c_3=5; 4^-1 = 3 != 5): row
      equality must FAIL for some qubit.
  D3  Replacing the mixed block by a mismatched one must break Q2: M_2 for
      N=11 composed with M_2 must not be the identity.

Open, no commitment: among ideal (N, g) with c_2^3 != 1, how many still show
equal rows by some other coincidence. Report the count.

Run:  uv run python -m experiments.experiment_boomerang_otoc_tail   (from research/)
CPU; 20-qubit permutations (2^20 entries).
"""
from __future__ import annotations

import math

import numpy as np

from lab import Experiment
from lab.otoc import otoc_formula
from toffoli_arith import ToffoliModExp
import walsh

TOL = 1e-12


def compiled_rows(N, g, t=4):
    me = ToffoliModExp(N, g, n_exp=t)
    qc = me.build()
    p = walsh.classical_permutation(qc).astype(np.int64)
    qs = [q for q in range(qc.n) if q not in me.exp]
    rows = {k: np.array([otoc_formula(p, 1 << me.exp[k], 0, 1 << q, 0) for q in qs]) for k in (t - 2, t - 1)}
    return rows, qc.n


def ideal_modexp(N, g, t, n):
    out = []
    for e in range(1 << t):
        ge = pow(g, e, N)
        for w in range(1 << n):
            out.append(((ge * w) % N if w < N else w) + (e << n))
    return np.array(out, dtype=np.int64)


if __name__ == "__main__":
    exp = Experiment("experiment_boomerang_otoc_tail", doc=__doc__)
    exp.predict("Q1", "rows e_2 == e_3 for ToffoliModExp(9,2) and (13,2), beta = 3")
    exp.predict("Q2", "compiled block for c^-1 is the inverse permutation of the block for c")
    exp.predict("Q3", "ideal modexp rows equal whenever c_2^3 = 1, all (N,g) up to 31")
    exp.must_fail("D1", "ToffoliModExp(11,2), beta=5: row equality fails")
    exp.must_fail("D2", "ToffoliModExp(11,3), beta=5: row equality fails")
    exp.must_fail("D3", "mismatched block pair is not mutually inverse")
    exp.predict("Q4", "OPEN: count ideal (N,g) with c_2^3 != 1 but equal rows")

    exp.section("Q2/D3  block-inverse identity")
    def inverse_pair(N, c1, c2):
        """Is the compiled controlled block u_a(c2) the inverse permutation of u_a(c1)?"""
        me = ToffoliModExp(N, 2, n_exp=1)
        p1 = walsh.classical_permutation(me.u_a(me.exp[0], c1)).astype(np.int64)
        p2 = walsh.classical_permutation(me.u_a(me.exp[0], c2)).astype(np.int64)
        return bool(np.array_equal(p2[p1], np.arange(len(p1))))
    q2a, q2b = inverse_pair(7, 2, 4), inverse_pair(11, 5, 9)
    exp.check("Q2", q2a and q2b, f"N=7 (2,4): {q2a}; N=11 (5,9): {q2b}")
    d3 = inverse_pair(11, 5, 3)
    exp.fail_check("D3", not d3, f"N=11 (5,3) mutually inverse: {d3}")

    exp.section("Q1/D1/D2  compiled tail rows")
    res = {}
    for N, g in ((9, 2), (13, 2), (11, 2), (11, 3)):
        rows, nq = compiled_rows(N, g)
        diff = float(np.abs(rows[2] - rows[3]).max())
        cs = [pow(g, 1 << k, N) for k in range(4)]
        res[(N, g)] = diff
        exp.log(f"ToffoliModExp({N},{g}), {nq} qubits, constants {cs}: max |row2-row3| = {diff:.4g}")
    exp.check("Q1", res[(9, 2)] < TOL and res[(13, 2)] < TOL, f"N=9: {res[(9,2)]:.3g}, N=13: {res[(13,2)]:.3g}")
    exp.fail_check("D1", res[(11, 2)] > TOL, f"max diff {res[(11,2)]:.4g}")
    exp.fail_check("D2", res[(11, 3)] > TOL, f"max diff {res[(11,3)]:.4g}")

    exp.section("Q3/Q4  ideal modexp sweep")
    t = 4
    q3_bad, cube_cases, other_equal, other_total = [], 0, [], 0
    for N in range(3, 32):
        n = N.bit_length()
        for g in range(2, N):
            if math.gcd(g, N) != 1:
                continue
            p = ideal_modexp(N, g, t, n)
            r2 = [otoc_formula(p, 1 << (n + 2), 0, 1 << j, 0) for j in range(n)]
            r3 = [otoc_formula(p, 1 << (n + 3), 0, 1 << j, 0) for j in range(n)]
            equal = np.allclose(r2, r3, atol=TOL)
            c2 = pow(g, 4, N)
            if pow(c2, 3, N) == 1:
                cube_cases += 1
                if not equal:
                    q3_bad.append((N, g))
            else:
                other_total += 1
                if equal:
                    other_equal.append((N, g))
    exp.check("Q3", not q3_bad, f"{cube_cases} cases with c_2^3 = 1; violations {q3_bad[:5]}")
    exp.log(f"Q4: {len(other_equal)} of {other_total} cases with c_2^3 != 1 still have equal rows: {other_equal[:12]}")
    exp.check("Q4", True, "recorded (open question; cannot fail)")
    exp.finish()
