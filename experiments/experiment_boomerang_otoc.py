"""Are OTOCs of reversible arithmetic boomerang counts, and do they see the r = beta*2^alpha invariant?

STATUS (after run 1, 2026-09-14): exits nonzero ON PURPOSE. Must-fail control
C5 did NOT fail: beta = 3 (N=7, g=3) also has identical tail rows, so P6 does
not isolate beta = 1. The mechanism (compiled blocks satisfy M_{c^-1} = M_c^-1,
giving equality IF c_{t-2}^3 = 1 or c_{t-2} = 1; converse unproved) was derived
afterwards and tested in experiment_boomerang_otoc_tail.
CORRECTIONS (review V7cf842d53f784f58), which leave the pre-registered text
below as registered:
  * The adder derivation's "when the b-part holds, the z-part holds" misses
    j = m-1, where the b-part always holds and every failure is a carry-out
    flip. Formula (A) is still exact; the repaired proof is in C101.
  * "W ... an involution": sigma is an involution, and W^2 = (-1)^(b.e) I.
  * P2 (adder6) cannot detect orientation; the anchor is PRESENT in P3. See claims/C101.md and notes/BO.
Do not re-scope C5 to make this file pass. Run 3 onward registers the same
predictions and controls with literal calls (the lint cannot see loop-registered
controls); runs 1-2 used an equivalent loop.

Context. C8/C25 tie diagonal Pauli-path cost to linear cryptanalysis and C99
ties off-diagonal cost to the difference table. Operator scrambling has not
been studied in this project. The out-of-time-order correlator
F = 2^-n Tr(W(t)^dag V^dag W(t) V) is the standard scrambling diagnostic.

DERIVATION (before measuring). U|y> = |pi(y)>, W = U^dag X^b Z^e U,
V = X^a Z^d with Z applied first. W|y> = s(y)|sigma(y)> with
sigma(y) = pi^-1(pi(y)^b), an involution, s(y) = (-1)^(e.pi(y)).
Following |y> through V, W, V^dag, W^dag returns to |y> iff
sigma(y^a) ^ a = sigma(y), and the phase is (-1)^(d.(y^sigma(y))) s(y) s(y^a):

    F = 2^-n sum_{y in B} (-1)^(d.(y^sigma(y)) ^ e.(pi(y)^pi(y^a))),
    B = {y : pi^-1(pi(y)^b) ^ pi^-1(pi(y^a)^b) = a}.                     (O)

|B| is exactly beta_pi(a,b), the Boomerang Connectivity Table of Cid et al.
(Boura & Canteaut, ToSC 2018(3), Definition 3, p.295). So the X-X OTOC is
F = beta_pi(a,b)/2^n.

Consequences read from that paper (body):
  * APN permutations: BCT = DDT for a,b != 0 (p.292, citing CHP+18).
  * Inverse map over F_{2^n}, n even: max_{a,b!=0} beta = 6 if n = 0 mod 4 and
    4 if n = 2 mod 4 (Proposition 6, p.303).
  * PRESENT S-box: beta(1,5) = 16 (Table 5 and Example 5, p.302-303).

Adder (derived). Consider b += a + c0 (mod 2^m) with dirty carry-in c0 and
carry-out z ^= carry, and c0, a restored. Let t = (a+b+c0) mod 2^m, which
is uniform and independent of a. Take V = X on a_i and W = U^dag X_{b_j} U.
The b-part of the boomerang condition holds iff shifting t by +-2^i leaves
bit j of t unchanged. When it holds, the z-part holds as well, because the
overflow of b +- 2^j depends only on that bit. Hence

    F(a_i, b_j) = 1 (j < i),   0 (j = i),   1 - 2^(i-j) (i < j < m).   (A)

Invariant (derived). Ideal modexp on (e, w): pi(e,w) = (e, g^e w mod N) for
w < N, identity for w >= N. If g^(2^k) = 1 mod N, flipping e_k does not
change pi's work output, so B is everything and F(X_{e_k}, X_{w_j}) = 1 for
every j. With r = beta*2^alpha this happens iff beta = 1 and k >= alpha.

Compiled tail symmetry (derived from C23/C46). For beta = 1, the compiled
blocks for exponent bits k >= alpha are the same controlled involution V and
are applied consecutively. The swap of two tail exponent qubits therefore
commutes with U and fixes every non-exponent Pauli, so
F(X_{e_k}, X_q) is equal for all tail k at every non-exponent qubit q.

PREDICTIONS, WRITTEN BEFORE MEASURING.
  P1  (O) equals dense traces for every label (a,d,b,e), 4^n x 4^n, on n=4
      S-boxes (inversion over F16, PRESENT, a random permutation, an affine
      map). For n=5 (Gold x^3, inversion F32) it holds on all X-X labels and on
      500 fixed-seed signed labels.
  P2  The 6-qubit Cuccaro adder built from its Clifford+T unitary: (O) with
      walsh.classical_permutation equals dense traces for every pair of
      single-qubit X/Y/Z operators (324 labels).
  P3  Boura-Canteaut: for APN (Gold x^3 F32, inversion F32), 2^n F(a,b) =
      DDT_pi(a,b) for all a,b != 0. Inversion F16 has max 2^n F = 6 and
      inversion F64 has max 4, both over a,b != 0 (X-X). PRESENT has
      2^n F(1,5) = 16. Affine: F = 1 for all X-X labels.
  P4  (A) holds exactly for every i, j < m on the compiled Cuccaro adders
      m = 2, 3, 4 (full space including c0, z).
  P5  Ideal modexp: F(X_{e_k}, X_{w_j}) = 1 for all j at k = 2, 3 for N=5 g=2
      (r=4); t=4 exponent bits, n=3 work bits.
  P6  ToffoliModExp(5, 2, n_exp=4) (beta=1, alpha=2): F(X_{e_2}, X_q) =
      F(X_{e_3}, X_q) for every non-exponent qubit q.

  C1  Orientation: (O) with pi^-1 in place of pi must disagree with dense on
      some X-X label of PRESENT (its BCT is not symmetric).
  C2  "2^n F = DDT" must fail somewhere on inversion F16 (a,b != 0).
  C3  (O) with the sign dropped must disagree with dense on some signed label.
  C4  Ideal modexp: F < 1 for some j at k = 0 and k = 1 (N=5, g=2), and at
      every k (N=7, g=3, r=6, beta=3).
  C5  ToffoliModExp(7, 3, n_exp=4) (beta=3): tail equality between e_2 and
      e_3 must fail for some q.

Run:  uv run python -m experiments.experiment_boomerang_otoc   (from research/)
CPU; largest permutation 2^17.
"""
from __future__ import annotations

import numpy as np

from lab import Experiment
from lab.otoc import (otoc_formula, otoc_dense_perm, otoc_dense_unitary, bct, ddt)
from lab.differential import power_map
from circuits import ripple_adder
from toffoli_arith import ToffoliModExp
import walsh

PRESENT = np.array([0xC, 5, 6, 0xB, 9, 0, 0xA, 0xD, 3, 0xE, 0xF, 8, 4, 7, 1, 2], dtype=np.int64)
TOL = 1e-9


def affine16(rng):
    while True:
        A = rng.integers(0, 2, size=(4, 4))
        if round(abs(np.linalg.det(A))) % 2 == 1:
            break
    ys = np.arange(16)
    bits = (ys[:, None] >> np.arange(4)) & 1
    return ((((bits @ A.T) % 2) * (1 << np.arange(4))).sum(axis=1) ^ 0b1011).astype(np.int64)


def ideal_modexp(N, g, t, n):
    out = []
    for e in range(1 << t):
        ge = pow(g, e, N)
        for w in range(1 << n):
            out.append(((ge * w) % N if w < N else w) + (e << n))
    p = np.array(out, dtype=np.int64)
    assert len(set(p.tolist())) == len(p)
    return p


if __name__ == "__main__":
    exp = Experiment("experiment_boomerang_otoc", doc=__doc__)
    exp.predict("P1", "(O) == dense traces: all labels n=4 S-boxes; X-X + 500 signed at n=5")
    exp.predict("P2", "(O) == dense Clifford+T unitary OTOC, 324 single-qubit labels, adder6")
    exp.predict("P3", "APN 2^nF = DDT; inv F16 max 6, inv F64 max 4; PRESENT F(1,5)=1; affine F=1")
    exp.predict("P4", "adder butterfly (A) exact for m=2,3,4 compiled Cuccaro")
    exp.predict("P5", "ideal modexp N=5 g=2: F(e_k, w_j)=1 for k=2,3, all j")
    exp.predict("P6", "ToffoliModExp(5,2,4): F(e_2,q) == F(e_3,q) for all non-exponent q")
    exp.must_fail("C1", "pi^-1 orientation disagrees on PRESENT X-X")
    exp.must_fail("C2", "2^nF = DDT fails on inversion F16")
    exp.must_fail("C3", "dropped sign disagrees on some signed label")
    exp.must_fail("C4", "ideal modexp F<1 somewhere at k=0,1 (N=5) and every k (N=7)")
    exp.must_fail("C5", "ToffoliModExp(7,3,4) tail equality fails")

    rng = np.random.default_rng(20260914)
    n4 = {"inv16": power_map(4, 14), "present": PRESENT,
          "rand16": rng.permutation(16).astype(np.int64), "affine16": affine16(rng)}
    n5 = {"gold32": power_map(5, 3), "inv32": power_map(5, 30)}

    exp.section("P1/C3  identity (O) against dense traces")
    bad, nlab, sign_mut = 0, 0, False
    for name, p in n4.items():
        for a in range(16):
            for d in range(16):
                for b in range(16):
                    for e in range(16):
                        f = otoc_formula(p, a, d, b, e)
                        g = otoc_dense_perm(p, a, d, b, e)
                        bad += abs(f - g) > TOL
                        nlab += 1
                        if (d or e) and not sign_mut and abs(otoc_formula(p, a, 0, b, 0) - g) > TOL:
                            sign_mut = True
    lab5 = rng.integers(0, 32, size=(500, 4))
    for name, p in n5.items():
        for a in range(32):
            for b in range(32):
                bad += abs(otoc_formula(p, a, 0, b, 0) - otoc_dense_perm(p, a, 0, b, 0)) > TOL
                nlab += 1
        for (a, d, b, e) in lab5:
            bad += abs(otoc_formula(p, a, d, b, e) - otoc_dense_perm(p, a, d, b, e)) > TOL
            nlab += 1
    exp.check("P1", bad == 0, f"{bad} mismatches over {nlab} labels")
    exp.fail_check("C3", sign_mut, "dropping the sign changes some n=4 signed label" if sign_mut else "never mattered")

    exp.section("P2  gate-level Clifford+T unitary, 6-qubit adder")
    add6, _ = ripple_adder(2)
    U6 = add6.to_unitary()
    p6 = walsh.classical_permutation(add6).astype(np.int64)
    singles = [(1 << q, 0) for q in range(6)] + [(0, 1 << q) for q in range(6)] + [(1 << q, 1 << q) for q in range(6)]
    bad2 = sum(abs(otoc_formula(p6, a, d, b, e) - otoc_dense_unitary(U6, a, d, b, e)) > TOL
               for (a, d) in singles for (b, e) in singles)
    exp.check("P2", bad2 == 0, f"{bad2} mismatches of {len(singles)**2}")

    exp.section("P3/C1/C2  Boura-Canteaut consequences")
    apn_ok = all(bct(p, a, b) == ddt(p, a, b) for p in n5.values() for a in range(1, 32) for b in range(1, 32))
    otoc_eq_bct = all(abs(otoc_formula(p, a, 0, b, 0) * len(p) - bct(p, a, b)) < TOL
                      for p in list(n4.values()) + list(n5.values()) for a in range(len(p)) for b in range(len(p)))
    inv64 = power_map(6, 62)
    m16 = max(round(otoc_formula(n4["inv16"], a, 0, b, 0) * 16) for a in range(1, 16) for b in range(1, 16))
    m64 = max(bct(inv64, a, b) for a in range(1, 64) for b in range(1, 64))
    pres = otoc_formula(PRESENT, 1, 0, 5, 0)
    aff = {otoc_formula(n4["affine16"], a, 0, b, 0) for a in range(16) for b in range(16)}
    exp.check("P3", apn_ok and otoc_eq_bct and m16 == 6 and m64 == 4 and pres == 1.0 and aff == {1.0},
              f"APN BCT=DDT {apn_ok}; 2^nF=BCT {otoc_eq_bct}; inv16 max {m16}; inv64 max {m64}; "
              f"PRESENT F(1,5)={pres}; affine {sorted(aff)}")
    c1 = any(abs(otoc_formula(PRESENT, a, 0, b, 0, orientation="inverse") - otoc_dense_perm(PRESENT, a, 0, b, 0)) > TOL
             for a in range(16) for b in range(16))
    exp.fail_check("C1", c1, "pi^-1 orientation disagrees on PRESENT" if c1 else "orientation never mattered")
    c2 = any(bct(n4["inv16"], a, b) != ddt(n4["inv16"], a, b) for a in range(1, 16) for b in range(1, 16))
    exp.fail_check("C2", c2, "BCT != DDT somewhere on inversion F16" if c2 else "BCT == DDT everywhere")

    exp.section("P4  adder butterfly (A)")
    p4bad, shown = 0, []
    for m in (2, 3, 4):
        qc, reg = ripple_adder(m)
        p = walsh.classical_permutation(qc).astype(np.int64)
        for i in range(m):
            for j in range(m):
                want = 1.0 if j < i else (0.0 if j == i else 1 - 2.0 ** (i - j))
                got = otoc_formula(p, 1 << reg["a"][i], 0, 1 << reg["b"][j], 0)
                p4bad += abs(got - want) > TOL
                if m == 4 and i == 0:
                    shown.append(f"F(a0,b{j})={got:g}")
    exp.check("P4", p4bad == 0, f"{p4bad} mismatches; m=4 row i=0: " + ", ".join(shown))

    exp.section("P5/C4  ideal modexp and the invariant")
    t, n = 4, 3
    p5 = ideal_modexp(5, 2, t, n)
    F5 = {k: [otoc_formula(p5, 1 << (n + k), 0, 1 << j, 0) for j in range(n)] for k in range(t)}
    exp.log("N=5 g=2: " + "; ".join(f"k={k}: {np.round(v, 4).tolist()}" for k, v in F5.items()))
    exp.check("P5", all(v == 1.0 for k in (2, 3) for v in F5[k]), "tail bits k=2,3")
    p7 = ideal_modexp(7, 3, t, n)
    F7 = {k: [otoc_formula(p7, 1 << (n + k), 0, 1 << j, 0) for j in range(n)] for k in range(t)}
    exp.log("N=7 g=3: " + "; ".join(f"k={k}: {np.round(v, 4).tolist()}" for k, v in F7.items()))
    c4 = all(min(F5[k]) < 1 for k in (0, 1)) and all(min(F7[k]) < 1 for k in range(t))
    exp.fail_check("C4", c4, "scrambling present wherever g^(2^k) != 1" if c4 else "some non-identity bit unscrambled")

    exp.section("P6/C5  compiled tail symmetry")
    def tail_rows(N, g):
        me = ToffoliModExp(N, g, n_exp=4)
        qc = me.build()
        p = walsh.classical_permutation(qc).astype(np.int64)
        qs = [q for q in range(qc.n) if q not in me.exp]
        rows = {k: np.array([otoc_formula(p, 1 << me.exp[k], 0, 1 << q, 0) for q in qs]) for k in range(4)}
        return rows, qc.n
    r5, nq5 = tail_rows(5, 2)
    exp.log(f"ToffoliModExp(5,2,4), {nq5} qubits: k=2 row {np.round(r5[2], 4).tolist()}")
    exp.log(f"                                  k=3 row {np.round(r5[3], 4).tolist()}")
    exp.check("P6", np.allclose(r5[2], r5[3], atol=TOL), f"max |diff| {np.abs(r5[2]-r5[3]).max():.3g}")
    r7, nq7 = tail_rows(7, 3)
    exp.log(f"ToffoliModExp(7,3,4): k=2 row {np.round(r7[2], 4).tolist()}")
    exp.log(f"                      k=3 row {np.round(r7[3], 4).tolist()}")
    exp.fail_check("C5", not np.allclose(r7[2], r7[3], atol=TOL), f"max |diff| {np.abs(r7[2]-r7[3]).max():.3g}")
    exp.finish()
