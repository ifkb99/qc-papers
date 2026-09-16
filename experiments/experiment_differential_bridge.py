"""Does differential cryptanalysis count off-diagonal Pauli paths the way linear cryptanalysis counts diagonal ones?

Context. C8/C12/C25 identify the PPS term count of a Z-type observable pulled
back through a basis permutation with Walsh sparsity, and import published
nonlinearity as a cost lower bound. F10 and C13 record that the identity
stops at X/Y-type observables: the pullback leaves the diagonal and the
propagator blew past its cap. This experiment tests a derived DIFFERENTIAL
counterpart that covers exactly those observables.

DERIVATION (before measuring). Let U|y> = |pi(y)>, a Pauli label (a,b) mean
X^a Z^b up to phase, and M = U^dag X^a Z^b U. Then

    M|y> = s(y) |sigma(y)>,  s(y) = (-1)^(b.pi(y)),  sigma(y) = pi^-1(pi(y)+a)

so the coefficient of X^c Z^d is

    alpha(c,d) = 2^-n  sum_{y in D_c} s(y) (-1)^(d.y),
    D_c = { y : y + sigma(y) = c }.

Substituting u = pi(y), |D_c| = DDT_{pi^-1}(a,c). Hence the X-support of the
pullback is exactly row a of the difference distribution table of pi^-1, and

    T(a,b) = sum_c | supp  W[ 1_{D_c} s ] |.                          (I)

sigma is an involution, so for a != 0 each D_c is a union of k = |D_c|/2
pairs {y, y+c}, and pi(y+c) = pi(y)+a. On one pair the transform is
(-1)^(d.y) s(y) (1 + (-1)^(d.c + b.a)), supported on the coset d.c = b.a.
Summing k pairs gives a coset of size 2^(n-1) times a signed sum of k terms:
  * k odd  -> never zero -> exactly 2^(n-1) terms;                    (K1)
  * k = 2  -> two independent affine constraints -> exactly 2^(n-2).   (K2)
Donoho-Stark on GF(2)^n gives |supp W| >= 2^n/|D_c|; with Cauchy-Schwarz
over sum_c |D_c| = 2^n and R_a = #{c : D_c nonempty} >= 2^n/delta,

    4^n/delta^2 <= sum_c 2^n/DDT(a,c) <= T(a,b) <= R_a 2^(n-1) <= 4^(n-1).   (B)

T(a,b) = 4^(n-1) for every a != 0 iff every nonzero DDT entry is 2, i.e.
iff pi is APN. For a differentially 4-uniform permutation, K1/K2 give

    T(a,b) = 2^(n-2) (2^n - 3 N4(a)),  N4(a) = #{c : DDT(a,c) = 4},    (F)

independent of b.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Identity (I): support set and |coefficients| from the DDT/Walsh formula
      equal an independent dense Pauli-transfer computation for EVERY Pauli
      (all 4^n labels) on every S-box fixture with n <= 6, and on the 6-qubit
      Cuccaro adder built from its gate-level Clifford+T unitary.
  P2  Identity (I) equals the existing gate-level PPS propagator's final
      support and magnitudes on the 8-qubit Cuccaro adder for the enumerated
      families X_q, Y_q (all q), X_q X_r (all q<r), X_q Z_r (all q != r).
  P3  a = 0 reduces to C8: T(0,b) equals Walsh sparsity of b.pi.
  P4  APN extremality: Gold x^3 and inversion over GF(2^3), GF(2^5) give
      T(a,b) = 16 and 256 respectively for every a != 0 and every b.
  P5  Formula (F): inversion over GF(2^4), GF(2^6) (delta = 4, N4 = 1 per row)
      gives T(a,b) = 52 and 976 respectively for every a != 0, every b.
  P6  Bounds (B) hold for every a != 0 and every b on every fixture,
      including one fixed-seed random permutation each at n = 4, 5, 6.
  P7  Class counts K1/K2: every D_c with k odd contributes exactly 2^(n-1),
      every D_c with k = 2 exactly 2^(n-2), on every fixture.
  P8  Affine permutation (random invertible GF(2) matrix plus constant,
      n = 5): T(a,b) = 1 for every label.

  C1  Wrong picture: using sigma'(y) = pi(pi^-1(y)+a) (Schroedinger rather
      than Heisenberg conjugation) must disagree with the dense reference on
      some label of the non-involutive Gold GF(2^5) and random n=5 fixtures.
  C2  Naive linear extrapolation "T(a,b) = Walsh sparsity of b.pi" must be
      wrong for some a != 0 label on the n=5 random fixture.
  C3  Dropping the sign s(y) must change the support set for some label with
      b.a odd (n=5 random fixture).
  C4  The unsquared mutant bound T >= 4^n/delta must be violated by the APN
      GF(2^5) fixture (predicted 512 > 256).

No commitment is made for rows with k >= 4 pairs (b-dependent cancellation
is possible there); they are covered only by (I) and (B).

Run:  uv run python -m experiments.experiment_differential_bridge   (from research/)
CPU only; n <= 8; dense references at n <= 6.
"""
from __future__ import annotations

import numpy as np

from lab import Experiment
from circuits import ripple_adder
from pps import propagate
import walsh

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P1", "formula (I) == dense Pauli transfer, all labels, n<=6 S-boxes and 6-qubit adder unitary")
exp.predict("P2", "formula (I) == gate-level PPS on 8-qubit adder, enumerated X/Y/XX/XZ families")
exp.predict("P3", "a=0: T(0,b) == Walsh sparsity of b.pi")
exp.predict("P4", "APN Gold/inverse GF(8), GF(32): T == 16 / 256 for all a!=0, all b")
exp.predict("P5", "inverse GF(16), GF(64): T == 52 / 976 for all a!=0, all b")
exp.predict("P6", "bounds (B) hold on every fixture, every a!=0, every b")
exp.predict("P7", "k odd classes contribute 2^(n-1), k=2 classes 2^(n-2), every fixture")
exp.predict("P8", "affine n=5: T == 1 for every label")
exp.must_fail("C1", "Schroedinger-direction sigma' disagrees with dense on a non-involution")
exp.must_fail("C2", "naive T = Walsh sparsity of b.pi is wrong for some a!=0")
exp.must_fail("C3", "dropping s(y) changes the support for some b.a odd label")
exp.must_fail("C4", "unsquared bound 4^n/delta violated by APN GF(32)")

SEED = 20260914
POLY = {3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011}


# ---------------------------------------------------------------- fixtures
def gf_mul(x: int, y: int, n: int) -> int:
    r = 0
    while y:
        if y & 1:
            r ^= x
        y >>= 1
        x <<= 1
        if (x >> n) & 1:
            x ^= POLY[n]
    return r


def gf_pow(x: int, e: int, n: int) -> int:
    r = 1
    while e:
        if e & 1:
            r = gf_mul(r, x, n)
        x = gf_mul(x, x, n)
        e >>= 1
    return r


def power_map(n: int, e: int) -> np.ndarray:
    p = np.array([gf_pow(x, e, n) if x else 0 for x in range(1 << n)], dtype=np.int64)
    assert len(set(p.tolist())) == 1 << n, f"x^{e} not a permutation of GF(2^{n})"
    return p


def random_perm(n: int, rng) -> np.ndarray:
    return rng.permutation(1 << n).astype(np.int64)


def affine_perm(n: int, rng) -> np.ndarray:
    while True:
        A = rng.integers(0, 2, size=(n, n))
        # rank over GF(2)
        M, r = A.copy(), 0
        for col in range(n):
            piv = next((i for i in range(r, n) if M[i, col]), None)
            if piv is None:
                continue
            M[[r, piv]] = M[[piv, r]]
            for i in range(n):
                if i != r and M[i, col]:
                    M[i] ^= M[r]
            r += 1
        if r == n:
            break
    const = int(rng.integers(0, 1 << n))
    ys = np.arange(1 << n)
    bits = (ys[:, None] >> np.arange(n)) & 1
    img = (bits @ A.T) % 2
    p = (img * (1 << np.arange(n))).sum(axis=1) ^ const
    assert len(set(p.tolist())) == 1 << n
    return p.astype(np.int64)


# ---------------------------------------------------------------- formula (I)
def fwht(v: np.ndarray) -> np.ndarray:
    """Unnormalised Walsh-Hadamard transform, written here (not walsh.py)."""
    v = v.astype(np.float64).copy()
    h = 1
    while h < len(v):
        v = v.reshape(-1, 2 * h)
        a, b = v[:, :h].copy(), v[:, h:].copy()
        v[:, :h], v[:, h:] = a + b, a - b
        v = v.reshape(-1)
        h *= 2
    return v


def parity(v: np.ndarray) -> np.ndarray:
    return np.bitwise_count(v) & 1


def formula_terms(perm, a, b, *, direction="heis", sign=True):
    """{(c,d): |alpha|} from identity (I), plus per-class detail."""
    n = int(len(perm)).bit_length() - 1
    ys = np.arange(1 << n, dtype=np.int64)
    inv = np.argsort(perm)
    if direction == "heis":
        sigma = inv[perm ^ a]
    else:                                  # C1 mutant
        sigma = perm[inv ^ a]
    s = np.where(parity(perm & b) == 1, -1.0, 1.0) if sign else np.ones(1 << n)
    cvec = ys ^ sigma
    out, classes = {}, []
    for c in np.unique(cvec):
        g = np.zeros(1 << n)
        mask = cvec == c
        g[mask] = s[mask]
        W = fwht(g)
        nz = np.nonzero(np.abs(W) > 0.5)[0]
        classes.append((int(c), int(mask.sum()), len(nz)))
        for d in nz:
            out[(int(c), int(d))] = abs(W[d]) / (1 << n)
    return out, classes


def ddt_row(perm, a):
    inv = np.argsort(perm)
    n = int(len(perm)).bit_length() - 1
    u = np.arange(1 << n)
    return np.bincount(inv[u ^ a] ^ inv[u], minlength=1 << n)


def delta_of(perm):
    n = int(len(perm)).bit_length() - 1
    return max(int(ddt_row(perm, a).max()) for a in range(1, 1 << n))


# ---------------------------------------------------------------- dense reference
_P1 = {(0, 0): np.eye(2), (1, 0): np.array([[0., 1.], [1., 0.]]),
       (0, 1): np.diag([1., -1.]), (1, 1): np.array([[0., 1.], [1., 0.]]) @ np.diag([1., -1.])}


def pauli_stack(n):
    """Real matrices X^x Z^z for all labels, index = x + (z << n), by explicit kron."""
    K, D = 4 ** n, 1 << n
    B = np.empty((K, D, D))
    for x in range(D):
        for z in range(D):
            m = np.eye(1)
            for q in reversed(range(n)):       # qubit 0 = least significant index bit
                m = np.kron(m, _P1[((x >> q) & 1, (z >> q) & 1)])
            B[x + (z << n)] = m
    # self-test of the ordering convention
    ys = np.arange(D)
    assert np.allclose(np.diag(B[1 << n]),
                       np.where(ys & 1, -1.0, 1.0))          # Z_0
    assert B[1 + 0][1, 0] == 1.0 and B[1][0, 1] == 1.0        # X_0 maps 0<->1
    return B


def dense_ptm_perm(perm, B):
    """PTM[i, j] = Tr(B_i^T U^T B_j U) / D with U|y> = |perm[y]>, via matrix products."""
    D = len(perm)
    U = np.zeros((D, D))
    U[perm, np.arange(D)] = 1.0
    K = B.shape[0]
    Q = np.einsum("ab,jbc,cd->jad", U.T, B, U, optimize=True)
    return (B.reshape(K, -1) @ Q.reshape(K, -1).T) / D


def dense_ptm_unitary(Uc, B):
    D = Uc.shape[0]
    K = B.shape[0]
    Q = np.einsum("ab,jbc,cd->jad", Uc.conj().T, B.astype(complex), Uc, optimize=True)
    return (B.reshape(K, -1).astype(complex) @ Q.reshape(K, -1).T) / D


def dense_terms(ptm_col, n):
    D = 1 << n
    nz = np.nonzero(np.abs(ptm_col) > 1e-9)[0]
    return {(int(i % D), int(i >> n)): float(abs(ptm_col[i])) for i in nz}


def same(t1, t2):
    if t1.keys() != t2.keys():
        return False
    return all(abs(t1[k] - t2[k]) < 1e-9 for k in t1)


# ---------------------------------------------------------------- build fixtures
rng = np.random.default_rng(SEED)
fixtures = {
    "gold3_gf8": power_map(3, 3),
    "inv_gf8": power_map(3, 6),
    "inv_gf16": power_map(4, 14),
    "gold3_gf32": power_map(5, 3),
    "inv_gf32": power_map(5, 30),
    "inv_gf64": power_map(6, 62),
    "rand4": random_perm(4, rng),
    "rand5": random_perm(5, rng),
    "rand6": random_perm(6, rng),
    "affine5": affine_perm(5, rng),
}
add6, _ = ripple_adder(2)
add8, _ = ripple_adder(3)
fixtures["adder6"] = walsh.classical_permutation(add6).astype(np.int64)
fixtures["adder8"] = walsh.classical_permutation(add8).astype(np.int64)
for name, p in fixtures.items():
    exp.log(f"{name:>11}: n={int(len(p)).bit_length()-1}  delta={delta_of(p)}")

# ---------------------------------------------------------------- P1 dense identity
exp.section("P1  formula (I) versus dense Pauli transfer, every label")
stacks = {}
p1_ok, p1_detail = True, []
T = {}          # (fixture, a, b) -> term count, reused below
for name, perm in fixtures.items():
    n = int(len(perm)).bit_length() - 1
    if n > 6:
        continue
    B = stacks.setdefault(n, pauli_stack(n))
    if name == "adder6":
        ptm = dense_ptm_unitary(add6.to_unitary(), B)
    else:
        ptm = dense_ptm_perm(perm, B)
    bad = 0
    for a in range(1 << n):
        for b in range(1 << n):
            f, _ = formula_terms(perm, a, b)
            T[(name, a, b)] = len(f)
            if not same(f, dense_terms(ptm[:, a + (b << n)], n)):
                bad += 1
    p1_detail.append(f"{name}:{bad}/{4**n}")
    p1_ok &= bad == 0
exp.check("P1", p1_ok, "mismatching labels " + ", ".join(p1_detail))

# ---------------------------------------------------------------- P2 gate-level PPS
exp.section("P2  formula (I) versus gate-level PPS on the 8-qubit adder")
perm8, n8 = fixtures["adder8"], 8
labels = [(1 << q, 0) for q in range(n8)] + [(1 << q, 1 << q) for q in range(n8)]
labels += [((1 << q) | (1 << r), 0) for q in range(n8) for r in range(q + 1, n8)]
labels += [(1 << q, 1 << r) for q in range(n8) for r in range(n8) if q != r]
bad, sizes = 0, []
for (a, b) in labels:
    f, _ = formula_terms(perm8, a, b)
    res = propagate(add8, {(a, b): 1.0}, delta=0.0, max_terms=4_000_000)
    assert not res.hit_cap
    g = {k: abs(v) for k, v in res.final_terms.items()}
    if not same(f, g):
        bad += 1
    sizes.append(len(f))
    T[("adder8", a, b)] = len(f)
exp.check("P2", bad == 0,
          f"{bad}/{len(labels)} labels mismatch; term counts {min(sizes)}..{max(sizes)}")

# ---------------------------------------------------------------- P3 a=0 is C8
exp.section("P3  a = 0 reduces to Walsh sparsity")
ok = True
for name, perm in fixtures.items():
    n = int(len(perm)).bit_length() - 1
    if n > 6:
        continue
    for b in range(1 << n):
        ws = int((np.abs(fwht(np.where(parity(perm & b) == 1, -1.0, 1.0))) > 0.5).sum())
        ok &= T[(name, 0, b)] == ws
exp.check("P3", ok, "all n<=6 fixtures, all b")

# ---------------------------------------------------------------- P4/P5 exact values
def all_offdiag(name, value):
    n = int(len(fixtures[name])).bit_length() - 1
    vals = {T[(name, a, b)] for a in range(1, 1 << n) for b in range(1 << n)}
    return vals == {value}, f"{name}: {sorted(vals)}"

exp.section("P4  APN extremality")
r = [all_offdiag(k, v) for k, v in
     [("gold3_gf8", 16), ("inv_gf8", 16), ("gold3_gf32", 256), ("inv_gf32", 256)]]
exp.check("P4", all(x for x, _ in r), "; ".join(d for _, d in r))

exp.section("P5  4-uniform closed form (F)")
r = [all_offdiag(k, v) for k, v in [("inv_gf16", 52), ("inv_gf64", 976)]]
exp.check("P5", all(x for x, _ in r), "; ".join(d for _, d in r))

# ---------------------------------------------------------------- P6/P7 bounds and classes
exp.section("P6/P7  bounds (B) and class counts K1/K2")
b_ok, k_ok, notes = True, True, []
tight = []
for name, perm in fixtures.items():
    n = int(len(perm)).bit_length() - 1
    dl = delta_of(perm)
    avals = range(1, 1 << n)
    bvals = range(1 << n) if n <= 6 else [0, 1, 0b10110101, (1 << n) - 1]
    if n > 6:
        avals = sorted({a for (a, _) in labels})
    for a in avals:
        row = ddt_row(perm, a)
        nzrow = row[row > 0]
        mid = int(np.sum((1 << n) / nzrow))
        up = len(nzrow) * (1 << (n - 1))
        for b in bvals:
            f, classes = formula_terms(perm, a, b)
            t = len(f)
            if not (4 ** n / dl ** 2 <= mid <= t <= up <= 4 ** (n - 1)):
                b_ok = False
                notes.append(f"{name} a={a} b={b}: {4**n/dl**2:.1f},{mid},{t},{up}")
            for (_, size, cnt) in classes:
                k = size // 2
                if k % 2 == 1 and cnt != 1 << (n - 1):
                    k_ok = False
                if k == 2 and cnt != 1 << (n - 2):
                    k_ok = False
    tight.append(f"{name}:4^n/d^2={4**n/dl**2:g}")
exp.check("P6", b_ok, "violations: " + ("; ".join(notes[:5]) if notes else "none"))
exp.log(" ".join(tight))
exp.check("P7", k_ok, "every k-odd class 2^(n-1), every k=2 class 2^(n-2)")

# ---------------------------------------------------------------- P8 affine
exp.section("P8  affine permutation")
vals = {T[("affine5", a, b)] for a in range(32) for b in range(32)}
exp.check("P8", vals == {1}, f"term counts {sorted(vals)}")

# ---------------------------------------------------------------- controls
exp.section("C1  Schroedinger-direction mutant")
failed, where = False, ""
for name in ("gold3_gf32", "rand5"):
    perm = fixtures[name]
    ptm = dense_ptm_perm(perm, stacks[5])
    for a in range(1, 32):
        for b in range(32):
            f, _ = formula_terms(perm, a, b, direction="schr")
            if not same(f, dense_terms(ptm[:, a + (b << 5)], 5)):
                failed, where = True, f"{name} a={a} b={b}"
                break
        if failed:
            break
    if failed:
        break
exp.fail_check("C1", failed, f"first disagreement {where or 'none'}")

exp.section("C2  naive linear extrapolation")
perm = fixtures["rand5"]
failed, where = False, ""
for a in range(1, 32):
    for b in range(32):
        ws = int((np.abs(fwht(np.where(parity(perm & b) == 1, -1.0, 1.0))) > 0.5).sum())
        if ws != T[("rand5", a, b)]:
            failed, where = True, f"a={a} b={b}: naive {ws} vs {T[('rand5', a, b)]}"
            break
    if failed:
        break
exp.fail_check("C2", failed, where or "naive agreed everywhere")

exp.section("C3  dropped sign")
failed, where = False, ""
for a in range(1, 32):
    for b in range(32):
        if parity(np.array([a & b]))[0] == 0:
            continue
        f1, _ = formula_terms(perm, a, b)
        f2, _ = formula_terms(perm, a, b, sign=False)
        if f1.keys() != f2.keys():
            failed, where = True, f"a={a} b={b}"
            break
    if failed:
        break
exp.fail_check("C3", failed, where or "sign never mattered")

exp.section("C4  unsquared mutant bound")
dl = delta_of(fixtures["gold3_gf32"])
mn = min(T[("gold3_gf32", a, b)] for a in range(1, 32) for b in range(32))
exp.fail_check("C4", mn < 4 ** 5 / dl, f"min T = {mn} vs mutant bound 4^n/delta = {4**5/dl:g}")

# ---------------------------------------------------------------- informational: F10 fixture
exp.section("info  F10 fixture (Toffoli modexp N=5 a=2 n_exp=1), derived counts only")
from toffoli_arith import ToffoliModExp
me = ToffoliModExp(5, 2, n_exp=1)
qc = me.build()
p14 = walsh.classical_permutation(qc).astype(np.int64)
for lab, (a, b) in [("X_x0", (1 << me.x[0], 0)), ("Y_x0", (1 << me.x[0], 1 << me.x[0]))]:
    f, classes = formula_terms(p14, a, b)
    ks = sorted({s // 2 for (_, s, _) in classes})
    exp.log(f"{lab}: T = {len(f)}  (R_a = {len(classes)}, pair counts k in {ks[:8]}"
            f"{'...' if len(ks) > 8 else ''}); F10 recorded only 'hit cap' at 60000")

exp.finish()
