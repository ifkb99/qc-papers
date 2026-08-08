"""TODO 11.1 -- AFFINE structures: the half of the structure test the rank
test cannot see.

A linear structure (g(y^w) = g(y)) confines the Walsh support to the
hyperplane <z,w> = 0; the existing rank test catches it. An AFFINE structure

    g(y ^ w) = g(y) ^ 1        for all y

confines the support to the COSET <z,w> = 1 instead -- which can still span
the full space, so the rank test is blind to it. Both cap density at 1/2.

COUNTING LEMMA (no measurement needed). Density > 1/2 implies NO linear or
affine structure exists. The broken reduction variants sit at 0.716-0.721, so
for them the answer to "is it another structure at higher order" is NO, by
counting alone. Running the finder on them is a tool-sanity check, not a
question. The genuine open target is the FUNCTION level, where the real r=6
table sits at exactly 1/2 on even t (NOTES SS I) -- prime suspect.

PREDICTIONS, WRITTEN BEFORE MEASURING. Derived from the construction:

  (i)   e ^ (2^t - 1) = 2^t - 1 - e, so complementing every exponent bit acts
        on c = e mod r as c -> (rho - 1 - c) mod r, where rho = 2^t mod r.
  (ii)  For N prime and r even, a^(r/2) = -1 (the only square root of 1
        besides 1), and for odd N, bit0(N - x) = 1 ^ bit0(x): the table is
        ANTIPODAL, h[c + r/2] = 1 ^ h[c].
  (iii) When additionally h[-c] = h[c] (i.e. bit_j(x) = bit_j(x^-1 mod N) on
        <a> -- an instance property), composing (i)+(ii) gives an exact
        structure at w = ALL-ONES whenever rho - 1 = -c - shift lands on the
        antipode or on the identity.

  P1  N=7 a=3 bit0 (r=6): h = [1,1,0,0,0,1] is antipodal AND symmetric.
      Even t has rho = 4 = r/2 + 1, so y -> y ^ 11...1 maps c -> 3 - c and
      h[3-c] = 1 ^ h[c]:  w = all-ones is an AFFINE structure (c=1).
      Odd t (rho = 2): no structure, density 1. This is the exact mechanism
      of the 0.500-on-even-t observation in SS I.
  P2  N=7 a=2 bit0 (r=3): h = [1,0,0], symmetric, r odd so no antipode.
      Even t has rho = 1, map is c -> -c, h[-c] = h[c]:  w = all-ones is a
      LINEAR structure (c=0). Same density-1/2 signature, different type.
  P3  NEGATIVE control, N=11 a=2 bit0 (r=10): at t = 0 mod 4, rho = 6 =
      r/2 + 1, the mapping is right BUT bit0 is not inversion-symmetric mod
      11 (2^-1 = 6 keeps bit0, 4^-1 = 3 flips it). Predict NO structure and
      density != 1/2 -- matches the 1.000 rows in SS I.
  P4  Tool sanity, forced by counting: v4/v5 broken variants (density .716)
      -> empty kernel; v0 baseline (density .473 > 1/4) -> exactly the known
      linear w = b_msb ^ anc and nothing affine; random f -> empty; a planted
      g(y) = f(y') ^ y_top -> w = top bit, c = 1.

If P1-P3 land, the function-level half-density is EXPLAINED and gets a
criterion (prime N, even r, inversion-symmetric bit table, t in the right
residue class). The circuit-level 0.716 is then definitively NOT a structure
of this kind and the hunt moves to the zero set (experiment_resid1).
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
from circuits import Circuit
import walsh

TOL = 1e-12


def gf2_rank_and_kernel(vecs, n):
    basis = []
    for v in vecs:
        cur = int(v)
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    rank = len(basis)
    piv = {b.bit_length() - 1: b for b in basis}
    kernel = []
    for free in range(n):
        if free in piv:
            continue
        w = 1 << free
        for p in sorted(piv):
            if (piv[p] & w).bit_count() % 2:
                w ^= (1 << p)
        kernel.append(w)
    return rank, kernel


def affine_structures(c, n):
    """All (w, eps) with <z,w> = eps for every support z.

    Translate the support by one member z0; the kernel of the difference set
    gives every w constant on the support, and eps = <z0,w> classifies it:
    eps = 0 -> linear structure, eps = 1 -> affine structure. The rank test
    alone only sees the eps = 0 kind.
    """
    zs = np.nonzero(np.abs(c) > TOL)[0]
    if zs.size == 0:
        return zs, 0, []
    z0 = int(zs[0])
    diffs = (zs ^ z0).tolist()
    rank, kernel = gf2_rank_and_kernel(diffs, n)
    out = [(w, (z0 & w).bit_count() & 1) for w in kernel]
    return zs, rank, out


def check_pointwise(g, w, eps):
    idx = np.arange(g.size, dtype=np.int64)
    return bool(np.array_equal(g[idx ^ w], g ^ eps))


def table_function(N, a, bit, t):
    r, v = 1, a % N
    while v != 1:
        v = (v * a) % N
        r += 1
    h = np.empty(r, dtype=np.int64)
    x = 1
    for cc in range(r):
        h[cc] = (x >> bit) & 1
        x = (x * a) % N
    return np.tile(h, (1 << t) // r + 1)[: 1 << t], r


def report(label, g, n, expect=None):
    chi = np.where(g == 1, -1.0, 1.0)
    c = walsh.wht(chi) / (1 << n)
    zs, rank, structs = affine_structures(c, n)
    dens = zs.size / (1 << n)
    print(f"  {label}")
    print(f"    n={n}  |supp|={zs.size}  density={dens:.6f}  diff-rank={rank}")
    if not structs:
        print("    no linear or affine structure (kernel empty)")
    for w, eps in structs:
        kind = "AFFINE (eps=1)" if eps else "linear (eps=0)"
        ok = check_pointwise(g, w, eps)
        wtxt = "all-ones" if w == (1 << n) - 1 else f"{w:#x}"
        print(f"    w = {wtxt:>10}  {kind}  pointwise g(y^w)=g(y)^eps: {ok}")
    if expect is not None:
        print(f"    expected: {expect}")
    return zs, structs


print("=" * 84)
print("P4a  Controls first: random (must find none), planted (must find w_top)")
print("=" * 84)
rng = np.random.default_rng(1)
g = rng.integers(0, 2, size=1 << 12).astype(np.int64)
report("random Boolean function, n=12", g, 12, "nothing")
print()
f = rng.integers(0, 2, size=1 << 11).astype(np.int64)
idx = np.arange(1 << 12, dtype=np.int64)
g = f[idx & ((1 << 11) - 1)] ^ ((idx >> 11) & 1)
report("planted g(y) = f(y') ^ y_top, n=12", g, 12, "w = 0x800, AFFINE")
print()

print("=" * 84)
print("P1/P2  Function level: the exact-1/2 rows of SS I")
print("=" * 84)
for N, a, bit, ts, expect in (
    (7, 3, 0, (12, 13, 14), "even t: w=all-ones AFFINE; odd t: none"),
    (7, 2, 0, (12, 13, 14), "even t: w=all-ones LINEAR; odd t: none"),
):
    for t in ts:
        g, r = table_function(N, a, bit, t)
        report(f"idealised bit fn N={N} a={a} bit{bit} r={r}, t={t}", g, t, expect)
        print()

print("=" * 84)
print("P3  Negative control: right mapping, broken symmetry")
print("=" * 84)
g, r = table_function(11, 2, 0, 12)
report(f"idealised bit fn N=11 a=2 bit0 r={r}, t=12  (t=0 mod 4)", g, 12,
       "nothing, density != 1/2")
print()

print("=" * 84)
print("P4b  Circuit level: forced by counting, run as tool sanity")
print("=" * 84)


class Variant(ToffoliModExp):
    def __init__(self, *a, mode="v0", **k):
        self.mode = mode
        super().__init__(*a, **k)

    def _wrap(self, qc):
        msb = self.b[self.m - 1]
        if self.mode == "v4":
            qc.toffoli(msb, self.t[0], self.anc)
        elif self.mode == "v5":
            qc.toffoli(msb, self.t[0], self.t[1])

    def cc_add_mod(self, qc, c1, c2, c):
        self._wrap(qc)
        super().cc_add_mod(qc, c1, c2, c)
        self._wrap(qc)


for mode, expect in (("v0", "one linear w (b_msb^anc), nothing affine"),
                     ("v4", "empty kernel (density > 1/2 forbids)"),
                     ("v5", "empty kernel (density > 1/2 forbids)")):
    me = Variant(N=5, a=2, n_exp=2, mode=mode)
    qc = me.build()
    perm = walsh.classical_permutation(qc)
    g = ((perm >> me.x[0]) & 1).astype(np.int64)
    report(f"modexp N=5 a=2 n_exp=2, variant {mode}", g, qc.n, expect)
    print()
