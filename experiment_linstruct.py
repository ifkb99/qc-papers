"""The unexplained factor of two: do modexp pullbacks have a LINEAR STRUCTURE?

Two independent observations never accounted for:
  * C7: circuit Walsh density converges to 0.498 -- approaching 1/2 FROM BELOW
    and never exceeding it -- where random Boolean functions give 1.000.
  * Step 6: the real modexp table at r=6 gives density exactly 0.500 on even t,
    where a random table of the same period gives 1.000.

Approaching 1/2 from below without crossing is the signature of a LINEAR
CONSTRAINT on the support. If every z in supp(g^) satisfies <z,w> = 0 for some
fixed w != 0, the support lies in a hyperplane of size 2^(n-1), capping density
at exactly 1/2.

Dually, g^(z) = 0 whenever <z,w> = 1 is equivalent to

    g(y XOR w) = g(y)   for all y,

i.e. **w is a linear structure of g** in the cryptanalytic sense -- a nonzero
translation the function is invariant under. Functions with linear structures
are considered weak precisely because that invariance is exploitable, so this
would also feed straight back into the C12 bridge.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  H1  The GF(2) rank of the support of a modexp pullback is n-1, not n. That
      forces density <= 1/2 and gives exactly one nontrivial w.
  H2  Random Boolean functions have full rank n (no constraint, density -> 1).
  H3  w is a MULTI-qubit pattern, not a single dead qubit -- the earlier
      live-variable check found every qubit live, which rules out a dead
      variable but not a multi-bit linear structure.
  H4  Direct verification: g(y XOR w) = g(y) pointwise.

If H1-H4 land, the density-1/2 ceiling is explained exactly, and modexp bit
functions are shown to carry a specific cryptanalytic weakness.
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
from circuits import ripple_adder
import walsh

TOL = 1e-12


def gf2_rank_and_kernel(vecs, n):
    """Rank over GF(2) of the given bitmask vectors, plus a basis of the
    orthogonal complement {w : <z,w> = 0 for all z}."""
    basis = []
    for v in vecs:
        cur = int(v)
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    rank = len(basis)
    # kernel of the map w -> (<z,w>)_z  == orthogonal complement of span(vecs)
    piv = {}
    for b in basis:
        piv[b.bit_length() - 1] = b
    kernel = []
    for free in range(n):
        if free in piv:
            continue
        w = 1 << free
        for p in sorted(piv):
            if (bin(piv[p] & w).count("1")) % 2:
                w ^= (1 << p)
        kernel.append(w)
    return rank, kernel


def analyse(label, c, n, gfun=None, reg=None):
    zs = np.nonzero(np.abs(c) > TOL)[0]
    rank, kernel = gf2_rank_and_kernel(zs.tolist(), n)
    dens = zs.size / (1 << n)
    print(f"  {label}")
    print(f"    n={n}  |support|={zs.size}  density={dens:.6f}  "
          f"GF(2) rank={rank}  (n-rank={n-rank})")
    if kernel:
        w = kernel[0]
        bits = [i for i in range(n) if (w >> i) & 1]
        tags = [f"{reg.get(i,'?')}{i}" for i in bits] if reg else bits
        print(f"    linear structure w = {w} -> qubits {bits}")
        if reg:
            print(f"                              = {tags}")
        ok = all((int(z) & w).bit_count() % 2 == 0 for z in zs.tolist())
        print(f"    every support z satisfies <z,w>=0 : {ok}")
        if gfun is not None:
            idx = np.arange(gfun.size, dtype=np.int64)
            inv = bool(np.array_equal(gfun, gfun[idx ^ w]))
            print(f"    H4  g(y XOR w) == g(y) pointwise  : {inv}")
        if len(kernel) > 1:
            print(f"    (kernel dimension {len(kernel)} -> density <= 2^-{len(kernel)})")
    else:
        print("    no linear structure (full rank)")
    return rank, kernel


print("=" * 80)
print("H1/H3/H4  Modular exponentiation pullbacks")
print("=" * 80)
for N, a, ne in ((5, 2, 1), (7, 6, 2), (7, 3, 2), (15, 7, 1)):
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    qc = me.build()
    perm = walsh.classical_permutation(qc)
    g = ((perm >> me.x[0]) & 1).astype(np.int64)
    c = walsh.pullback_coefficients(qc, me.x[0])
    reg = {}
    for q_ in me.b: reg[q_] = "b"
    for q_ in me.t: reg[q_] = "t"
    for q_ in me.x: reg[q_] = "x"
    reg[me.c0] = "c0"; reg[me.anc] = "anc"
    for q_ in me.exp: reg[q_] = "e"
    analyse(f"modexp N={N} a={a} n_exp={ne}, obs Z_x0", c, qc.n, g, reg)
    print()

print("=" * 80)
print("H2  Controls")
print("=" * 80)
rng = np.random.default_rng(0)
for n in (12, 14):
    gg = rng.integers(0, 2, size=1 << n).astype(np.int64)
    c = walsh.wht(np.where(gg == 1, -1.0, 1.0)) / (1 << n)
    analyse(f"random Boolean function, n={n}", c, n, gg)
    print()

add, lay = ripple_adder(4)
c = walsh.pullback_coefficients(add, lay["b"][2])
perm = walsh.classical_permutation(add)
g = ((perm >> lay["b"][2]) & 1).astype(np.int64)
analyse("ripple adder, obs Z_b2", c, add.n, g)
