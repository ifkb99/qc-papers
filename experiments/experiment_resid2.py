"""TODO 11 endgame -- the residue is a CONDITIONAL linear structure.

experiment_resid1 left three facts wanting one mechanism: (i) one injected
monomial msb&t0 gives density ~0.72 whether injected once or twice (v4, v5,
v45), (ii) a second independent monomial (v3x) jumps to ~0.98, (iii) the v4
trajectory climbs 0.7196 -> 0.7417, which is 98.9% of 3/4 at the same size
where the baseline sits at 98.7% of its proved cap 1/2.

READING THE CONSTRUCTION. The t register is restored between reductions
(Cuccaro MAJ/UMA restores the operand; loads are unloaded), so at every wrap
point t holds its INITIAL value. Hence on the half-space {y : t0(y) = 0} the
v4/v5/v45 wraps are literally the identity and the circuit IS the baseline
there. The baseline's linear structure w = b_msb ^ anc (proved, L2) therefore
survives as a CONDITIONAL structure:

    g_v4(y ^ w) = g_v4(y)   for every y with t0(y) = 0.        (*)

Walsh consequence: for any z with <z,w> = 1 the even half cancels pairwise, so
G^(z) is carried by the t0=1 slice alone -- and its zeros there are the zeros
of the 14-bit slice function's own spectrum. If that slice carries a linear/
affine structure of its own, the odd coset is capped at half density and the
TOTAL density is capped at (1 + 1/2)/2 = 3/4. v3x has wraps inert only on a
QUARTER space, which forces no exact zeros at all -- consistent with 0.98.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  T1  (*) holds pointwise: 0 violations on the t0=0 half. MUST-FAIL control:
      violations on the t0=1 half are nonzero (else the structure would be
      total and density <= 1/2, contradicting 0.716).
  T2  Zero localisation: support restricted to the coset <z,w> = 0 is near
      full (generic ~0.99+); the deficit lives almost entirely in <z,w> = 1.
  T3  The t0=0 slice of g_v4 equals the t0=0 slice of g_v0 exactly.
  T4  The t0=1 slice, as a 14-bit function, itself has a linear or affine
      structure (this is what caps the odd coset at 1/2 and the total at 3/4;
      it is the one link not yet forced). If found, verify pointwise and
      check the support's coset split against it.
  T5  v3x control: no conditional structure at t0=0 (its wrap C fires there),
      and no large coset asymmetry w.r.t. w.
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
import walsh

TOL = 1e-12


class Variant(ToffoliModExp):
    MODES = {"v0": (), "v4": ("A",), "v45": ("A", "B"), "v3x": ("A", "B", "C")}

    def __init__(self, *a, mode="v0", **k):
        self.mode = mode
        super().__init__(*a, **k)

    def _wrap(self, qc):
        msb = self.b[self.m - 1]
        for tag in self.MODES[self.mode]:
            if tag == "A":
                qc.toffoli(msb, self.t[0], self.anc)
            elif tag == "B":
                qc.toffoli(msb, self.t[0], self.t[1])
            elif tag == "C":
                qc.toffoli(msb, self.t[1], self.anc)

    def cc_add_mod(self, qc, c1, c2, c):
        self._wrap(qc)
        super().cc_add_mod(qc, c1, c2, c)
        self._wrap(qc)


def gf2_rank_and_kernel(vecs, n):
    basis = []
    for v in vecs:
        cur = int(v)
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
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
    return len(basis), kernel


def affine_structures(g, n):
    chi = np.where(g == 1, -1.0, 1.0)
    c = walsh.wht(chi) / (1 << n)
    zs = np.nonzero(np.abs(c) > TOL)[0]
    z0 = int(zs[0])
    _, kernel = gf2_rank_and_kernel((zs ^ z0).tolist(), n)
    return zs, [(w, (z0 & w).bit_count() & 1) for w in kernel]


def slice_fn(g, n, bit, val):
    """Restrict g to {y : y_bit = val} as a function of n-1 bits."""
    idx = np.arange(1 << (n - 1), dtype=np.int64)
    lo = idx & ((1 << bit) - 1)
    hi = (idx >> bit) << (bit + 1)
    return g[hi | (np.int64(val) << bit) | lo]


N, a, ne = 5, 2, 2
me = Variant(N=N, a=a, n_exp=ne, mode="v4")
qc = me.build()
n = qc.n
perm = walsh.classical_permutation(qc)
g4 = ((perm >> me.x[0]) & 1).astype(np.int64)
g0 = ((walsh.classical_permutation(Variant(N=N, a=a, n_exp=ne, mode="v0").build())
       >> me.x[0]) & 1).astype(np.int64)

w = (1 << me.b[me.m - 1]) | (1 << me.anc)          # b_msb ^ anc
t0 = me.t[0]
idx = np.arange(1 << n, dtype=np.int64)
even = ((idx >> t0) & 1) == 0

print(f"modexp N={N} a={a} n_exp={ne}, q={n}; w = b{me.b[me.m-1]}^anc{me.anc} "
      f"= {w:#x}; t0 = qubit {t0}")
print()
print("T1  Conditional structure g(y^w) = g(y) on the t0=0 half")
viol_even = int(np.count_nonzero((g4[idx ^ w] != g4) & even))
viol_odd = int(np.count_nonzero((g4[idx ^ w] != g4) & ~even))
print(f"    violations on t0=0 half : {viol_even}   (predict 0)")
print(f"    violations on t0=1 half : {viol_odd}   (must-fail control, predict >0)")

print()
print("T3  t0=0 slice of v4 == t0=0 slice of v0")
s4_0 = slice_fn(g4, n, t0, 0)
s0_0 = slice_fn(g0, n, t0, 0)
print(f"    identical: {bool(np.array_equal(s4_0, s0_0))}   (predict True)")
s4_1 = slice_fn(g4, n, t0, 1)
s0_1 = slice_fn(g0, n, t0, 1)
print(f"    (t0=1 slices identical: {bool(np.array_equal(s4_1, s0_1))} -- "
      f"expect False, the wraps act there)")

print()
print("T2  Where do the zeros live?  support split by <z,w>")
c = walsh.wht(np.where(g4 == 1, -1.0, 1.0)) / (1 << n)
zs = np.nonzero(np.abs(c) > TOL)[0]
par = np.array([(int(z) & w).bit_count() & 1 for z in zs.tolist()])
half = 1 << (n - 1)
d_even = int(np.count_nonzero(par == 0)) / half
d_odd = int(np.count_nonzero(par == 1)) / half
print(f"    density on <z,w>=0 : {d_even:.6f}   (predict ~1)")
print(f"    density on <z,w>=1 : {d_odd:.6f}   (predict = 2*0.716 - even)")

print()
print("T4  Does the t0=1 slice carry its own structure?  (the 3/4 link)")
w14 = ((w & ((1 << t0) - 1)) | ((w >> (t0 + 1)) << t0))   # w in 14-bit coords
zs1, structs = affine_structures(s4_1, n - 1)
print(f"    slice density {zs1.size/(1<<(n-1)):.6f}, "
      f"structures found: {len(structs)}")
i14 = np.arange(1 << (n - 1), dtype=np.int64)
for ww, eps in structs:
    ok = bool(np.array_equal(s4_1[i14 ^ ww], s4_1 ^ eps))
    same = "  == w restricted!" if ww == w14 else ""
    print(f"    w' = {ww:#x}  eps={eps}  pointwise: {ok}{same}")

print()
print("T5  v3x control: quarter-space inertness forces nothing")
g3 = ((walsh.classical_permutation(Variant(N=N, a=a, n_exp=ne, mode="v3x").build())
       >> me.x[0]) & 1).astype(np.int64)
viol_even3 = int(np.count_nonzero((g3[idx ^ w] != g3) & even))
c3 = walsh.wht(np.where(g3 == 1, -1.0, 1.0)) / (1 << n)
zs3 = np.nonzero(np.abs(c3) > TOL)[0]
par3 = np.array([(int(z) & w).bit_count() & 1 for z in zs3.tolist()])
print(f"    violations of (*) on t0=0 half : {viol_even3}   (predict >0)")
print(f"    density on <z,w>=0 : {int(np.count_nonzero(par3==0))/half:.6f}")
print(f"    density on <z,w>=1 : {int(np.count_nonzero(par3==1))/half:.6f}"
      f"   (predict both ~0.98, no asymmetry)")

print()
print("BONUS  same anatomy at N=7 a=3 (the trajectory instance)")
me7 = Variant(N=7, a=3, n_exp=2, mode="v4")
qc7 = me7.build()
n7 = qc7.n
g7 = ((walsh.classical_permutation(qc7) >> me7.x[0]) & 1).astype(np.int64)
w7 = (1 << me7.b[me7.m - 1]) | (1 << me7.anc)
t07 = me7.t[0]
i7 = np.arange(1 << n7, dtype=np.int64)
even7 = ((i7 >> t07) & 1) == 0
print(f"    violations of (*) on t0=0 half : "
      f"{int(np.count_nonzero((g7[i7 ^ w7] != g7) & even7))}")
c7 = walsh.wht(np.where(g7 == 1, -1.0, 1.0)) / (1 << n7)
zs7 = np.nonzero(np.abs(c7) > TOL)[0]
par7 = np.array([(int(z) & w7).bit_count() & 1 for z in zs7.tolist()])
h7 = 1 << (n7 - 1)
print(f"    density on <z,w>=0 : {int(np.count_nonzero(par7==0))/h7:.6f}")
print(f"    density on <z,w>=1 : {int(np.count_nonzero(par7==1))/h7:.6f}")
s7_1 = slice_fn(g7, n7, t07, 1)
zs71, structs7 = affine_structures(s7_1, n7 - 1)
print(f"    t0=1 slice: density {zs71.size/(1<<(n7-1)):.6f}, "
      f"structures {len(structs7)}")
for ww, eps in structs7:
    ok = bool(np.array_equal(s7_1[np.arange(1 << (n7-1), dtype=np.int64) ^ ww],
                             s7_1 ^ eps))
    print(f"    w' = {ww:#x}  eps={eps}  pointwise: {ok}")
