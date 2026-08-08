"""Claims regression -- cheap executable reproductions of the headline
results in CLAIMS.md, pinned to the exact logged numbers. Suite 7 of the
correctness gate.

Until this file existed, only the *infrastructure* had regression protection;
the science itself did not. Each block names the claim it re-verifies, uses a
small instance that runs in seconds, and includes the control where the
original had one. A future refactor that silently breaks a result fails here
loudly.

Not covered here: C25/C26 (the cryptanalytic-bound import -- kept out of the
automated gate by request; `experiments/experiment_crypto.py` remains the
record), C19 (a citation, not a computation), and anything needing q > 17.
"""
from __future__ import annotations
import os

import numpy as np

os.environ.setdefault("LAB_NO_CACHE", "1")   # regression must recompute

from circuits import ripple_adder
from lab import (build_modexp, support, sparsity, peak_pps, fn_support,
                 find_structures, quadrant_counts, coset_split,
                 bit_table, tile, order, v2_split)
import walsh


def t(name, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    assert ok, name


print("[C8] PPS term count = Walsh sparsity, identical support")
me = build_modexp(N=5, a=2, n_exp=1)
qc = me.build()
zs_w = support(qc, me.x[0])
res = peak_pps(qc, 1 << me.x[0])
zs_p = np.sort(np.array(list(res.final_terms), dtype=np.int64))
t("walsh sparsity = 3086", zs_w.size == 3086)
t("pps final = 3086", len(res.final_terms) == 3086)
t("supports identical (not just counts)", np.array_equal(zs_w, zs_p))

print("[C10] adder collapse is affineness: low bit sparsity 1")
add, lay = ripple_adder(4)
t("4-bit adder Z_b0 sparsity = 1", sparsity(add, lay["b"][0]) == 1)
t("4-bit adder Z_b2 sparsity = 10 (control: higher bit not affine)",
  sparsity(add, lay["b"][2]) == 10)

print("[C15/C18-lite] constancy for beta=1, growth for beta>1 (fix N, vary a)")
t("N=7 a=6 (r=2): 15549 at n_exp=2 and 3",
  sparsity(build_modexp(N=7, a=6, n_exp=2).build(), 8) == 15549
  and sparsity(build_modexp(N=7, a=6, n_exp=3).build(), 8) == 15549)
t("N=7 a=3 (r=6) control: 15539 -> 30712",
  sparsity(build_modexp(N=7, a=3, n_exp=2).build(), 8) == 15539
  and sparsity(build_modexp(N=7, a=3, n_exp=3).build(), 8) == 30712)

print("[C21] onset at n_exp = v2(r)+1: alpha=2 grows once, then locks")
t("N=5 a=2 (r=4): 15493 -> 32143 -> 32143 at n_exp=2,3,4",
  [sparsity(build_modexp(N=5, a=2, n_exp=k).build(), 8) for k in (2, 3, 4)]
  == [15493, 32143, 32143])

print("[C24] support confined to z_I in {0, 1_I} (halves pinned)")
me = build_modexp(N=7, a=6, n_exp=3)
alpha, _ = v2_split(order(6, 7))
zs = support(me.build(), me.x[0])
mask = 0
for eq in me.exp[alpha:]:
    mask |= 1 << eq
zi = zs & mask
n0 = int(np.count_nonzero(zi == 0))
n1 = int(np.count_nonzero(zi == mask))
t("halves 7770/7779, other = 0",
  (n0, n1) == (7770, 7779) and n0 + n1 == zs.size)

print("[C30] the 1/2 ceiling is the linear structure w = b_msb^anc")
me = build_modexp(N=5, a=2, n_exp=1)
qc = me.build()
zs = support(qc, me.x[0])
w = (1 << me.b[me.m - 1]) | (1 << me.anc)
ev, od = coset_split(zs, w)
t("support in the hyperplane <z,w>=0", od == 0)
g = ((walsh.classical_permutation(qc) >> me.x[0]) & 1).astype(np.int64)
idx = np.arange(g.size, dtype=np.int64)
t("g(y^w) = g(y) pointwise", bool(np.array_equal(g[idx ^ w], g)))

print("[C32/C34] one nonlinear monomial: 23464, empty (msb=1,anc=0) quadrant")
me4 = build_modexp(N=5, a=2, n_exp=2, wraps=("msb_t0_anc",))
zs4 = support(me4.build(), me4.x[0])
t("v4 support = 23464", zs4.size == 23464)
quad = quadrant_counts(zs4, me4.b[me4.m - 1], me4.anc)
t("quadrant (1,0) exactly empty", quad[(1, 0)] == 0)
t("density above 1/2 (linear structure destroyed)",
  zs4.size / (1 << me4.n_qubits) > 0.5)

print("[C35] breaking the structure keeps C15 constancy (23488)")
t("v4 N=7 a=6: 23488 at n_exp=2 and 3",
  sparsity(build_modexp(N=7, a=6, n_exp=2, wraps=("msb_t0_anc",)).build(), 8)
  == 23488
  and sparsity(build_modexp(N=7, a=6, n_exp=3, wraps=("msb_t0_anc",)).build(),
               8) == 23488)

print("[C33] function-level all-ones structures, with negative control")
_, st = find_structures(tile(bit_table(7, 3, 0), 12), 12)
t("r=6 even t: all-ones AFFINE", st == [((1 << 12) - 1, 1)])
_, st = find_structures(tile(bit_table(7, 2, 0), 12), 12)
t("r=3 even t: all-ones LINEAR", st == [((1 << 12) - 1, 0)])
zs, st = find_structures(tile(bit_table(11, 2, 0), 12), 12)
t("r=10 control: no structure, density 1", st == [] and zs.size == 4096)

print("[F12-lite] function-level dichotomy")
t("N=15 a=7 (r=4): sparsity 4 at t=12 and t=16",
  fn_support(tile(bit_table(15, 7, 0), 12)).size == 4
  and fn_support(tile(bit_table(15, 7, 0), 16)).size == 4)

print("\nALL TESTS PASSED")
