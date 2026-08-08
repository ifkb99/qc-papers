"""Engine tests for lab/ -- every extracted helper is pinned to a number
already logged in NOTES.md/CLAIMS.md before the extraction, so a refactor
that changes behaviour fails against the historical record, not against
itself. Part of the correctness gate (suite 6 of 7).
"""
from __future__ import annotations
import os

import numpy as np

os.environ["LAB_NO_CACHE"] = "1"          # engine tests never touch the cache

from lab import (Experiment, rank_kernel, structures, find_structures,
                 check_structure, coset_split, quadrant_counts, slice_fn,
                 order, v2_split, ord2, carmichael, bit_table, tile,
                 build_modexp, verify_correctness, WRAPS,
                 random_boolean, random_table, planted_structure,
                 fn_spectrum, fn_support, support, sparsity)
import lab.measure as measure


def t(name, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    assert ok, name


print("[1] modarith vs known values")
t("order(6,7)=2", order(6, 7) == 2)
t("order(3,7)=6", order(3, 7) == 6)
t("order(7,15)=4", order(7, 15) == 4)
t("v2_split(12)=(2,3)", v2_split(12) == (2, 3))
t("v2_split(16)=(4,1)", v2_split(16) == (4, 1))
t("ord2(9)=6", ord2(9) == 6)
t("ord2(1)=0", ord2(1) == 0)
t("carmichael: 15->4, 21->6, 51->16, 143->60",
  [carmichael(x) for x in (15, 21, 51, 143)] == [4, 6, 16, 60])
t("bit_table(7,3,0) = [1,1,0,0,0,1]",
  bit_table(7, 3, 0).tolist() == [1, 1, 0, 0, 0, 1])
t("tile length and content", tile(np.array([1, 0, 1]), 4).tolist()
  == [1, 0, 1] * 5 + [1])

print("[2] gf2: rank/kernel against experiment_linstruct logged values")
# modexp N=5 a=2 n_exp=1 (q=14): |supp|=3086, rank 13, w = {b3, anc12}
me = build_modexp(N=5, a=2, n_exp=1)
qc = me.build()
zs_modexp = support(qc, me.x[0])
t("N=5 modexp |supp| = 3086", zs_modexp.size == 3086)
rank, kern = rank_kernel(zs_modexp.tolist(), qc.n)
t("rank 13, defect 1", (rank, len(kern)) == (13, 1))
w0 = (1 << me.b[me.m - 1]) | (1 << me.anc)
t("kernel element = b_msb^anc", kern[0] == w0)
ev, od = coset_split(zs_modexp, w0)
t("support entirely in <z,w>=0", od == 0 and ev == 3086)

print("[3] gf2: affine finder on planted structures (C33 machinery)")
g, w, eps = planted_structure(12, rng=1, eps=1)
zs, st = find_structures(g, 12)
t("planted affine found", st == [(w, 1)])
g, w, eps = planted_structure(12, rng=2, eps=0)
zs, st = find_structures(g, 12)
t("planted linear (dead var) found", st == [(w, 0)])
g = random_boolean(12, rng=3)
zs, st = find_structures(g, 12)
t("random function: no structure (must-fail control)", st == [])

print("[4] gf2: function-level all-ones structures (C33 pinned)")
g = tile(bit_table(7, 3, 0), 12)
zs, st = find_structures(g, 12)
t("r=6 t=12: density exactly 1/2", zs.size == 2048)
t("r=6 t=12: w=all-ones AFFINE", st == [((1 << 12) - 1, 1)])
g13 = tile(bit_table(7, 3, 0), 13)
_, st13 = find_structures(g13, 13)
t("r=6 t=13: none (odd t)", st13 == [])
g3 = tile(bit_table(7, 2, 0), 12)
_, st3 = find_structures(g3, 12)
t("r=3 t=12: w=all-ones LINEAR", st3 == [((1 << 12) - 1, 0)])

print("[5] slices and quadrants against experiment_resid2 logged values")
me4 = build_modexp(N=5, a=2, n_exp=2, wraps=("msb_t0_anc",))
qc4 = me4.build()
zs4 = support(qc4, me4.x[0])
t("v4 |supp| = 23464", zs4.size == 23464)
q = quadrant_counts(zs4, me4.b[me4.m - 1], me4.anc)
t("v4 quadrant (msb=1,anc=0) EMPTY", q[(1, 0)] == 0)
t("v4 quadrants sum to |supp|", sum(q.values()) == 23464)
import walsh as _walsh
g4 = ((_walsh.classical_permutation(qc4) >> me4.x[0]) & 1).astype(np.int64)
s1 = slice_fn(g4, qc4.n, me4.t[0], 1)
w14 = 1 << (me4.b[me4.m - 1])          # msb < t0, position unchanged in slice
t("v4 t0=1 slice structure = msb alone (pinned)",
  check_structure(s1, w14, 0))
t("slice_fn shape", s1.size == 1 << (qc4.n - 1))

print("[6] variants: registry, correctness gate, C35 constancy")
t("known wraps registered",
  set(WRAPS) >= {"msb_t0_anc", "msb_t0_t1", "msb_t1_anc"})
t("v4 computes a^e mod N", verify_correctness(me4))
t("baseline computes a^e mod N",
  verify_correctness(build_modexp(N=5, a=2, n_exp=2)))
try:
    build_modexp(N=5, a=2, n_exp=1, wraps=("nope",))
    t("unknown wrap rejected", False)
except ValueError:
    t("unknown wrap rejected", True)
s2 = sparsity(build_modexp(N=7, a=6, n_exp=2, wraps=("msb_t0_anc",)).build(),
              build_modexp(N=7, a=6, n_exp=2).x[0])
s3 = sparsity(build_modexp(N=7, a=6, n_exp=3, wraps=("msb_t0_anc",)).build(),
              build_modexp(N=7, a=6, n_exp=3).x[0])
t("C35: v4 constancy 23488 at n_exp=2,3", s2 == 23488 and s3 == 23488)

print("[7] measure: cache round-trip (content-addressed)")
os.environ["LAB_NO_CACHE"] = "0"
key_path = measure.CACHE_DIR / f"supp_{measure._key(qc, me.x[0])}.npy"
if key_path.exists():
    key_path.unlink()
zs_a = support(qc, me.x[0])
t("cache file created", key_path.exists())
zs_b = support(qc, me.x[0])
t("cache round-trip identical", np.array_equal(zs_a, zs_b))
t("cached == uncached", np.array_equal(zs_a, zs_modexp))
os.environ["LAB_NO_CACHE"] = "1"

print("[8] nulls: trap defaults")
h = random_table(3, rng=0)
t("random_table rejects constants (r=3)",
  all(0 < int(random_table(3, rng=s).sum()) < 3 for s in range(20)))
t("fn_spectrum normalised (Parseval)",
  abs((fn_spectrum(random_boolean(10, rng=4)) ** 2).sum() - 1) < 1e-9)

print("[9] harness: protocol enforcement")
exp = Experiment("selftest", exit_on_fail=False)
exp.predict("P1", "trivially true")
exp.check("P1", True)
exp.predict("LATE", "declared after measuring")     # must warn
exp.check("LATE", True)
ok = exp.finish()
t("passing experiment reports ok", ok)
t("derive-then-test violation warned",
  any("AFTER measurement" in w for w in exp.warnings))
t("missing must-fail control warned",
  any("must-fail" in w for w in exp.warnings))
exp2 = Experiment("selftest2", exit_on_fail=False)
exp2.must_fail("C1", "control")
exp2.fail_check("C1", False)                        # control did NOT fail
t("vacuous control detected", not exp2.finish())

print("\nALL TESTS PASSED")
