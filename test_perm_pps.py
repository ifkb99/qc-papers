"""Verify permutation-native PPS against Walsh and against rotation-level PPS."""
import numpy as np
from circuits import ripple_adder, Circuit
from toffoli_arith import ToffoliModExp
from pps import propagate
from perm_pps import propagate_perm
import walsh

FAILED=[]
def check(name, cond, extra=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' -- '+extra) if extra else ''}")
    if not cond: FAILED.append(name)

print("\n[A] agrees with Walsh (final coefficients) and stays diagonal")
for label, qc, tq in (("3-bit adder Z_b0", *((lambda t: (t[0], t[1]['b'][0]))(ripple_adder(3)))),
                      ("4-bit adder Z_b2", *((lambda t: (t[0], t[1]['b'][2]))(ripple_adder(4))))):
    c = walsh.pullback_coefficients(qc, tq)
    r = propagate_perm(qc, 1 << tq, delta=0.0)
    v = np.zeros(1 << qc.n)
    for z, val in r.final_terms.items(): v[z] = val
    check(f"{label}: coeffs match Walsh", np.max(np.abs(v-c)) < 1e-9,
          f"maxerr {np.max(np.abs(v-c)):.2e}")

print("\n[B] agrees with rotation-level PPS on <O>")
for N,a,ne in ((5,2,1),(7,3,1)):
    me = ToffoliModExp(N=N,a=a,n_exp=ne); qc = me.build(); q = me.x[0]
    rp = propagate_perm(qc, 1<<q, delta=0.0)
    rr = propagate(qc, {(0,1<<q):1.0}, delta=0.0, max_terms=4_000_000)
    check(f"modexp N={N}: <O> matches ({rp.expectation:+.6f})",
          abs(rp.expectation-rr.expectation) < 1e-9,
          f"{rp.expectation} vs {rr.expectation}")
    check(f"modexp N={N}: final count matches", len(rp.final_terms)==rr.n_terms[-1],
          f"{len(rp.final_terms)} vs {rr.n_terms[-1]}")

print("\n[C] PEAK MEMORY: permutation-native vs rotation-level")
print(f"  {'circuit':>20} {'final':>8} {'rot N_max':>10} {'perm N_max':>11} {'saving':>8}")
rows=[]
for nb in (3,4,5):
    add, lay = ripple_adder(nb)
    tq = lay['b'][0]
    rr = propagate(add, {(0,1<<tq):1.0}, delta=0.0, max_terms=4_000_000)
    rp = propagate_perm(add, 1<<tq, delta=0.0)
    print(f"  {f'{nb}-bit adder Z_b0':>20} {len(rp.final_terms):8d} {rr.n_max:10d} "
          f"{rp.n_max:11d} {rr.n_max/max(rp.n_max,1):7.1f}x")
for N,a,ne in ((5,2,1),(7,3,1),(15,7,1)):
    me = ToffoliModExp(N=N,a=a,n_exp=ne); qc = me.build(); q = me.x[0]
    rr = propagate(qc, {(0,1<<q):1.0}, delta=0.0, max_terms=4_000_000)
    rp = propagate_perm(qc, 1<<q, delta=0.0)
    print(f"  {f'modexp N={N}':>20} {len(rp.final_terms):8d} {rr.n_max:10d} "
          f"{rp.n_max:11d} {rr.n_max/max(rp.n_max,1):7.1f}x")

print("\n" + ("ALL TESTS PASSED" if not FAILED else f"FAILURES: {FAILED}"))
