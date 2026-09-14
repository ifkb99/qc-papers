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

print("\n[D] exact |+> contraction only after a qubit's final reverse use")
qc = Circuit(2).cnot(0, 1).cnot(0, 1)
traced = propagate_perm(qc, 2, trace_plus=[0])
check("reused control: defer contraction until both CNOTs are propagated",
      traced.final_terms == {2: 1.0} and traced.expectation == 1.0
      and traced.trace_events == [(2, 1, 1, 1)])
premature = propagate_perm(Circuit(2).cnot(0, 1), 2).final_terms
check("control: pruning after only the final CNOT gives the WRONG zero operator",
      not {z: c for z, c in premature.items() if not (z & 1)})
check("unused |+> qubit in the observable contracts to zero",
      propagate_perm(Circuit(2), 1, trace_plus=[0]).final_terms == {})
check("unused |+> qubit outside the observable leaves it unchanged",
      propagate_perm(Circuit(2), 2, trace_plus=[0]).final_terms == {2: 1.0})
split = Circuit(3).toffoli(0, 1, 2)
traced = propagate_perm(split, 4, trace_plus=[0])
check("peak includes the within-gate support before contraction",
      traced.n_max == 4 and traced.trace_events == [(1, 1, 4, 2)])
capped = propagate_perm(split, 4, trace_plus=[0], max_terms=2)
check("contraction must not hide a pre-contraction term-cap violation",
      capped.hit_cap and capped.n_max == 4 and not capped.trace_events)

rng = np.random.default_rng(20260909)
max_error = 0.0
max_state_error = 0.0
nontrivial = 0
import statevec
for _ in range(24):
    n = 5
    qc = Circuit(n)
    for _ in range(14):
        a, b, c = map(int, rng.choice(n, 3, replace=False))
        kind = int(rng.integers(3))
        if kind == 0:
            qc.x(a)
        elif kind == 1:
            qc.cnot(a, b)
        else:
            qc.toffoli(a, b, c)
    tq = int(rng.integers(n))
    qs = list(map(int, rng.choice(n, int(rng.integers(1, n)), replace=False)))
    mask = sum(1 << q for q in qs)
    exact = walsh.pullback_coefficients(qc, tq)
    ids = np.arange(1 << n)
    exact[(ids & mask) != 0] = 0
    traced = propagate_perm(qc, 1 << tq, trace_plus=qs)
    got = np.zeros(1 << n)
    for z, value in traced.final_terms.items():
        got[z] = value
    max_error = max(max_error, float(np.max(np.abs(got - exact))))
    nontrivial += 0 < np.count_nonzero(exact) < (1 << n)
    psi = np.where((ids & ~mask) == 0, 2. ** (-len(qs) / 2), 0.)
    out = statevec.run(qc, psi.astype(complex))
    expectation = float(np.dot(np.abs(out) ** 2, 1 - 2 * ((ids >> tq) & 1)))
    max_state_error = max(max_state_error, abs(expectation - traced.expectation))
check("24 random circuits: the entire reduced coefficient vector matches Walsh",
      max_error < 1e-12 and nontrivial > 0,
      f"max error {max_error:.2e}; {nontrivial} nonzero examples")
check("independent statevector check with actual |+>/|0> inputs",
      max_state_error < 1e-10, f"max error {max_state_error:.2e}")
try:
    propagate_perm(Circuit(2), 1, trace_plus=[2])
except ValueError:
    check("reject an out-of-range contracted qubit", True)
else:
    check("reject an out-of-range contracted qubit", False)

print("\n" + ("ALL TESTS PASSED" if not FAILED else f"FAILURES: {FAILED}"))
raise SystemExit(bool(FAILED))
