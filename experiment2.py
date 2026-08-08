"""PARTIALLY RETRACTED. Section 0 (gate verification) and section 1 (the adder
collapses to one term) still stand. Sections 2 and 3 -- the "Shor sandwich" and
the claim that the QFT is the bottleneck -- are RETRACTED: the toy is not a
proxy for Shor (see experiment3.py header) and the numbers came from a
bit-reversed QFT.

The mechanism paragraph below is also superseded. The adder collapses to one
term because its low output bit is XOR-AFFINE (Walsh sparsity 1), not because
the circuit is a permutation -- permutation-ness only guarantees the pullback
stays diagonal, it says nothing about how many terms that costs. See NOTES.md F9.

---

Why arithmetic alone is easy for PPS, and what Shor's structure does to that.

First run showed something I did not predict: a pure ripple-carry adder with a
Z observable collapses to ONE Pauli term under PPS, exactly, no matter how
small delta gets. The reason is structural:

  A classical reversible circuit is a permutation matrix. Conjugating a
  diagonal operator by a permutation gives another diagonal operator. Z-type
  Pauli strings span the diagonals, so Z-type strings map to Z-type strings.
  Every intermediate branch into an X/Y-type string must cancel exactly.

That protection is what Shor's algorithm destroys: the QFT turns Z-type into
X-type, so the observable arrives at the arithmetic block in the one form that
does branch. This script tests that.
"""
from __future__ import annotations
import numpy as np
from circuits import Circuit, ripple_adder
from pps import propagate, exact_expectation

FAIL = []
def check(name, cond, extra=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' -- ' + extra) if extra and not cond else ''}")
    if not cond: FAIL.append(name)


print("=" * 76)
print("0. VERIFY the new gates")
print("=" * 76)
qc = Circuit(2); qc.cphase(0, 1, np.pi / 3)
tgt = np.diag([1, 1, 1, np.exp(1j * np.pi / 3)]).astype(complex)
U = qc.to_unitary(); ph = U[0, 0] / tgt[0, 0]
check("controlled-phase decomposition", np.allclose(U, ph * tgt, atol=1e-9))

m = 3
qc = Circuit(m); qc.qft(list(range(m)))
F = np.array([[np.exp(2j * np.pi * i * j / 2 ** m) for j in range(2 ** m)]
              for i in range(2 ** m)]) / np.sqrt(2 ** m)
U = qc.to_unitary()
# QFT up to the usual bit-reversal and a global phase
# with the swap network (default) the QFT is the plain DFT, sign +1
ph = U[0, 0] / F[0, 0]
check("QFT == plain DFT (swap network on)", np.allclose(U, ph * F, atol=1e-9))

qc = Circuit(m); qc.qft(list(range(m))); qc.qft(list(range(m)), inverse=True)
check("QFT then inverse QFT == identity",
      np.allclose(qc.to_unitary(), qc.to_unitary()[0, 0] * np.eye(2 ** m), atol=1e-9))


def zfrac(terms):
    tot = sum(v * v for v in terms.values())
    z = sum(v * v for (x, _), v in terms.items() if x == 0)
    return z / tot if tot > 0 else 0.0


def run(label, qc, obs, deltas=(0.0, 1e-6, 1e-4, 1e-2)):
    print(f"\n  {label}")
    print(f"    {'delta':>9} {'N_max':>8} {'N_final':>9} {'<O>':>13} {'err vs exact':>13}")
    try:
        truth = exact_expectation(qc, obs)
    except Exception:
        truth = None
    for d in deltas:
        r = propagate(qc, obs, delta=d, max_terms=600_000)
        err = f"{abs(r.expectation - truth):13.2e}" if truth is not None else " " * 13
        cap = " (CAP)" if r.hit_cap else ""
        print(f"    {d:9.0e} {r.n_max:8d} {r.n_terms[-1] if r.n_terms else 0:9d} "
              f"{r.expectation:13.8f} {err}{cap}")
    return truth


print()
print("=" * 76)
print("1. MECHANISM: does a permutation circuit preserve Z-type strings?")
print("=" * 76)
add, lay = ripple_adder(3)
obs_z = {(0, 1 << lay["b"][0]): 1.0}
r = propagate(add, obs_z, delta=0.0, max_terms=600_000)
final = {k: v for k, v in zip([None], [None])}  # placeholder
res = propagate(add, obs_z, delta=0.0, max_terms=600_000)
print(f"  adder(3-bit), observable Z_b0")
print(f"    non-Clifford (T) gates in circuit : {add.n_nonclifford()}")
print(f"    peak Pauli terms during propagation: {res.n_max}")
print(f"    final Pauli terms                 : {res.n_terms[-1]}")
print(f"    <O> = {res.expectation:.10f}  (exact {exact_expectation(add, obs_z):.10f})")
print("    -> intermediate branching is real, but it all cancels at the end.")

print("\n  Same circuit, X-type observable (X_b0) -- no diagonal protection:")
obs_x = {(1 << lay["b"][0], 0): 1.0}
res_x = propagate(add, obs_x, delta=0.0, max_terms=600_000)
print(f"    peak Pauli terms: {res_x.n_max}   final terms: {res_x.n_terms[-1]}")
print(f"    <O> = {res_x.expectation:.10f}  (exact {exact_expectation(add, obs_x):.10f})")

print()
print("=" * 76)
print("2. THE SHOR SANDWICH:  H layer  ->  arithmetic  ->  inverse QFT")
print("=" * 76)
nb = 3
_, lay = ripple_adder(nb)
areg = lay["a"]

# (a) arithmetic only
c_arith, _ = ripple_adder(nb)

# (b) H layer then arithmetic
c_h = Circuit(2 * nb + 2)
for q in areg: c_h.h(q)
c_h.gates += ripple_adder(nb)[0].gates

# (c) full sandwich: H, arithmetic, inverse QFT on the a register
c_shor = Circuit(2 * nb + 2)
for q in areg: c_shor.h(q)
c_shor.gates += ripple_adder(nb)[0].gates
c_shor.qft(areg, inverse=True)

obs = {(0, 1 << areg[0]): 1.0}
for label, qc in (("(a) arithmetic only            ", c_arith),
                  ("(b) H layer + arithmetic       ", c_h),
                  ("(c) H + arithmetic + inverse QFT", c_shor)):
    run(label + f"  [T gates={qc.n_nonclifford()}]", qc, obs)

print()
print("=" * 76)
print("3. WHERE THE PROTECTION BREAKS: Z-weight of the observable during (c)")
print("=" * 76)
for label, qc in (("arithmetic only", c_arith), ("full sandwich", c_shor)):
    terms = dict(obs)
    zf = []
    from pauli import commutes, i_sigma_p
    for sigma, theta in reversed(qc.gates):
        c, s = np.cos(theta), np.sin(theta)
        if abs(s) < 1e-12: continue
        new = {}
        for P, coeff in terms.items():
            if commutes(P, sigma):
                new[P] = new.get(P, 0.0) + coeff; continue
            if abs(c) > 1e-12: new[P] = new.get(P, 0.0) + c * coeff
            Q, sg = i_sigma_p(sigma, P)
            new[Q] = new.get(Q, 0.0) + s * sg * coeff
        terms = {k: v for k, v in new.items() if abs(v) > 1e-13}
        zf.append(zfrac(terms))
    zf = np.array(zf)
    print(f"  {label:18s} fraction of weight on Z-type strings: "
          f"start={zf[0]:.3f} min={zf.min():.3f} end={zf[-1]:.3f}  "
          f"(mean {zf.mean():.3f})")
print("""
  Z-type weight staying at 1.0 means the observable never leaves the diagonal
  subspace and PPS is effectively free. Weight leaking away means the branching
  is real and truncation has to start making choices.""")

print("\n" + ("ALL CHECKS PASSED" if not FAIL else f"FAILURES: {FAIL}"))
