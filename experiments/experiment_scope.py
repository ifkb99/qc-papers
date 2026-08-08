"""Two quick questions that bound the Walsh thesis before it gets built on.

Q1 (scope). C8 was verified for single-qubit Z observables. Where exactly does
    it stop? Prediction: it holds for ANY Z-type observable (the pullback of a
    diagonal operator under a permutation is diagonal) and fails for X/Y-type
    (the pullback is not diagonal, so it has no Walsh description).

Q2 (C5). With the theta=pi bug fixed, does coefficient truncation still produce
    estimates that violate |<O>| <= 1? That was the only truncation claim left
    standing, and its earlier evidence used wrong <O> values.
"""
from __future__ import annotations
import numpy as np
from circuits import ripple_adder
from toffoli_arith import ToffoliModExp
from modexp import ModExp
from pps import propagate
import walsh, statevec as sv


def z_pullback_sparsity(qc, zmask, perm=None):
    """Walsh sparsity for a general Z-type observable Z^zmask."""
    if perm is None:
        perm = (walsh.classical_permutation(qc) if qc.is_classical()
                else walsh.permutation_via_statevector(qc))
    n = qc.n
    par = np.bitwise_count(perm & zmask) & 1
    c = walsh.wht(np.where(par == 1, -1.0, 1.0)) / (1 << n)
    return c


print("=" * 76)
print("Q1  SCOPE: which observables does the Walsh identity cover?")
print("=" * 76)
me = ToffoliModExp(5, 2, n_exp=1)
qc = me.build()
perm = walsh.classical_permutation(qc)
print(f"  Toffoli modexp N=5 a=2, {qc.n} qubits, {len(qc.gates)} gates\n")
print(f"  {'observable':>22} {'type':>6} {'walsh':>7} {'pps':>7} {'nonZ':>6} {'match':>6}")

cases = [("Z_x0", (0, 1 << me.x[0])),
         ("Z_x1", (0, 1 << me.x[1])),
         ("Z_x0 Z_x1", (0, (1 << me.x[0]) | (1 << me.x[1]))),
         ("Z_x0 Z_x1 Z_x2", (0, (1 << me.x[0]) | (1 << me.x[1]) | (1 << me.x[2]))),
         ("X_x0", (1 << me.x[0], 0)),
         ("Y_x0", (1 << me.x[0], 1 << me.x[0]))]

for label, P in cases:
    x, z = P
    kind = "Z" if x == 0 else ("X" if z == 0 else "Y")
    # X/Y-type observables blow past any cap; a small one suffices to show that
    # the pullback leaves the diagonal, which is the whole point.
    cap = 2_000_000 if kind == "Z" else 60_000
    r = propagate(qc, {P: 1.0}, delta=0.0, max_terms=cap)
    nonz = sum(1 for (xx, _) in r.final_terms if xx != 0)
    if kind == "Z":
        c = z_pullback_sparsity(qc, z, perm)
        ws = int((np.abs(c) > 1e-9).sum())
        match = "YES" if ws == r.n_terms[-1] else "no"
    else:
        ws, match = -1, "n/a"
    wtxt = "-" if ws < 0 else str(ws)
    tag = " (hit cap)" if r.hit_cap else ""
    print(f"  {label:>22} {kind:>6} {wtxt:>7} {r.n_terms[-1]:7d} {nonz:6d} "
          f"{match:>6}{tag}", flush=True)

print("""
  Z-type: pullback stays diagonal, Walsh predicts the count exactly.
  X/Y-type: pullback leaves the diagonal (nonZ > 0), so there is no Walsh
  description -- the identity is specific to computational-basis observables.
  That is the honest scope: it covers what you actually measure in Shor
  (bits of the output register), and nothing more.""")

print()
print("=" * 76)
print("Q2  C5: does truncation still violate |<O>| <= 1 after the bugfix?")
print("=" * 76)
print(f"  {'circuit':>26} {'delta':>8} {'N_max':>8} {'<O>':>10} {'true':>5} "
      f"{'|err|':>9} {'admissible':>11}")
targets = [("Toffoli modexp N=5", ToffoliModExp(5, 2, n_exp=2)),
           ("Fourier modexp N=5", ModExp(5, 2, n_exp=2)),
           ("Fourier modexp N=7", ModExp(7, 3, n_exp=2))]
violations = 0
total = 0
for label, m in targets:
    c = m.build(); q = m.x[0]
    j, _ = sv.peak(sv.run(c))
    true = (-1.0) ** ((j >> q) & 1)
    for d in (0.0, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1):
        r = propagate(c, {(0, 1 << q): 1.0}, delta=d, max_terms=2_000_000)
        adm = abs(r.expectation) <= 1 + 1e-9
        total += 1
        violations += (not adm)
        print(f"  {label:>26} {d:8.0e} {r.n_max:8d} {r.expectation:+10.5f} "
              f"{true:+5.0f} {abs(r.expectation-true):9.2e} "
              f"{'yes' if adm else 'VIOLATED':>11}", flush=True)
    print()
print(f"  {violations}/{total} truncated runs produced inadmissible estimates "
      f"(|<O>| > 1).")
print("""  The check costs nothing and is one-sided: it can prove a run invalid,
  never valid. Worth reporting as a practitioner's guard regardless of how the
  Walsh thesis develops.""")
