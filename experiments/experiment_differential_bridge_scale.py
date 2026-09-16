"""Does identity (I) of C99 survive at 14 qubits, where F10 only saw a cap hit?

Context. experiment_differential_bridge verified C99's DDT identity for n <= 6
against dense matrices and for 8 qubits against gate-level PPS. Its logged
out-of-sample prediction for F10's Toffoli modexp (N=5, a=2, n_exp=1; 14
qubits) could not be settled by the dictionary propagator: Clifford+T
intermediates exceeded 6M terms (DB). This run settles it with a dense GPU
reference: U and X^a Z^b are materialised as 16384 x 16384 float64 matrices,
multiplied with cuBLAS and projected onto the Pauli basis by a batched Walsh
transform. The reference never forms sigma, D_c or a DDT.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P9   The GPU reference (lab.differential.dense_terms_gpu) equals the formula
       support on every label of the random n=6 permutation (seed 20260914,
       same draw order as the bridge experiment) and of the 6-qubit adder.
       This validates the new reference before it is trusted at 14 qubits.
  P10  F10 fixture, diagonal labels, recorded in F10 before this project
       had C99: Z_x0 -> 3086, Z_x1 -> 2926, Z_x0 Z_x1 -> 2848 terms.
  P11  F10 fixture, off-diagonal labels, predicted in
       out/differential_bridge/run1.log: X_x0 -> 849,836 and
       Y_x0 -> 849,442 terms, with full (c, d) support sets equal to the
       formula's.

  C5   The Schroedinger-direction mutant sigma'(y) = pi(pi^-1(y)^a) must give
       a support set different from the GPU reference for X_x0.

Run:  uv run python -m experiments.experiment_differential_bridge_scale   (from research/)
GPU: one RTX A4500 (device 1), peak ~10 GiB; CPU otherwise.
"""
from __future__ import annotations

import time

import numpy as np

from lab import Experiment
from lab.differential import formula_counts, dense_terms_gpu, formula_terms
from circuits import ripple_adder
from toffoli_arith import ToffoliModExp
import walsh

exp = Experiment("experiment_differential_bridge_scale", doc=__doc__)
exp.predict("P9", "GPU dense reference == formula on every label of rand6 and adder6")
exp.predict("P10", "F10 fixture: Z_x0 3086, Z_x1 2926, Z_x0Z_x1 2848")
exp.predict("P11", "F10 fixture: X_x0 849836, Y_x0 849442, support sets equal formula")
exp.must_fail("C5", "Schroedinger-direction set differs from GPU reference for X_x0")


def same_sets(f, g):
    return f.keys() == g.keys() and all(np.array_equal(f[c], g[c]) for c in f)


def total(t):
    return sum(len(v) for v in t.values())


exp.section("P9  validate the GPU reference at n = 6")
rng = np.random.default_rng(20260914)
rng.permutation(16), rng.permutation(32)           # rand4, rand5 draws
rand6 = rng.permutation(64).astype(np.int64)
add6, _ = ripple_adder(2)
fx = {"rand6": rand6, "adder6": walsh.classical_permutation(add6).astype(np.int64)}
bad = {}
t0 = time.time()
for name, perm in fx.items():
    bad[name] = sum(not same_sets(formula_counts(perm, a, b), dense_terms_gpu(perm, a, b))
                    for a in range(64) for b in range(64))
exp.check("P9", all(v == 0 for v in bad.values()),
          f"mismatching labels {bad} of 4096 each ({time.time()-t0:.0f}s)")

exp.section("P10/P11  F10 fixture, 14 qubits")
me = ToffoliModExp(5, 2, n_exp=1)
perm = walsh.classical_permutation(me.build()).astype(np.int64)
x0, x1 = me.x[0], me.x[1]
diag = {"Z_x0": ((0, 1 << x0), 3086), "Z_x1": ((0, 1 << x1), 2926),
        "Z_x0Z_x1": ((0, (1 << x0) | (1 << x1)), 2848)}
got = []
ok10 = True
for lab_, ((a, b), want) in diag.items():
    t0 = time.time()
    g = dense_terms_gpu(perm, a, b)
    got.append(f"{lab_}={total(g)} ({time.time()-t0:.0f}s)")
    ok10 &= total(g) == want
exp.check("P10", ok10, ", ".join(got))

off = {"X_x0": ((1 << x0, 0), 849_836), "Y_x0": ((1 << x0, 1 << x0), 849_442)}
got, ok11, dense_x = [], True, None
for lab_, ((a, b), want) in off.items():
    t0 = time.time()
    g = dense_terms_gpu(perm, a, b)
    f = formula_counts(perm, a, b)
    eq = same_sets(f, g)
    got.append(f"{lab_}: dense {total(g)}, formula {total(f)}, sets equal {eq} "
               f"({time.time()-t0:.0f}s)")
    ok11 &= total(g) == want and eq
    if lab_ == "X_x0":
        dense_x = g
exp.check("P11", ok11, "; ".join(got))

exp.section("C5  Schroedinger-direction mutant at 14 qubits")
terms, _ = formula_terms(perm, 1 << x0, 0, direction="schr")
mut: dict[int, list[int]] = {}
for (c, d) in terms:
    mut.setdefault(c, []).append(d)
mut = {c: np.sort(np.array(v, dtype=np.int64)) for c, v in mut.items()}
exp.fail_check("C5", not same_sets(mut, dense_x),
               f"mutant {total(mut)} terms over {len(mut)} X-parts vs reference "
               f"{total(dense_x)} over {len(dense_x)}")

exp.finish()
