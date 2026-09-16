"""Does the shrinking ROBDD growth step of the u_a pullback persist at n = 8 (29 qubits)?

Context. experiment_heisenberg_bdd run_v2 (board run Rdb87a87850f64ccc, log
out/heisenberg_bdd/run_v2.log, report out/heisenberg_bdd/report_v1.json)
measured, for ToffoliModExp(N, 2, n_exp=1).u_a(ctrl, 2), qubit x0, full dirty
space, MSB-first interleaved order, median plain ROBDD log2 steps 2.239,
2.162, 2.037, 1.629 (n = 4..7) against median Walsh log2 steps near 3, and
Walsh/ROBDD ratios R = 4.38, 9.63, 17.85, 35.28, 91.79 (n = 3..7). Its
registered P5a failed at n = 4 only (Walsh step 3.375). Nothing at n = 8 has
been observed. Same object, order, reducer and Walsh route as that run; the
helpers below are copied verbatim so that frozen source stays untouched.

Two hypotheses, no commitment between them:
  H-table: B grows like 2^n poly(n) (table scale in N), so the step tends
           to 1 plus a shrinking polynomial term.
  H-square: B grows like 2^(2n) with small-n transients, so the step
           returns towards 2.

PREDICTIONS, WRITTEN BEFORE MEASURING (n = 8 fixtures N in {129, 163, 197,
251}, all coprime to 2; n = 7 fixtures repeated from run_v2).
  Q1  median plain ROBDD log2 step n=7 -> 8 lies in [1.0, 2.0].   [empirical]
  Q2  median Walsh log2 step n=7 -> 8 lies in [2.9, 3.1].        [empirical]
  Q3  2 <= R_8 / R_7 <= 4.                                       [empirical]
  Q4  the n = 7 rows reproduce run_v2 exactly (walsh and both ROBDD counts
      for all eight N).                                          [determinism]

  C1  MUST FAIL. Random reversible circuit with u_a's gate multiset (seed
      20260915 + n, N = 2^(n-1) + 1, as in run_v2) at n = 6 and 7: its ratio
      growth R_7 / R_6 must stay below 2 (the arithmetic effect vanishes).

Run (from research/):
  PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 \
      --with 'numpy<2.5' --with 'cupy-cuda13x>=14.1.1' python -u \
      -m experiments.experiment_heisenberg_bdd_n8 > out/heisenberg_bdd/n8_run_v1.log 2>&1
Needs a CUDA card for q >= 20 (29-qubit replay and exact Walsh count) and
roughly 20 GiB of host memory for the 2^29-entry reduction.
"""
from __future__ import annotations

import gc
import math
import random
import statistics
import sys
import time
from pathlib import Path

import numpy as np

import accel
import walsh
from circuits import Circuit
from lab.bdd_count import output_bit_table, robdd_size
from lab.harness import Experiment
from toffoli_arith import ToffoliModExp

OUT = Path("out/heisenberg_bdd")
RUN_V2_N7 = {  # (walsh, robdd_plain, robdd_complement) from run_v2.log
    65: (16666608, 184735, 155052), 73: (16677172, 190388, 161166),
    81: (16677746, 184369, 160059), 89: (16676342, 186502, 160567),
    97: (16666468, 170020, 141076), 105: (16672812, 178842, 150277),
    113: (16662972, 172938, 138358), 121: (16658080, 173566, 135394),
}
RANDOM_SEED = 20260915

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)
exp.predict("Q1", "median ROBDD log2 step n=7->8 in [1.0, 2.0]")
exp.predict("Q2", "median Walsh log2 step n=7->8 in [2.9, 3.1]")
exp.predict("Q3", "2 <= R_8/R_7 <= 4")
exp.predict("Q4", "n=7 rows reproduce run_v2 exactly")
exp.must_fail("C1", "random same-multiset circuit: R_7/R_6 < 2")


def permutation(qc: Circuit) -> np.ndarray:
    if qc.n >= 20 and accel.HAVE_GPU:
        return accel.classical_permutation(qc)
    return walsh.classical_permutation(qc)


def walsh_count(qc: Circuit, j: int, perm: np.ndarray) -> int:
    """Final PPS support (C8): exact integer route on the GPU for q >= 20,
    walsh.py's float transform with its default tolerance below that."""
    if qc.n >= 20 and accel.HAVE_GPU:
        return int(accel.pullback_stats(qc, j, perm=perm)["count"])
    return int(np.count_nonzero(np.abs(walsh.pullback_coefficients(qc, j, perm=perm)) > 1e-12))


def ua_fixture(N: int):
    me = ToffoliModExp(N, 2, n_exp=1)
    qc = me.u_a(me.exp[0], 2)
    order = [me.exp[0], me.c0, me.anc, me.t[me.n], me.b[me.n]]
    for i in range(me.n - 1, -1, -1):
        order += [me.x[i], me.t[i], me.b[i]]
    return me, qc, order


def measure(qc: Circuit, j: int, order) -> tuple[int, int, int]:
    perm = permutation(qc)
    w = walsh_count(qc, j, perm)
    tab = output_bit_table(perm, j)
    del perm
    gc.collect()
    bp = robdd_size(tab, qc.n, order)
    bc = robdd_size(tab, qc.n, order, "complement")
    del tab
    gc.collect()
    return w, bp, bc


if not accel.HAVE_GPU:
    raise SystemExit("this experiment needs the CUDA backend for q >= 20")

rows: list[dict] = []
per_n: dict[int, list[tuple[int, int]]] = {7: [], 8: []}
q4 = True
exp.section("Q1-Q4  u_a(ctrl, 2), x0, full dirty space, n = 7 and 8")
for n, Ns in ((7, sorted(RUN_V2_N7)), (8, [129, 163, 197, 251])):
    for N in Ns:
        t0 = time.time()
        me, qc, order = ua_fixture(N)
        assert me.n == n
        w, bp, bc = measure(qc, me.x[0], order)
        if n == 7:
            q4 &= (w, bp, bc) == RUN_V2_N7[N]
        per_n[n].append((w, bp))
        row = dict(fixture="u_a", n=n, N=N, q=qc.n, gates=len(qc.logical),
                   walsh=w, robdd_plain=bp, robdd_complement=bc,
                   seconds=round(time.time() - t0, 1))
        rows.append(row)
        exp.log(row)
        sys.stdout.flush()

W = {n: statistics.median(w for w, _ in v) for n, v in per_n.items()}
B = {n: statistics.median(b for _, b in v) for n, v in per_n.items()}
R = {n: W[n] / B[n] for n in per_n}
for n in per_n:
    exp.log(f"n={n}: median W={W[n]} (log2 {math.log2(W[n]):.3f}) "
            f"median B={B[n]} (log2 {math.log2(B[n]):.3f}) R={R[n]:.2f}")
bstep = math.log2(B[8]) - math.log2(B[7])
wstep = math.log2(W[8]) - math.log2(W[7])
exp.check("Q1", 1.0 <= bstep <= 2.0, f"ROBDD log2 step {bstep:.3f}")
exp.check("Q2", 2.9 <= wstep <= 3.1, f"Walsh log2 step {wstep:.3f}")
exp.check("Q3", 2 <= R[8] / R[7] <= 4, f"R_8/R_7 = {R[8] / R[7]:.3f}")
exp.check("Q4", q4, "n=7 walsh/plain/complement equal to run_v2 for all eight N")

exp.section("C1  random reversible circuit with u_a's gate multiset, n = 6, 7")
Rr = {}
for n in (6, 7):
    me, qc, order = ua_fixture(2 ** (n - 1) + 1)
    rng = random.Random(RANDOM_SEED + n)
    rq = Circuit(qc.n)
    for op in qc.logical:
        getattr(rq, op[0])(*rng.sample(range(qc.n), len(op) - 1))
    w, bp, bc = measure(rq, me.x[0], order)
    Rr[n] = w / bp
    row = dict(fixture="random_multiset", n=n, q=rq.n, walsh=w, robdd_plain=bp,
               robdd_complement=bc, ratio=Rr[n])
    rows.append(row)
    exp.log(row)
exp.fail_check("C1", Rr[7] / Rr[6] < 2,
               f"random R_7/R_6 = {Rr[7] / Rr[6]:.3f} (u_a R_8/R_7 = {R[8] / R[7]:.3f})")

exp.finish(report_path=OUT / "n8_report_v1.json", rows=rows,
           metadata=dict(task="Ta294b7d88a4c4015", attempt="A46391e1cf2954646",
                         gpu=accel.HAVE_GPU))
