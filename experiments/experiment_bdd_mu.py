"""Out of sample: is the exponent-level ROBDD width of compiled modexp governed by mu(beta), the order of 2 in (Z/beta)^x/{+-1}?

Context. run_v1 of experiment_bdd_beta_v1 (board run R5e2b4f1518e74ec6)
found exponent-level widths equal to 2^h (k-h+1) for beta = 3 and to
2^h (F_(k-h+3) - 1) for beta = 5. derivation_v1.md (same attempt directory)
then proved, for any group and any schedule with a_i = a_(i-mu)^(+-1), the
bound |S_j| <= 2|S_(j-1)| - |S_(j-1-mu)|, and applied it with C101's
full-space block inverse to get width(e_k) <= L(k) = 2^k for k < h and
min(2^k, 2^h D_mu(k-h)) otherwise. Here D_mu(j) = 2^j for j <= mu and
D_mu(j) = 2 D_mu(j-1) - D_mu(j-1-mu). None of the fixtures below has been
run in any form.

Object, order and reducer are exactly those of experiment_bdd_beta_v1:
x0 through ToffoliModExp(N, a, n_exp=t).build(), full space; exponent bits
e_0..e_{t-1} first, then work MSB-first interleaved; plain ROBDD via
lab/bdd_count.py.

Fixtures (all n = 5, q_w = 19, t = 1..9, q <= 28):
  N=19 a=2   r=18 beta=9  alpha=1  mu=3 via c_(k+3) = c_k^-1
  N=29 a=2   r=28 beta=7  alpha=2  mu=3 via c_(k+3) = c_k (period only)
  N=31 a=3   r=30 beta=15 alpha=1  mu=4 via c_(k+4) = c_k (period only)
  N=23 a=2   r=11 beta=11 alpha=0  mu=5 via c_(k+5) = c_k^-1
  N=25 a=2   r=20 beta=5  alpha=2  mu=2 via c_(k+2) = c_k^-1

PREDICTIONS, WRITTEN BEFORE MEASURING.
  M1  every exponent-level width <= L(k) with h = alpha, all fixtures, all t.  [proved]
  M2  equality: every exponent-level width at t = 9 equals L(k).
      [conjecture generalizing run_v1's equality to new beta; can fail]
  M3  both mu = 3 fixtures (beta = 9 via inverse, beta = 7 via period) have the
      same tail width sequence D_3(j) = 1,2,4,8,15,28,52,96, scaled by 2^h. [from M2]
  M4  Walsh support ratio W(9)/W(5) in [12, 20] for every fixture (doubling per
      exponent bit for beta > 1, as in C15 and run_v1).           [empirical]

  C1  MUST FAIL. Tribonacci compression must vanish where mu = 5: for N=23 the
      width at k=6 must exceed D_3(6) = 52.
  C2  MUST FAIL. Fibonacci compression must vanish where mu = 4: for N=31 some
      width at k <= 8 must exceed 2^h D_2(k-h).
  C3  MUST FAIL. A mutant recursion that ignores the relation (D(j) = 2^j) must
      disagree with the measured widths somewhere for N=19.

Run (from research/):
  mkdir -p out/bdd_mu && \
  PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=. uv run --no-project \
    --python 3.12 --with 'numpy<2.5' --with 'cupy-cuda13x>=14.1.1' python -u \
    experiments/experiment_bdd_mu.py \
    > out/bdd_mu/mu_run.log 2>&1
"""
from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

import numpy as np

import accel
import walsh
from lab.bdd_count import output_bit_table, robdd_size
from lab.harness import Experiment
from toffoli_arith import ToffoliModExp

OUT = Path("out/bdd_mu")
exp = Experiment("experiment_bdd_mu_v1", doc=__doc__)
exp.predict("M1", "exponent-level width <= L(k) for all fixtures and t")
exp.predict("M2", "exponent-level widths at t=9 equal L(k)")
exp.predict("M3", "beta=9 and beta=7 (mu=3) share tail sequence 1,2,4,8,15,28,52,96 scaled by 2^h")
exp.predict("M4", "Walsh W(9)/W(5) in [12,20] for every fixture")
exp.must_fail("C1", "N=23 (mu=5) width at k=6 > 52")
exp.must_fail("C2", "N=31 (mu=4) some width exceeds the mu=2 bound")
exp.must_fail("C3", "D(j)=2^j mutant disagrees with N=19 widths")

FIXTURES = [(19, 2, 9, 1, 3), (29, 2, 7, 2, 3), (31, 3, 15, 1, 4), (23, 2, 11, 0, 5), (25, 2, 5, 2, 2)]
TMAX = 9


def D(mu: int, j: int) -> int:
    d = [2 ** i for i in range(mu + 1)]
    while len(d) <= j:
        i = len(d)
        d.append(2 * d[i - 1] - d[i - 1 - mu])
    return d[j]


def L(k: int, h: int, mu: int) -> int:
    return 2 ** k if k < h else min(2 ** k, 2 ** h * D(mu, k - h))


def order_and_alpha(N: int, a: int):
    r, v = 1, a % N
    while v != 1:
        v, r = v * a % N, r + 1
    beta, alpha = r, 0
    while beta % 2 == 0:
        beta, alpha = beta // 2, alpha + 1
    return r, beta, alpha


def mu_of(beta: int) -> int:
    o = 1
    while pow(2, o, beta) != 1:
        o += 1
    return next((k for k in range(1, o + 1) if pow(2, k, beta) == beta - 1), o)


def work_order(me):
    order = [me.c0, me.anc, me.t[me.n], me.b[me.n]]
    for i in range(me.n - 1, -1, -1):
        order += [me.x[i], me.t[i], me.b[i]]
    return order


def permutation(qc):
    if qc.n >= 20 and accel.HAVE_GPU:
        return accel.classical_permutation(qc)
    return walsh.classical_permutation(qc)


def walsh_count(qc, j, perm):
    if qc.n >= 20 and accel.HAVE_GPU:
        return int(accel.pullback_stats(qc, j, perm=perm)["count"])
    return int(np.count_nonzero(np.abs(walsh.pullback_coefficients(qc, j, perm=perm)) > 1e-12))


rows: list[dict] = []
exp.section("fixture parameters recomputed from N, a")
for N, a, beta, h, mu in FIXTURES:
    r, b, al = order_and_alpha(N, a)
    ok = (b, al, mu_of(b)) == (beta, h, mu)
    exp.log(f"N={N} a={a}: r={r} beta={b} alpha={al} mu={mu_of(b)} (declared {beta},{h},{mu}) {ok}")
    assert ok

exp.section("M1-M4  sweep")
m1 = m2 = m3 = m4 = True
last: dict[tuple, dict] = {}
walsh_by: dict[tuple, dict[int, int]] = {}
for N, a, beta, h, mu in FIXTURES:
    walsh_by[(N, a)] = {}
    for t in range(1, TMAX + 1):
        t0 = time.time()
        me = ToffoliModExp(N, a, n_exp=t)
        qc = me.build()
        j = me.x[0]
        order = list(me.exp) + work_order(me)
        perm = permutation(qc)
        w = walsh_count(qc, j, perm)
        tab = output_bit_table(perm, j)
        del perm
        gc.collect()
        B, levels = robdd_size(tab, qc.n, order, per_level=True)
        del tab
        gc.collect()
        widths = levels[:t]
        bound = [L(k, h, mu) for k in range(t)]
        m1 &= all(x <= y for x, y in zip(widths, bound))
        walsh_by[(N, a)][t] = w
        row = dict(N=N, a=a, beta=beta, alpha=h, mu=mu, t=t, q=qc.n, walsh=w, robdd=B,
                   exp_levels=widths, bound=bound, seconds=round(time.time() - t0, 1))
        rows.append(row)
        last[(N, a)] = row
        exp.log(row)
        sys.stdout.flush()
    m2 &= last[(N, a)]["exp_levels"] == last[(N, a)]["bound"]
    ratio = walsh_by[(N, a)][9] / walsh_by[(N, a)][5]
    m4 &= 12 <= ratio <= 20
    exp.log(f"N={N} a={a}: widths {last[(N, a)]['exp_levels']} bound {last[(N, a)]['bound']} "
            f"Walsh W9/W5={ratio:.2f}")

tail9 = [x // 2 for x in last[(19, 2)]["exp_levels"][1:]]
tail7 = [x // 4 for x in last[(29, 2)]["exp_levels"][2:]]
target = [1, 2, 4, 8, 15, 28, 52, 96]
m3 &= tail9 == target[:len(tail9)] and tail7 == target[:len(tail7)]
exp.log(f"beta=9 tail/2^h {tail9}; beta=7 tail/2^h {tail7}")
exp.check("M1", m1, "all widths within the proved bound")
exp.check("M2", m2, "t=9 widths equal the bound for every fixture")
exp.check("M3", m3, "mu=3 fixtures share the tribonacci-type tail sequence")
exp.check("M4", m4, "Walsh W9/W5 in [12,20] for every fixture")

exp.section("C1-C3  controls")
w23 = last[(23, 2)]["exp_levels"]
exp.fail_check("C1", w23[6] > D(3, 6), f"N=23 width at k=6 = {w23[6]} vs D_3(6) = {D(3, 6)}")
w31 = last[(31, 3)]["exp_levels"]
fib = [L(k, 1, 2) for k in range(len(w31))]
exp.fail_check("C2", any(x > y for x, y in zip(w31, fib)), f"N=31 widths {w31} vs mu=2 bound {fib}")
w19 = last[(19, 2)]["exp_levels"]
exp.fail_check("C3", w19 != [2 ** k for k in range(len(w19))], f"N=19 widths {w19} vs 2^k")

exp.finish(report_path=OUT / "mu_report_v1.json", rows=rows,
           metadata=dict(task="T58ddf9e611124d65", attempt="A54e702253b2c46ca", gpu=accel.HAVE_GPU))
