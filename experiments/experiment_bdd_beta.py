"""Does the ROBDD of the modexp pullback key on beta the way PPS (Walsh support) does?

Context. Paper B / C15: for x0 pulled back through ToffoliModExp(N, a,
n_exp=t).build() on the full space, Walsh support is constant in t when
beta = 1 and grows exponentially when beta > 1. TODO13 asks whether a third,
structurally unrelated store keys on the same invariant. C102 (task
Ta294b7d88a4c4015, submitted, unreviewed) introduced the truth-table ROBDD
reducer lab/bdd_count.py with brute-force and dd cross-checks. Board task
T58ddf9e611124d65 (proposal mode).

Variable order, top first: e_0, e_1, ..., e_{t-1}, then the work register
[c0, anc, t_n, b_n, x_{n-1}, t_{n-1}, b_{n-1}, ..., x_0, t_0, b_0].
Plain ROBDD counts (two terminals, no complement edges). The node count at
level v is the number of distinct subfunctions at that cut which depend on v.

DERIVATION (before measuring). The exponent register is only a control, so
f(e, w) = bit x0 of P_t(w), with P_k = M_{c_{k-1}}^{e_{k-1}} ... M_{c_0}^{e_0},
where M_c is the compiled block u_a(., c) on the full dirty work space and
c_k = a^(2^k) mod N. After e_0..e_{k-1} are fixed, the subfunction of the
remaining variables is (e_rest, w) -> f(e_rest, P_k(w)); it is determined by
the permutation P_k. So the width at exponent level k is at most the number
of distinct P_k, and the work part is the shared ROBDD of at most
s_t = #{P_t} functions of the q_w work bits, whose width at work level l is at
most min(s_t 2^l, 2^(2^(q_w - l))). Write WorkCap(q_w, s) for that sum plus
two terminals.

  beta = 1 (r = 2^alpha). c_k = 1 for k >= alpha, and C23: every identity
  block applies the same involution V, so P_t = V^(parity of tail) H with H
  from the head. The subfunction classes at each tail level are (head
  prefix, parity so far). Both parities are reachable at every tail level
  after the first (at level alpha only parity 0 is), and the dependence and
  distinctness of these classes do not depend on k or t, and the work-part
  function set is the same for every t >= alpha + 1. Hence the node count is
  identical at every tail exponent level k = alpha+1..t-1, and B(t) is
  exactly affine in t for t >= alpha + 1.

  beta = 3. After the head, c_{k+1} = c_k^2 = c_k^-1 (c_k^3 = 1), and C101
  proves M_{c^-1} = M_c^-1 as full-space permutations. So for k >= alpha,
  P_k = M_c^s H with H one of at most 2^alpha head products and s an
  alternating signed sum of k - alpha bits, which takes at most k - alpha + 1
  values. With L(k) = 2^k for k < alpha and min(2^k, 2^alpha (k - alpha + 1))
  otherwise, the width at exponent level k is at most L(k), and
  B(t) <= sum_{k<t} L(k) + WorkCap(q_w, L(t)), polynomial in t.

  beta = 5 (N = 11). After the head the schedule cycles c, d, c^-1, d^-1:
  two generators. If M_c and M_d commute on the full work space, the same
  argument bounds widths by O(k^2); if not, no bound follows. Open, no
  commitment between H-comm and H-free.

PREDICTIONS, WRITTEN BEFORE MEASURING.
Fixtures: N=7, a in {6 (r=2, beta=1, alpha=1), 2 (r=3, beta=3, alpha=0),
3 (r=6, beta=3, alpha=1)}, t = 1..16 (q = 13 + t); N=11, a in {10 (r=2,
beta=1), 2 (r=10, beta=5), 3 (r=5, beta=5)}, t = 1..12; N=13, a in {5 (r=4,
beta=1, alpha=2), 2 (r=12, beta=3, alpha=2), 3 (r=3, beta=3, alpha=0)},
t = 1..12 (q = 16 + t). Head length alpha is read off the schedule as the
first k after which the beta-structure holds (listed per fixture below).
  B1  beta = 1: exponent-level counts identical over tail levels
      k = alpha+1..t-1, and B(t+1) - B(t) constant for t >= alpha + 1. [derived]
  B3  beta = 3: every exponent-level count <= L(k), with h the head length
      in place of alpha, and B(t) <= sum_{k<t} L(k) + WorkCap(q_w, L(t)). [derived]
  B5  beta = 5: report whether M_c, M_d commute on the full work space and
      the exponent-level counts; no prediction.                       [open]
  W1  N=7, a=6: Walsh support = 15549 for every t = 2..16 (C15's perm_pps
      record at t = 2, 4, 6, 8; C24 proves independence of tail length). [ledger]
  W2  Block-inverse identity M_{c^-1} = M_c^-1 holds as full work-space
      permutations for (N, c) = (7, 2), (13, 3), (11, 4), (11, 5).     [C101]
  S1  For every beta = 3 fixture, B(t_max) / B(t_max - 4) <= 2, while for
      Walsh the same ratio is >= 8.   [expected from the B3 bound's slow growth
      and C15's Walsh growth; not itself derived]

  C1  MUST FAIL. Polynomial growth must vanish in the Walsh basis: for N=7,
      a=2, Walsh(16) / Walsh(8) >= 2^6.
  C2  MUST FAIL. Corrupting 1% of the N=7, a=2, t=12 truth table (seeded
      positions) must destroy the small ROBDD: size >= 10x unperturbed.
  C3  MUST FAIL. Mutant block-inverse check: M_2 composed with M_2 (not
      M_4) must not be the identity at N=7.

Run (from research/):
  mkdir -p out/bdd_beta && \
  PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=. uv run --no-project \
    --python 3.12 --with 'numpy<2.5' --with 'cupy-cuda13x>=14.1.1' python -u \
    experiments/experiment_bdd_beta.py \
    > out/bdd_beta/run.log 2>&1
"""
from __future__ import annotations

import gc
import math
import sys
import time
from pathlib import Path

import numpy as np

import accel
import walsh
from lab.bdd_count import output_bit_table, robdd_size
from lab.harness import Experiment
from toffoli_arith import ToffoliModExp

OUT = Path("out/bdd_beta")
exp = Experiment("experiment_bdd_beta_v1", doc=__doc__)
exp.predict("B1", "beta=1: tail exponent-level counts identical; B(t+1)-B(t) constant for t>=alpha+1")
exp.predict("B3", "beta=3: exponent-level count <= min(2^k, 2^h(k-h+1)); B(t) <= polynomial bound")
exp.predict("B5", "beta=5: commutation and widths reported, no prediction")
exp.predict("W1", "N=7 a=6 Walsh = 15549 for t=2..16")
exp.predict("W2", "M_{c^-1} = M_c^-1 on full work space for (7,2),(13,3),(11,4),(11,5)")
exp.predict("S1", "beta=3: B(tmax)/B(tmax-4) <= 2 while Walsh ratio >= 8")
exp.must_fail("C1", "N=7 a=2 Walsh(16)/Walsh(8) >= 2^6")
exp.must_fail("C2", "1% corrupted N=7 a=2 t=12 table: ROBDD >= 10x")
exp.must_fail("C3", "M_2 o M_2 is not the identity at N=7")

# (N, a, beta, head length h, t range)
FIXTURES = [
    (7, 6, 1, 1, range(1, 17)), (7, 2, 3, 0, range(1, 17)), (7, 3, 3, 1, range(1, 17)),
    (11, 10, 1, 1, range(1, 13)), (11, 2, 5, 1, range(1, 13)), (11, 3, 5, 0, range(1, 13)),
    (13, 5, 1, 2, range(1, 13)), (13, 2, 3, 2, range(1, 13)), (13, 3, 3, 0, range(1, 13)),
]


def permutation(qc):
    if qc.n >= 20 and accel.HAVE_GPU:
        return accel.classical_permutation(qc)
    return walsh.classical_permutation(qc)


def walsh_count(qc, j, perm):
    if qc.n >= 20 and accel.HAVE_GPU:
        return int(accel.pullback_stats(qc, j, perm=perm)["count"])
    return int(np.count_nonzero(np.abs(walsh.pullback_coefficients(qc, j, perm=perm)) > 1e-12))


def work_order(me):
    order = [me.c0, me.anc, me.t[me.n], me.b[me.n]]
    for i in range(me.n - 1, -1, -1):
        order += [me.x[i], me.t[i], me.b[i]]
    return order


def level_cap(k: int, h: int) -> int:
    return 2 ** k if k < h else min(2 ** k, 2 ** h * (k - h + 1))


def work_cap(qw: int, s: int) -> int:
    return sum(min(s * (1 << l), 2 ** (2 ** (qw - l)) if qw - l <= 5 else s * (1 << l))
               for l in range(qw)) + 2


def block_perm(N: int, c: int) -> np.ndarray:
    """Full work-space permutation of the compiled block u_a(ctrl, c) with ctrl = 1."""
    me = ToffoliModExp(N, 2 if N != 2 else 1, n_exp=1)
    qc = me.u_a(me.exp[0], c)
    perm = walsh.classical_permutation(qc)
    qw = me.exp[0]
    assert qw == qc.n - 1
    mask = (1 << qw) - 1
    return (perm[(1 << qw) + np.arange(1 << qw)] & mask).astype(np.int64)


rows: list[dict] = []

exp.section("W2, C3, B5  compiled block permutations on the full work space")
w2 = True
for N, c in [(7, 2), (13, 3), (11, 4), (11, 5)]:
    M, Mi = block_perm(N, c), block_perm(N, pow(c, -1, N))
    ok = np.array_equal(Mi[M], np.arange(M.size))
    w2 &= ok
    exp.log(f"N={N} c={c}: M_c^-1 o M_c == id: {ok}")
exp.check("W2", w2, "block-inverse identity on all four pairs")
M2 = block_perm(7, 2)
exp.fail_check("C3", not np.array_equal(M2[M2], np.arange(M2.size)), "M_2 o M_2 != id at N=7")
commute = {}
for N, (c, d) in [(11, (4, 5)), (11, (9, 4))]:
    Mc, Md = block_perm(N, c), block_perm(N, d)
    diff = int(np.count_nonzero(Mc[Md] != Md[Mc]))
    commute[(N, c, d)] = diff
    exp.log(f"N={N}: M_{c} M_{d} vs M_{d} M_{c}: {diff} of {Mc.size} work states differ")

exp.section("B1, B3, B5, W1, S1  sweep over t")
b1 = b3 = w1 = s1 = True
series: dict[tuple, dict[int, tuple[int, int]]] = {}
corrupt_input = None
for N, a, beta, h, trange in FIXTURES:
    series[(N, a)] = {}
    for t in trange:
        t0 = time.time()
        me = ToffoliModExp(N, a, n_exp=t)
        qc = me.build()
        j = me.x[0]
        order = list(me.exp) + work_order(me)
        assert sorted(order) == list(range(qc.n))
        perm = permutation(qc)
        w = walsh_count(qc, j, perm)
        tab = output_bit_table(perm, j)
        del perm
        gc.collect()
        B, levels = robdd_size(tab, qc.n, order, per_level=True)
        if (N, a, t) == (7, 2, 12):
            corrupt_input = (tab.copy(), qc.n, order, B)
        del tab
        gc.collect()
        exp_levels = levels[:t]
        qw = qc.n - t
        if beta == 3:
            bound_levels = all(cnt <= level_cap(k, h) for k, cnt in enumerate(exp_levels))
            bound_total = sum(level_cap(k, h) for k in range(t)) + work_cap(qw, level_cap(t, h))
            b3 &= bound_levels and B <= bound_total
        else:
            bound_levels = bound_total = None
        if N == 7 and a == 6 and t >= 2:
            w1 &= w == 15549
        series[(N, a)][t] = (w, B)
        row = dict(N=N, a=a, beta=beta, head=h, t=t, q=qc.n, walsh=w, robdd=B,
                   exp_levels=exp_levels, work_nodes=B - sum(exp_levels),
                   level_bound_ok=bound_levels, total_bound=bound_total,
                   seconds=round(time.time() - t0, 1))
        rows.append(row)
        exp.log(row)
        sys.stdout.flush()
    ts = sorted(series[(N, a)])
    if beta == 1:
        tail_counts = [rows[-1]["exp_levels"][k] for k in range(h + 1, ts[-1])]
        diffs = [series[(N, a)][t + 1][1] - series[(N, a)][t][1] for t in ts if t >= h + 1 and t + 1 in series[(N, a)]]
        ok = len(set(tail_counts)) == 1 and len(set(diffs)) == 1
        b1 &= ok
        exp.log(f"N={N} a={a}: tail level counts {sorted(set(tail_counts))}, B increments {sorted(set(diffs))}")
    if beta == 3:
        tm = ts[-1]
        rb = series[(N, a)][tm][1] / series[(N, a)][tm - 4][1]
        rw = series[(N, a)][tm][0] / series[(N, a)][tm - 4][0]
        s1 &= rb <= 2 and rw >= 8
        exp.log(f"N={N} a={a}: B({tm})/B({tm - 4}) = {rb:.3f}, Walsh ratio {rw:.1f}")
    if beta == 5:
        exp.log(f"N={N} a={a}: B by t {[series[(N, a)][t][1] for t in ts]}")

exp.check("B1", b1, "identical tail level counts and constant increments for every beta=1 fixture")
exp.check("B3", b3, "level and total bounds hold for every beta=3 fixture and t")
exp.check("B5", True, f"recorded; commutation differences {commute}")
exp.check("W1", w1, "N=7 a=6 Walsh == 15549 for t=2..16")
exp.check("S1", s1, "beta=3 ROBDD ratio <= 2 and Walsh ratio >= 8 over the last four t")

exp.section("C1, C2  controls")
wr = series[(7, 2)][16][0] / series[(7, 2)][8][0]
exp.fail_check("C1", wr >= 2 ** 6, f"N=7 a=2 Walsh(16)/Walsh(8) = {wr:.1f}")
tab, n, order, B = corrupt_input
rng = np.random.default_rng(20260915)
idx = rng.choice(tab.size, size=tab.size // 100, replace=False)
tab[idx] ^= 1
Bc = robdd_size(tab, n, order)
exp.fail_check("C2", Bc >= 10 * B, f"corrupted {Bc} vs clean {B} ({Bc / B:.1f}x)")

exp.finish(report_path=OUT / "report_v1.json", rows=rows,
           metadata=dict(task="T58ddf9e611124d65", attempt="A54e702253b2c46ca",
                         commutation=str(commute), gpu=accel.HAVE_GPU))
