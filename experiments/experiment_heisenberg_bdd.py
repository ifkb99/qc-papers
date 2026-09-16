"""Does a canonical decision diagram store the pulled-back arithmetic observable asymptotically more compactly than PPS (Walsh support)?

Context. C8: PPS storage for Z_j through a basis permutation is the Walsh
sparsity of f(y) = bit j of U(y), on the full space. A reduced ordered BDD of
the SAME f is canonical and exact. TODO54 / task Ta294b7d88a4c4015. Scratch
pilots (board message M741d61dd419947ef) used dd complement-edge counts for
n <= 6; the plain two-terminal counts measured here and every n = 7 value
were never observed before this file was written.

DERIVATIONS (before this run). Cuccaro adder b += t + c0 mod 2^m as built by
ToffoliModExp._add_t_into_b (t, c0 restored); f = top sum bit
s = b_{m-1} ^ t_{m-1} ^ c_{m-1}, carry c_k into position k.
  * Walsh. In +-1 form C_{k+1} = (B_k + T_k + C_k - B_k T_k C_k)/2 with C_k
    free of B_k, T_k, so the four monomial families are disjoint:
    S_{k+1} = 2 + 2 S_k, S_0 = 1, S_k = 3 2^k - 2. B_{m-1} T_{m-1} multiplies
    by a character on fresh variables, so sparsity(f) = 3 2^(m-1) - 2.
  * ROBDD, order c0, b0, t0, ..., b_{m-1}, t_{m-1}. Level c0: 1 node. For
    i < m-1, level b_i: the subfunctions are indexed by the incoming carry
    (2, distinct, both depend on b_i); level t_i: after (c,b_i) = (0,0)/(1,1)
    the carry is fixed (no node), after (0,1)/(1,0) it equals t_i (one shared
    node). Level b_{m-1}: 2 nodes; level t_{m-1}: t ^ parity, 2 nodes; two
    terminals. Plain total 3m + 4. With complemented edges the b_{m-1} pair,
    the t_{m-1} pair and the terminals each merge: 3m + 1.
  * Register-block order b_0..b_{m-1}, t_0..t_{m-1}, c0: after all b the
    subfunction of (t, c0) is "bit m-1 of b + t + c0", a half-length cyclic
    interval shifted by -b; distinct b give distinct non-constant
    subfunctions, so the ROBDD has at least 2^m nodes.
  * Affine pieces. p_i = b_i ^ t_i. {c_k = 1} is the disjoint union over the
    highest generating position i in {-1 (c0 = 1), 0..k-1 (b_i = t_i = 1)} of
    {generate at i, p_l = 1 for i < l < k}: k + 1 affine pieces; {c_k = 0}
    likewise with kills: k + 1 pieces. Splitting s on p_{m-1} gives a
    disjoint affine decomposition with 2m pieces.

PREDICTIONS, WRITTEN BEFORE MEASURING.
  P1  Adder Walsh sparsity = 3 2^(m-1) - 2, m = 1..10.            [derived]
  P2  Adder interleaved ROBDD = 3m + 4 plain, 3m + 1 complement, m = 1..10.
                                                                  [derived]
  P3  Adder register-block ROBDD >= 2^m plain, m = 1..10.         [derived]
  P4  lab.affine_pieces on the adder, no merging, reproduces f exactly with
      <= 2m pieces when a Toffoli branches on op[2] (the propagate wire after
      MAJ's CNOTs). Branching on op[1] gives 2^m.  [the 2m bound is derived;
      both exact counts reproduce the scratch pilot, not a fresh test]
  P5  Walsh vs plain ROBDD of ToffoliModExp(N, 2, n_exp=1).u_a(ctrl, 2),
      qubit x0, full dirty space (q = 3n + 5), order
      [exp, c0, anc, t_n, b_n, x_{n-1}, t_{n-1}, b_{n-1}, ..., x_0, t_0, b_0].
      Fixtures: every odd N in [2^(n-1), 2^n) for n = 3..6; N in
      {65, 73, 81, 89, 97, 105, 113, 121} for n = 7. Medians over N:
        P5a  log2 W_n - log2 W_{n-1} in [2.8, 3.2] for n = 4..7.
        P5b  log2 B_n - log2 B_{n-1} in [1.0, 2.7] for n = 4..7, and
             the n = 7 step is <= 2.2.
        P5c  R_n = W_n / B_n strictly increases over n = 3..7, R_7 >= 3 R_5.
      [empirical, extrapolated from pilot trends; not derived]
  P6  Independent canonical route: dd (0.6.0) gate-by-gate composition gives
      the same final size as the complement-edge reducer, for the adder
      m = 1..8 and u_a at every n = 3, 4 fixture.                 [derived]
  P7  Peaks. dd composition peak / final <= 1.35 for every n = 3, 4, 5
      fixture; PPS (perm_pps) peak / final Walsh >= 2 for every n = 3, 4
      fixture.  [pilot-informed for N = 7, 11, 13, 23; others unobserved]
  P8  Every plain count B satisfies Bc <= B <= 2 (Bc - 1) + 2, where Bc is the
      complement count (each negation class holds at most two plain nodes
      and there are at most two terminals).                       [derived]

  C1  MUST FAIL. Random reversible circuit with u_a's exact gate multiset
      (qubits drawn uniformly, seed 20260915 + n), same q, observable index
      and order, n = 3..6 with N = 2^(n-1) + 1: the growth R_6 / R_3 of its
      Walsh/ROBDD ratio must stay below 3 (the arithmetic effect vanishes).
  C2  MUST FAIL. Adder register-block order: the linear-size effect must
      vanish, i.e. size(m = 10) >= 2^10.
  C3  MUST FAIL. Mutant reducer without redundant-node elimination must
      violate the P2 formula for some m in 2..10.
  C4  MUST FAIL. Mutant adder with the first MAJ Toffoli removed (and its UMA
      partner) must violate the P1 formula for some m in 2..10.

Run (from research/):
  uv run --with dd==0.6.0 python -m experiments.experiment_heisenberg_bdd \
      > out/heisenberg_bdd/run_v1.log 2>&1
Uses accel (GPU) for q >= 20 when available; walsh.py below that.
"""
from __future__ import annotations

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
from lab import affine_pieces as AP
from lab.bdd_count import output_bit_table, robdd_size
from lab.harness import Experiment
from perm_pps import propagate_perm
from toffoli_arith import ToffoliModExp

OUT = Path("out/heisenberg_bdd")
exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P1", "adder Walsh sparsity = 3*2^(m-1)-2, m=1..10")
exp.predict("P2", "adder interleaved ROBDD = 3m+4 plain, 3m+1 complement, m=1..10")
exp.predict("P3", "adder register-block ROBDD >= 2^m, m=1..10")
exp.predict("P4", "affine pieces exact; op[2] branch <= 2m, op[1] branch = 2^m")
exp.predict("P5a", "u_a median Walsh log2 step in [2.8,3.2], n=4..7")
exp.predict("P5b", "u_a median ROBDD log2 step in [1.0,2.7], n=4..7; n=7 step <= 2.2")
exp.predict("P5c", "R_n = W/B strictly increasing n=3..7 and R_7 >= 3 R_5")
exp.predict("P6", "dd composition final size == complement reducer (adder m<=8, u_a n=3,4)")
exp.predict("P7", "dd peak/final <= 1.35 (n=3,4,5); PPS peak/final >= 2 (n=3,4)")
exp.predict("P8", "Bc <= B <= 2(Bc-1)+2 for every u_a fixture")
exp.must_fail("C1", "random same-multiset circuit: R_6/R_3 < 3")
exp.must_fail("C2", "adder register-block order at m=10 has >= 2^10 nodes")
exp.must_fail("C3", "reducer without redundancy elimination violates P2 formula")
exp.must_fail("C4", "adder without first MAJ/UMA Toffoli violates P1 formula")

rows: list[dict] = []


def cuccaro(m: int, drop_first_toffoli: bool = False) -> Circuit:
    """b = 0..m-1, t = m..2m-1, c0 = 2m; the exact MAJ/UMA gates of
    ToffoliModExp._add_t_into_b. The mutant replaces the first MAJ and the
    last UMA by their CNOT parts only."""
    qc = Circuit(2 * m + 1)
    b, t, c0 = list(range(m)), list(range(m, 2 * m)), 2 * m
    if drop_first_toffoli:
        qc.cnot(t[0], b[0]); qc.cnot(t[0], c0)
    else:
        ToffoliModExp._maj(qc, c0, b[0], t[0])
    for i in range(1, m):
        ToffoliModExp._maj(qc, t[i - 1], b[i], t[i])
    for i in range(m - 1, 0, -1):
        ToffoliModExp._uma(qc, t[i - 1], b[i], t[i])
    if drop_first_toffoli:
        qc.cnot(t[0], c0); qc.cnot(c0, b[0])
    else:
        ToffoliModExp._uma(qc, c0, b[0], t[0])
    return qc


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


def quasi_reduced_size(table, n, order) -> int:
    """C3 mutant: merge equal pairs but never eliminate lo == hi nodes."""
    from lab.bdd_count import _arrange
    ids = _arrange(np.asarray(table).astype(np.int8), n, order).astype(np.int64)
    total, nxt = 2, 2
    for _ in range(n):
        pairs = ids.reshape(-1, 2)
        key = pairs[:, 0] * (1 << 32) + pairs[:, 1]
        uniq, inv = np.unique(key, return_inverse=True)
        ids = inv + nxt
        nxt += uniq.size
        total += uniq.size
    return total


def dd_compose(qc: Circuit, j: int, order):
    """Independent canonical route: dd BDD of the pulled-back bit, composed
    gate by gate in reverse; returns (final nodes, peak nodes)."""
    from dd.autoref import BDD
    bdd = BDD()
    bdd.configure(reordering=False)
    bdd.declare(*[f"y{q}" for q in order])
    var = {q: bdd.var(f"y{q}") for q in range(qc.n)}
    f = var[j]
    peak = len(f)
    for op in reversed(qc.logical):
        tq = op[-1]
        if f"y{tq}" not in bdd.support(f):
            continue
        if op[0] == "x":
            sub = ~var[tq]
        elif op[0] == "cnot":
            sub = bdd.apply("xor", var[tq], var[op[1]])
        else:
            sub = bdd.apply("xor", var[tq], var[op[1]] & var[op[2]])
        f = bdd.let({f"y{tq}": sub}, f)
        peak = max(peak, len(f))
    return len(f), peak


def ua_fixture(N: int):
    me = ToffoliModExp(N, 2, n_exp=1)
    qc = me.u_a(me.exp[0], 2)
    order = [me.exp[0], me.c0, me.anc, me.t[me.n], me.b[me.n]]
    for i in range(me.n - 1, -1, -1):
        order += [me.x[i], me.t[i], me.b[i]]
    return me, qc, order


# ---------------------------------------------------------------- adder ----
exp.section("P1-P4, C2-C4  Cuccaro adder top sum bit")
p1 = p2 = p3 = p4 = True
c2_block = None
c3_violated = c4_violated = False
for m in range(1, 11):
    qc = cuccaro(m)
    perm = permutation(qc)
    j, n = m - 1, qc.n
    tab = output_bit_table(perm, j)
    inter = [2 * m] + [q for i in range(m) for q in (i, m + i)]
    block = list(range(2 * m + 1))
    w = walsh_count(qc, j, perm)
    bp, bc = robdd_size(tab, n, inter), robdd_size(tab, n, inter, "complement")
    bb = robdd_size(tab, n, block)
    quasi = quasi_reduced_size(tab, n, inter)
    p1 &= w == 3 * 2 ** (m - 1) - 2
    p2 &= bp == 3 * m + 4 and bc == 3 * m + 1
    p3 &= bb >= 2 ** m
    if m >= 2 and quasi != 3 * m + 4:
        c3_violated = True
    mut = cuccaro(m, drop_first_toffoli=True)
    wm = walsh_count(mut, j, permutation(mut))
    if m >= 2 and wm != 3 * 2 ** (m - 1) - 2:
        c4_violated = True
    row = dict(fixture="adder", m=m, walsh=w, robdd_plain=bp, robdd_complement=bc,
               robdd_block=bb, quasi_reduced=quasi, mutant_walsh=wm)
    if m <= 8:
        sec, _ = AP.pullback(qc, j, "second", merge=False)
        fir, _ = AP.pullback(qc, j, "first", merge=False)
        want = tab.astype(np.int64)
        exact = (np.array_equal(AP.evaluate(sec, n), want)
                 and np.array_equal(AP.evaluate(fir, n), want))
        p4 &= exact and len(sec) <= 2 * m and len(fir) == 2 ** m
        ddf, ddp = dd_compose(qc, j, inter)
        row.update(pieces_second=len(sec), pieces_first=len(fir),
                   pieces_exact=exact, dd_final=ddf, dd_peak=ddp)
    if m == 10:
        c2_block = bb
    rows.append(row)
    exp.log(row)
exp.check("P1", p1, "Walsh = 3*2^(m-1)-2 for m=1..10")
exp.check("P2", p2, "interleaved 3m+4 plain / 3m+1 complement for m=1..10")
exp.check("P3", p3, "register-block >= 2^m for m=1..10")
exp.check("P4", p4, "affine pieces exact; op[2] <= 2m, op[1] = 2^m for m=1..8")
exp.fail_check("C2", c2_block is not None and c2_block >= 2 ** 10,
               f"block order m=10: {c2_block} nodes vs interleaved {3 * 10 + 4}")
exp.fail_check("C3", c3_violated, "quasi-reduced counts differ from 3m+4")
exp.fail_check("C4", c4_violated, "mutant adder Walsh differs from 3*2^(m-1)-2")
adder_dd_ok = all(r["dd_final"] == r["robdd_complement"]
                  for r in rows if r["fixture"] == "adder" and "dd_final" in r)

# ------------------------------------------------------------------ u_a ----
exp.section("P5-P8  u_a(ctrl, 2), x0, full dirty space")
FIXTURES = {n: [N for N in range(2 ** (n - 1) + 1, 2 ** n, 2)] for n in range(3, 7)}
FIXTURES[7] = [65, 73, 81, 89, 97, 105, 113, 121]
p6 = adder_dd_ok
p7 = True
p8 = True
per_n: dict[int, list[tuple[int, int]]] = {}
for n, Ns in FIXTURES.items():
    per_n[n] = []
    for N in Ns:
        t0 = time.time()
        me, qc, order = ua_fixture(N)
        assert me.n == n
        j = me.x[0]
        perm = permutation(qc)
        w = walsh_count(qc, j, perm)
        tab = output_bit_table(perm, j)
        del perm
        bp = robdd_size(tab, qc.n, order)
        bc = robdd_size(tab, qc.n, order, "complement")
        p8 &= bc <= bp <= 2 * (bc - 1) + 2
        row = dict(fixture="u_a", n=n, N=N, q=qc.n, gates=len(qc.logical),
                   walsh=w, robdd_plain=bp, robdd_complement=bc)
        if n <= 5:
            ddf, ddp = dd_compose(qc, j, order)
            row.update(dd_final=ddf, dd_peak=ddp)
            if n <= 4:
                p6 &= ddf == bc
            p7 &= ddp / ddf <= 1.35
        if n <= 4:
            pps = propagate_perm(qc, 1 << j, max_terms=8_000_000)
            row.update(pps_final=len(pps.final_terms), pps_peak=pps.n_max)
            p7 &= len(pps.final_terms) == w and pps.n_max / w >= 2
        row["seconds"] = round(time.time() - t0, 1)
        rows.append(row)
        per_n[n].append((w, bp))
        exp.log(row)
        sys.stdout.flush()

W = {n: statistics.median(w for w, _ in v) for n, v in per_n.items()}
B = {n: statistics.median(b for _, b in v) for n, v in per_n.items()}
R = {n: W[n] / B[n] for n in per_n}
for n in per_n:
    exp.log(f"n={n}: median W={W[n]} (log2 {math.log2(W[n]):.3f}) "
            f"median B={B[n]} (log2 {math.log2(B[n]):.3f}) R={R[n]:.2f}")
wsteps = {n: math.log2(W[n]) - math.log2(W[n - 1]) for n in range(4, 8)}
bsteps = {n: math.log2(B[n]) - math.log2(B[n - 1]) for n in range(4, 8)}
exp.check("P5a", all(2.8 <= s <= 3.2 for s in wsteps.values()),
          f"Walsh log2 steps { {k: round(v, 3) for k, v in wsteps.items()} }")
exp.check("P5b", all(1.0 <= s <= 2.7 for s in bsteps.values()) and bsteps[7] <= 2.2,
          f"ROBDD log2 steps { {k: round(v, 3) for k, v in bsteps.items()} }")
exp.check("P5c", all(R[n] > R[n - 1] for n in range(4, 8)) and R[7] >= 3 * R[5],
          f"R { {k: round(v, 2) for k, v in R.items()} }")
exp.check("P6", p6, "dd final == complement reducer (adder m<=8, u_a n=3,4)")
exp.check("P7", p7, "dd peak/final <= 1.35 (n<=5); PPS final == Walsh and peak/final >= 2 (n<=4)")
exp.check("P8", p8, "Bc <= B <= 2(Bc-1)+2 on every u_a fixture")

# -------------------------------------------------------------- control ----
exp.section("C1  random reversible circuit with u_a's gate multiset")
Rr = {}
for n in range(3, 7):
    me, qc, order = ua_fixture(2 ** (n - 1) + 1)
    rng = random.Random(20260915 + n)
    rq = Circuit(qc.n)
    for op in qc.logical:
        getattr(rq, op[0])(*rng.sample(range(qc.n), len(op) - 1))
    perm = permutation(rq)
    j = me.x[0]
    w = walsh_count(rq, j, perm)
    bp = robdd_size(output_bit_table(perm, j), rq.n, order)
    Rr[n] = w / bp
    row = dict(fixture="random_multiset", n=n, q=rq.n, walsh=w, robdd_plain=bp, ratio=Rr[n])
    rows.append(row)
    exp.log(row)
growth = Rr[6] / Rr[3]
exp.fail_check("C1", growth < 3, f"random ratio growth R_6/R_3 = {growth:.2f} "
               f"(u_a: {R[6] / R[3]:.2f})")

exp.finish(report_path=OUT / "report_v1.json", rows=rows,
           metadata=dict(task="Ta294b7d88a4c4015", attempt="A46391e1cf2954646",
                         gpu=accel.HAVE_GPU))
