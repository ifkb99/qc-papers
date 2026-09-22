"""Does a virtual reset of clean scratch (C105's clear rule, TODO 66) give the exact
clean-input spectrum of Toffoli modexp while cutting PPS support?

Context. C45/C104 remove the exponent register from PPS memory; what remains is
exponential in the work width q_w = 3n + 4. C50 and C106 Theorem C put the excess on
dirty-scratch inputs. C105 proves a reset pulls back as Z^z -> Z^(z & ~mask). In
`ToffoliModExp`, t, c0, anc are |0> at every cc_add_mod boundary and b is |0> at
every u_a block boundary when x < N (the inverse cmult returns b = x_orig to 0 only
for x_orig < N). Inserting clears there is a "virtual reset". Derivation: TODO 66.

Objects. f = bit x0 of the circuit on all 2^q inputs. S = {b = t = c0 = anc = 0},
coordinates (x, exp). R0: clear at the input boundary only. R1: also at every block
boundary (mask b|t|c0|anc). R2: R1 plus t|c0|anc at every cc_add_mod boundary.
g = f with scratch zeroed after each block (classical replay, the independent route).

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  R0's final spectrum equals the normalized WHT of f|S exactly (subcube
      restriction identity).
  P2  R1's final spectrum equals the WHT of g (replay with block resets), exactly.
  P3  The functions read back from R1's and R2's spectra equal f on every input of
      S whose x at block-0 entry (x_in XOR 1, after build()'s X) is < N. Off that
      set R2 need not equal R1 (anc can be dirty at cc_add_mod boundaries).
  P4  Disagreements between R1's function and f|S occur only at x >= N with at
      least two exponent bits set (the first set bit leaves b dirty; only a later
      controlled block moves it into x).
  P5  With C45 contraction of the exponent qubits, every R1 block-boundary clear
      leaves at most 2^n terms.
  P6  The Shor-input expectation (x = 1, exponent |+>, scratch |0>) is unchanged by
      R1 and R2.
  EXPLORATION (no commitment): peak retained terms n_max for baseline, C45, R1+C45
      and R2+C45, and final support sizes. The within-block peak is not derived.

  C1  must fail: clearing t|c0 in the middle of the first `b -= N` ladder of the
      last block (t holds carries there) must disagree with f on some x < N input.
  C2  must fail: R1 must disagree with f|S somewhere at x >= N, popcount(e) >= 2
      (vacuous at t = 1, which is why it is resolved at t >= 2 only).

Run:  uv run python -m experiments.experiment_scratch_reset     (from research/)
Report: out/scratch_reset/report.json
"""
from __future__ import annotations

import numpy as np

import walsh
from circuits import Circuit
from perm_pps import propagate_perm
from toffoli_arith import ToffoliModExp
from lab import Experiment

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P1", "R0 final spectrum == WHT(f|S)")
exp.predict("P2", "R1 final spectrum == WHT(g), g = replay with block resets")
exp.predict("P3", "R1 and R2 functions == f on S with x < N (x at block-0 entry)")
exp.predict("P4", "R1 != f|S only at x >= N with popcount(e) >= 2")
exp.predict("P5", "C45 + R1: every block-boundary clear leaves <= 2^n terms")
exp.predict("P6", "Shor-input expectation unchanged by R1 and R2")
exp.must_fail("C1", "mid-ladder clear of t|c0 disagrees with f on some x < N")
exp.must_fail("C2", "R1 disagrees with f|S at some x >= N, popcount(e) >= 2 (t >= 2)")

CASES = [(7, 2, t) for t in (1, 2, 3)] + [(5, 2, t) for t in (1, 2, 3)]


def mask_of(qs) -> int:
    return sum(1 << q for q in qs)


def build(me: ToffoliModExp):
    """Rebuild me.build() while recording block and cc_add_mod boundaries."""
    qc = Circuit(me.n_qubits)
    qc.x(me.x[0])
    blocks, adds, mid = [len(qc.logical)], [], None
    for i, eq in enumerate(me.exp):
        c = pow(me.a, 1 << i, me.N)
        for mult, sign in ((c, +1), (pow(c, -1, me.N), -1)):
            pieces = []
            for j, xq in enumerate(me.x):
                piece = Circuit(me.n_qubits)
                me.cc_add_mod(piece, eq, xq, (mult * (1 << j)) % me.N)
                pieces.append(piece)
            if sign < 0:
                pieces = [p.inverse() for p in reversed(pieces)]
            for p in pieces:
                if sign > 0 and i == len(me.exp) - 1 and mid is None:
                    # first cc_add_mod of the last block: add_c, then b -= N
                    add_len = len(me._add_const(mult % me.N, (eq, me.x[0])).logical)
                    sub = me._add_const(me.N).inverse()
                    mid = len(qc.logical) + add_len + len(sub.logical) // 2
                qc.extend(p)
                adds.append(len(qc.logical))
            if sign > 0:
                for xq, bq in zip(me.x, me.b[: me.n]):
                    qc.cswap(eq, xq, bq)
                adds.append(len(qc.logical))
        blocks.append(len(qc.logical))
    assert qc.logical == me.build().logical, "rebuilt circuit differs from build()"
    return qc, blocks, sorted(set(adds) - set(blocks)), mid


def compact(me):
    """Coordinates of S: x bits then exponent bits."""
    return list(me.x) + list(me.exp)


def spectrum_on(coords, terms, n_total):
    """Dense vector over compact coordinates; raise if a key leaves them."""
    allowed = mask_of(coords)
    v = np.zeros(1 << len(coords))
    for z, c in terms.items():
        assert not z & ~allowed, f"key {z:#x} outside S coordinates"
        k = sum(((z >> q) & 1) << i for i, q in enumerate(coords))
        v[k] = c
    return v


def function_from(v):
    """0/1 truth table from a normalized +-1 spectrum (inverse WHT)."""
    chi = walsh.wht(v)
    assert np.allclose(np.abs(chi), 1.0, atol=1e-9), "spectrum is not a sign function"
    return (chi < 0).astype(np.int8)


def inputs_of(coords):
    """Full-space basis labels of the S inputs, in compact order."""
    k = np.arange(1 << len(coords), dtype=np.int64)
    y = np.zeros_like(k)
    for i, q in enumerate(coords):
        y |= ((k >> i) & 1) << q
    return k, y


def replay_with_resets(qc, cuts, mask, y):
    y = y.copy()
    prev = 0
    for k in cuts + [len(qc.logical)]:
        seg = Circuit(qc.n)
        seg.logical = qc.logical[prev:k]
        y = walsh.classical_images(seg, y)
        if k < len(qc.logical):
            y &= ~mask
        prev = k
    return y


rows = []
for N, a, t in CASES:
    me = ToffoliModExp(N=N, a=a, n_exp=t)
    qc, blocks, adds, mid = build(me)
    q, n, obs = me.n_qubits, me.n, 1 << me.x[0]
    scratch_b = mask_of(me.b + me.t + [me.c0, me.anc])
    scratch_a = mask_of(me.t + [me.c0, me.anc])
    coords = compact(me)
    kc, yS = inputs_of(coords)
    # x entering block 0 is x_in ^ 1: build() applies X(x0) first.
    xval = (kc & ((1 << n) - 1)) ^ 1
    evals = kc >> n
    pop = np.array([bin(int(e)).count("1") for e in evals])
    tag = f"N={N} a={a} t={t} q={q}"
    exp.section(tag)

    f_S = (walsh.classical_images(qc, yS) >> me.x[0]) & 1
    ref_fS = walsh.wht(np.where(f_S == 1, -1.0, 1.0)) / f_S.size
    g = (replay_with_resets(qc, blocks[1:-1], scratch_b, yS) >> me.x[0]) & 1
    ref_g = walsh.wht(np.where(g == 1, -1.0, 1.0)) / g.size

    r1_resets = {k: scratch_b for k in blocks[:-1]} | {0: scratch_b}
    r2_resets = r1_resets | {k: scratch_a for k in adds}
    base = propagate_perm(qc, obs)
    r0 = propagate_perm(qc, obs, reset_before={0: scratch_b})
    r1 = propagate_perm(qc, obs, reset_before=r1_resets)
    r2 = propagate_perm(qc, obs, reset_before=r2_resets)
    traced = propagate_perm(qc, obs, trace_plus=me.exp)
    r1t = propagate_perm(qc, obs, trace_plus=me.exp, reset_before=r1_resets)
    r2t = propagate_perm(qc, obs, trace_plus=me.exp, reset_before=r2_resets)
    bad = propagate_perm(qc, obs, reset_before={mid: mask_of(me.t + [me.c0]), 0: scratch_b})

    v0, v1, v2 = (spectrum_on(coords, r.final_terms, q) for r in (r0, r1, r2))
    e0 = float(np.max(np.abs(v0 - ref_fS)))
    exp.check("P1", e0 < 1e-12, f"{tag}: max |R0 - WHT(f|S)| = {e0:.1e}")
    e1 = float(np.max(np.abs(v1 - ref_g)))
    exp.check("P2", e1 < 1e-12, f"{tag}: max |R1 - WHT(g)| = {e1:.1e}")
    h1, h2 = function_from(v1), function_from(v2)
    low = xval < N
    exp.check("P3", bool(np.all(h1[low] == f_S[low]) and np.all(h2[low] == f_S[low])),
              f"{tag}: R1 and R2 == f on {int(low.sum())} x<N inputs;"
              f" R2 != R1 on {int((h2 != h1).sum())} (x>=N, not predicted)")
    diff = h1 != f_S
    outside = diff & ~((xval >= N) & (pop >= 2))
    exp.check("P4", not outside.any(),
              f"{tag}: {int(diff.sum())} disagreements, {int(outside.sum())} outside x>=N & |e|>=2")
    after = [ev[3] for ev in r1t.reset_events if ev[1] == scratch_b]
    exp.check("P5", max(after) <= (1 << n),
              f"{tag}: block-boundary clears leave {sorted(set(after))} terms (bound {1 << n})")
    ok6 = abs(r1t.expectation - traced.expectation) < 1e-12 and \
        abs(r2t.expectation - traced.expectation) < 1e-12
    exp.check("P6", ok6, f"{tag}: <O> {traced.expectation:+.12f} / R1 {r1t.expectation:+.12f}"
                         f" / R2 {r2t.expectation:+.12f}")

    hb = function_from(spectrum_on(coords, bad.final_terms, q))
    exp.fail_check("C1", bool(np.any(hb[low] != f_S[low])),
                   f"{tag}: mid-ladder clear disagrees on {int((hb[low] != f_S[low]).sum())}"
                   f" of {int(low.sum())} x<N inputs")
    if t >= 2:
        exp.fail_check("C2", bool(diff.any()),
                       f"{tag}: R1 differs from f|S on {int(diff.sum())} inputs")

    row = dict(N=N, a=a, t=t, q=q, n=n,
               nmax_base=base.n_max, nmax_c45=traced.n_max,
               nmax_r1=r1.n_max, nmax_r1_c45=r1t.n_max, nmax_r2_c45=r2t.n_max,
               final_base=len(base.final_terms), final_r1=len(r1.final_terms),
               final_c45=len(traced.final_terms), final_r1_c45=len(r1t.final_terms),
               disagreements_r1=int(diff.sum()))
    rows.append(row)
    exp.log("n_max  base {nmax_base}  C45 {nmax_c45}  R1 {nmax_r1}  R1+C45 {nmax_r1_c45}"
            "  R2+C45 {nmax_r2_c45}".format(**row))
    exp.log("final  base {final_base}  R1 {final_r1}  C45 {final_c45}  R1+C45 {final_r1_c45}"
            .format(**row))

exp.finish(report_path="out/scratch_reset/report.json", rows=rows)
