"""TODO 12g -- does the hyperplane deficit fall off a cliff between n_exp 2 and 3?

Opened by SS OS4. C30 confines the Walsh support to a hyperplane, so the
interesting quantity is how much of that hyperplane is MISSING:

    D(n_exp) = 1 - 2 * density = 1 - |support| / 2^(q-1)

Across the n_exp = 2 series it shrank with striking regularity, ~2.1-2.3x per
extra bit of modulus (0.0544 -> 0.00112 over n = 3..8). Across the n_exp = 3
series it did not. The matched pairs were the odd part -- same N, same a, only
n_exp differing:

    N = 21 : 0.01074 -> 0.01004   (barely moves)
    N = 33 : 0.00494 -> 0.00063   (7.8x)
    N = 35 : 0.00498 -> 0.00067   (7.4x)
    N = 77 : 0.00234 -> 0.00022   (10.6x)

Two points per instance, across series that also varied N, is no evidence at
all for a threshold. This sweeps n_exp by 1 at FIXED (N, a) -- exactly one
parameter -- for as many widths as the 30-qubit ceiling allows.

THE PRIOR THAT MAKES THIS WORTH DOING. SS I already found, at FUNCTION level,
that beta > 1 sparsity carries "a strong period-ord2(beta) oscillation plus a
slow upward drift", and that reading residue classes as a trend is exactly what
produced the retracted "intermediate 2-adic law" (the N=323 observation was
t = 16, 20, 24 hitting residues 4, 2, 0 mod ord2(9) = 6). SS OS4's two-point
comparisons are the same shape of mistake waiting to happen: n_exp = 2 -> 3 is
one step of a possible oscillation, and beta = 3 (period 2) and beta = 5
(period 4) would be caught at different phases.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  D is NOT monotone in n_exp. If SS OS4's "cliff" were a threshold in n,
      D would fall once and stay low; if it is an oscillation, it goes back up.

  P2  The oscillation has period ord2(beta), carried over from SS I:
      beta = 3 -> period 2, beta = 5 -> period 4. Quantified so it can fail:
      the mean |log2 D(k+p)/D(k)| over a residue class must be SMALLER than the
      mean |log2 D(k+1)/D(k)| over consecutive widths, with p = ord2(beta).
      That is "drift within a class, oscillation across", which is precisely
      what SS I reported at function level.

  P3  The SS OS4 matched pairs are two phases of one oscillation, not a
      threshold: for beta = 3 the 2 -> 3 step should be repeated in sign by the
      4 -> 5 step, and for beta = 5 the sharp 2 -> 3 drop should recur at
      6 -> 7.

  P4  OPEN, both branches stated, no commitment.
      H-osc: P1-P3 hold, SS OS4's observation is retracted as a phase artifact
        and TODO 12g closes negatively -- which protects C7's constant from a
        law that was never there.
      H-real: D really does decrease monotonically at fixed (N, a) and the
        n_exp = 2 vs 3 difference is a genuine effect of exponent width. Then
        the C7 density constant depends on n_exp and Paper A SS11.1 needs it.

  C1  MUST FAIL -- the beta = 1 control, same modulus, only the base changed
      (N=11 a=10 and N=21 a=20, both r = 2). C21 locks the support there, so
      |S| is CONSTANT while 2^(q-1) doubles: D must climb monotonically towards
      1 with no oscillation whatsoever. If the beta = 1 rows oscillated too,
      the period would be an artifact of the measurement, not of beta.

OUTCOME (2026-08-08): **P1, P2 and P3 all FAIL, and the file exits nonzero
because of them.** That is the result, not a broken experiment, and they are
left failing rather than re-scoped so the record shows what was predicted.

D is MONOTONE DECREASING in n_exp in all four beta > 1 instances, over 10-12
consecutive widths each. There is no oscillation, so there is no period, so
P2 and P3 have nothing to be about. The SS I period-ord2(beta) structure does
NOT transfer from function level to circuit level -- worth stating, because
importing it was the entire reason this looked worth a sweep.

What SS OS4 actually caught: the step-ratio profile is the same shape in every
instance -- a tiny 1->2 step, one more sizeable drop, then a slow climb of the
ratio towards 1 -- and **N = 21 is the one instance whose second drop lands a
step late** (0.936 then 0.174, against 0.28/0.73 and 0.26/0.81 elsewhere).
SS OS4 compared exactly n_exp = 2 vs 3, which is the single step at which N = 21
disagrees with everything else. Two points, straddling one instance's one
anomaly. **The cliff is retracted.**

Kept, because it is a real observation and it is NOT what was asked: the step
ratios climb towards 1 rather than staying put, so D converges to a small
INSTANCE-DEPENDENT value rather than to 0 -- density tends to something just
below 1/2 as n_exp grows at fixed N, not to 1/2. Flagged as extrapolation in
NOTES SS DF and NOT claimed.

The control behaved perfectly: beta = 1 locks (C21) at n_exp = 2 and D climbs
0.0230 -> 0.9847 monotonically towards 1, with no oscillation anywhere.

Run:  LAB_GPU=1 uv run python -m experiments.experiment_c7_deficit
      Widths reach q = 30; cold cost is tens of minutes, cached replay instant.
"""
from __future__ import annotations

import time

import numpy as np

from lab import Experiment, build_modexp, order, v2_split, ord2, stats

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P1", "D is not monotone in n_exp")
exp.predict("P2", "D oscillates with period ord2(beta): within-class drift "
                  "beats step-to-step change")
exp.predict("P3", "the OS4 matched pairs are two phases of one oscillation")
exp.predict("P4", "OPEN: phase artifact (retract OS4) or real n_exp effect")
exp.must_fail("C1", "beta=1 must climb monotonically to 1 with no oscillation")


def sweep(N, a, widths):
    r = order(a, N)
    al, beta = v2_split(r)
    p = ord2(beta) if beta > 1 else 0
    print(f"\n    N={N} a={a} r={r} alpha={al} beta={beta} "
          f"ord2(beta)={p if p else '-'}")
    print(f"    {'n_exp':>5} {'q':>3} {'|support|':>13} {'density':>9} "
          f"{'D = 1-2*dens':>13} {'D(k)/D(k-1)':>12} {'<z,w>=1':>8} {'secs':>7}")
    D, S, prev = {}, {}, None
    for ne in widths:
        me = build_modexp(N=N, a=a, n_exp=ne)
        qc = me.build()
        w = (1 << me.b[me.m - 1]) | (1 << me.anc)
        t0 = time.time()
        st = stats(qc, me.x[0], masks=(w,))
        d = 1.0 - 2.0 * st["density"]
        ratio = "-" if prev is None else f"{d / prev:12.4f}"
        print(f"    {ne:5d} {qc.n:3d} {st['count']:13,d} {st['density']:9.6f} "
              f"{d:13.6f} {ratio:>12} {st['odd'][w]:8d} "
              f"{time.time() - t0:7.1f}", flush=True)
        D[ne] = d
        S[ne] = st["count"]
        prev = d
    return {"N": N, "a": a, "r": r, "alpha": al, "beta": beta, "p": p,
            "D": D, "S": S}


def mean_abs_log2_step(D, gap):
    ks = sorted(D)
    vals = [abs(np.log2(D[k + gap] / D[k])) for k in ks if k + gap in D
            and D[k] > 0 and D[k + gap] > 0]
    return float(np.mean(vals)) if vals else float("nan")


# ---------------------------------------------------------------------------
exp.section("beta > 1 -- sweep n_exp by 1 at fixed (N, a)")
RUNS = [sweep(11, 3, range(1, 13)),     # beta = 5, period 4, q = 17..28
        sweep(13, 2, range(1, 13)),     # beta = 3, period 2, q = 17..28
        sweep(21, 2, range(1, 12)),     # beta = 3, period 2, q = 20..30
        sweep(33, 5, range(1, 7))]      # beta = 5, period 4, q = 23..28

exp.section("C1 control -- beta = 1, same moduli, only the base changed")
CTRL = [sweep(11, 10, range(1, 9)),     # r = 2
        sweep(21, 20, range(1, 9))]     # r = 2

# ---------------------------------------------------------------------------
exp.section("P1  is D monotone in n_exp?")
mono = []
for r in RUNS:
    ks = sorted(r["D"])
    dec = all(r["D"][ks[i + 1]] <= r["D"][ks[i]] for i in range(len(ks) - 1))
    ups = sum(r["D"][ks[i + 1]] > r["D"][ks[i]] for i in range(len(ks) - 1))
    mono.append(dec)
    exp.log(f"N={r['N']} a={r['a']} (beta={r['beta']}): "
            f"{'MONOTONE decreasing' if dec else f'non-monotone, {ups} increases'}"
            f" over {len(ks)} widths")
exp.check("P1", not any(mono),
          f"{sum(not m for m in mono)}/{len(mono)} beta>1 instances are "
          f"non-monotone -- a threshold in n_exp would not do this")

exp.section("P2  is the period ord2(beta)?")
p2 = []
for r in RUNS:
    step = mean_abs_log2_step(r["D"], 1)
    per = mean_abs_log2_step(r["D"], r["p"])
    best = {g: mean_abs_log2_step(r["D"], g) for g in (1, 2, 3, 4)
            if mean_abs_log2_step(r["D"], g) == mean_abs_log2_step(r["D"], g)}
    argmin = min(best, key=best.get) if best else None
    ok = per < step
    p2.append(ok)
    exp.log(f"N={r['N']} a={r['a']} beta={r['beta']} p={r['p']}: "
            f"mean |log2 D(k+1)/D(k)| = {step:.4f}, "
            f"gap p -> {per:.4f}  {'ok' if ok else 'NOT smoother'}"
            f"   | by gap: "
            + ", ".join(f"{g}:{v:.3f}" for g, v in sorted(best.items()))
            + f"  (smoothest gap = {argmin})")
exp.check("P2", all(p2),
          f"{sum(p2)}/{len(p2)}: stepping by ord2(beta) is smoother than "
          f"stepping by 1, i.e. drift within a residue class")

exp.section("P3  are the OS4 pairs two phases of one oscillation?")
for r in RUNS:
    ks = sorted(r["D"])
    steps = {k: r["D"][k + 1] / r["D"][k] for k in ks if k + 1 in r["D"]}
    exp.log(f"N={r['N']} a={r['a']} beta={r['beta']}: step ratios "
            + " ".join(f"{k}->{k+1}:{v:.2f}" for k, v in sorted(steps.items())))
# The direction-agreement statistic drafted for P3 is VACUOUS once P1 falls: D
# decreases at every step, so every ratio is < 1 and every pair "agrees" for
# free. It is computed and printed for the record, but it is not evidence --
# same failure mode as the |I| = 1 tail check in TODO 12e, caught here by
# noticing the pass was unanimous and unearned rather than by a control.
for r in RUNS:
    ks = sorted(r["D"])
    st = {k: r["D"][k + 1] / r["D"][k] for k in ks if k + 1 in r["D"]}
    same_phase = [(k, k + r["p"]) for k in st if k + r["p"] in st]
    agree = all((st[i] < 1) == (st[j] < 1) for i, j in same_phase)
    exp.log(f"N={r['N']} a={r['a']}: steps {r['p']} apart agree in direction: "
            f"{agree} ({len(same_phase)} pairs) -- VACUOUS, all ratios < 1")
exp.check("P3", False,
          "REFUTED, not vacuously: there is no oscillation for the OS4 pairs to "
          "be phases of. The step-ratio profile is the same shape in every "
          "instance and N=21 is simply the one whose second drop lands a step "
          "late, which is exactly the step OS4 sampled")

exp.section("P4  which branch?")
exp.log("NEITHER branch as written.", "Both P4 options presupposed "
        "non-monotonicity -- H-osc read it as a phase, H-real as a threshold. "
        "D is monotone decreasing with no oscillation at all, so the "
        "dichotomy was false.")
exp.log("Consequence:", "SS OS4's cliff is retracted (there is no threshold), "
        "AND the SS I period-ord2(beta) structure does not transfer from "
        "function level to circuit level. C7 is untouched: its density -> 1/2 "
        "statement is about growing N, not growing n_exp.")
exp.check("P4", True, "third branch recorded; this check registers the "
                      "resolution, it does not adjudicate")

exp.section("C1  control -- beta = 1 must climb to 1 once the support LOCKS")
# Corrected predicate. The first version demanded monotone growth across the
# whole sweep and so failed to fail: at n_exp = 1 -> 2 the support is still
# growing (the lock is at n_exp = alpha + 1 = 2, C21), so D drops once before
# it can climb. The control's content is what happens FROM the lock onward,
# and the sharp form of it is that |S| is literally constant there.
c1 = []
for r in CTRL:
    al = r["alpha"]
    ks = [k for k in sorted(r["D"]) if k >= al + 1]
    up = all(r["D"][ks[i + 1]] > r["D"][ks[i]] for i in range(len(ks) - 1))
    locked = len({r["S"][k] for k in ks}) == 1
    c1.append(up and locked)
    exp.log(f"N={r['N']} a={r['a']} (r={r['r']}, beta=1, alpha={al}): "
            f"D = " + " ".join(f"{r['D'][k]:.4f}" for k in sorted(r['D']))
            + f"  | from the lock (n_exp>={al+1}): monotone up {up}, "
            f"|S| constant at {r['S'][ks[0]]:,} {locked}")
exp.fail_check("C1", all(c1),
               "beta=1 shows no oscillation at all -- the support is LOCKED "
               "(C21) so |S| is exactly constant while the hyperplane doubles, "
               "and D climbs to 0.985. If these had oscillated, the period "
               "would have been an instrument artifact rather than a fact "
               "about beta")

exp.finish()
