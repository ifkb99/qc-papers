---
code: DF
title: "TODO 12g CLOSED NEGATIVELY: there is no cliff, and §OS4's CLIFF OBSERVATION is RETRACTED"
outcome: retracted
claims: [C1, C7, C21, C30]
todo: [12g]
---
# DF — TODO 12g CLOSED NEGATIVELY: there is no cliff, and §OS4's CLIFF OBSERVATION is RETRACTED

`experiments/experiment_c7_deficit.py`. **Three of its four predictions fail
and the file exits nonzero; that is the result.** No new claim; §OS4's
*hyperplane-deficit cliff* observation is withdrawn and TODO 12g closes.

> **Scope of this retraction — read before following a pointer to §OS4.** What
> dies here is one observation *inside* §OS4: the claimed cliff in
> `D = 1 − 2·density` between n_exp = 2 and 3. **§OS4's C7 result is
> untouched** — the 30-qubit circuit series, density 0.4728 → 0.49944, pooled
> slope 1.0062 bits/qubit, and the bit-for-bit exact-integer reproduction all
> stand, and `CLAIMS.md`'s C7 row correctly cites §OS4 for them. Do not read
> the header as retracting the section.

**What was asked.** §OS4 noticed that the missing fraction of the C30
hyperplane, `D = 1 − 2·density`, dropped 8–11× between n_exp = 2 and 3 for
N = 33, 35 and 77 but barely moved for N = 21 (0.01074 → 0.01004). Two points
per instance, in series that also varied N. This swept n_exp by 1 at **fixed
(N, a)** — one parameter — over 10 to 12 consecutive widths per instance.

```
  n_exp     1        2        3        4        5        6       ...     12
 N=11  .528687  .023872  .006561  .004818  .003883  .003302  ...  .002125   β=5
 N=13  .530884  .023491  .006027  .004881  .004114  .003696  ...  .002187   β=3
 N=33  .506047  .004931  .000630  .000468  .000359  .000284               β=5
 N=21  .512184  .010733  .010044  .001746  .001472  .001164  ...  .000613   β=3
 ctrl  .524780  .022964  .511482  .755741  .877871  .938935  ...  .984734   β=1
```

(N = 11 and 13 run to n_exp = 12 / q = 28, N = 21 to n_exp = 11 / q = 30 with
|S| = 536,542,026, N = 33 to n_exp = 6, control to n_exp = 8.)

**P1 (D is not monotone) — REFUTED, 0/4.** D decreases monotonically at every
one of ~40 consecutive steps. **P2 (period ord₂(β)) and P3 — vacuous once P1
falls: there is no oscillation, so there is no period.** The §I
period-ord₂(β) structure is a **function-level** phenomenon and does **not**
transfer to circuit level. Importing it was the whole reason this looked worth
a sweep, so that transfer failing is the useful part.

**What §OS4 actually caught.** The step-ratio profile is the same shape
everywhere — a tiny 1→2 step, one further sizeable drop, then a slow climb of
the ratio towards 1 — and **N = 21 is the single instance whose second drop
lands one step late**:

```
  step        1->2   2->3   3->4   4->5   5->6
  N=11        .045   .275   .734   .806   .850
  N=13        .044   .257   .810   .843   .898
  N=33        .010   .128   .743   .767   .791
  N=21        .021   .936   .174   .843   .791     <- second drop delayed by 1
```

§OS4 compared exactly n_exp = 2 vs 3 — the one step at which N = 21 disagrees
with everything else. **The cliff was two points straddling one instance's one
anomaly. Retracted.** Why N = 21 is late is not explained and is not pursued:
it is not α (N = 33 shares α = 1 and behaves normally) and not β (N = 13
shares β = 3 and behaves normally).

**The control is the clean part.** β = 1 at the same moduli (N = 11 a = 10 and
N = 21 a = 20, both r = 2) locks at n_exp = 2 exactly as C21 says, so |S| is
frozen — at 128,062 and 1,037,435 respectively, *bit-identical across six
widths* — while the hyperplane doubles and D climbs monotonically to 0.985. No
oscillation anywhere, which is what makes "no oscillation" in the β > 1 rows a
statement about the circuits rather than about the instrument.

**Two of my own predicates were wrong, in opposite directions, and both are
worth recording.**

- **C1 failed to fail on the first run.** I demanded monotone growth across the
  whole sweep, but the lock is at n_exp = α + 1 = 2, so D legitimately *drops*
  once (0.5248 → 0.0230) before it can climb. The control's content is what
  happens *from the lock onward*; the predicate included a pre-lock step.
  Sharpened while fixing it: the strong form is that |S| is **exactly
  constant** past the lock, not merely that D rises — which is a cleaner
  statement of C21 than the one being tested.
- **P3 passed vacuously.** Its statistic was "do steps p apart agree in
  direction" — but once D is monotone every ratio is < 1, so every pair agrees
  for free, and it returned 4/4 unanimous. Caught by noticing the pass was
  unanimous and unearned. P3 is now recorded as **refuted**, which is its
  actual status: there is no oscillation for the §OS4 pairs to be phases of.
  Same failure mode as the |I| = 1 tail check in §OS3, one week's worth of
  lesson apart — **a statistic that cannot come out any other way is not
  evidence, and monotone data makes sign tests vacuous.**

**One observation kept, and explicitly NOT claimed.** The step ratios climb
towards 1 rather than settling, so D looks to converge to a small
*instance-dependent* constant rather than to 0 — i.e. as n_exp grows at fixed
N, density tends to something just below ½, not to ½. Fitting the successive
differences of the N = 11 row (which decay by a near-constant ~0.71) puts
D∞ ≈ 0.0020, density ≈ 0.4990. **This is extrapolation from a converging
sequence, which is exactly the move §I punished, so it is logged and not
believed.** It does not threaten C7, whose density → ½ statement is about
growing N and is measured, not extrapolated. If anyone wants it: the honest
route is a function-level calculation of the n_exp → ∞ limit, not more widths.
