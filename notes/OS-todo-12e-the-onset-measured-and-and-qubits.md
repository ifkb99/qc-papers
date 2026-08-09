---
code: OS
title: "TODO 12e: THE ONSET MEASURED AT α = 3 AND 4 (and C7 to 30 qubits)"
date: 2026-08-08
outcome: record
claims: [C2, C7, C21, C23, C24, C30, C37, C38, C43]
todo: [12e]
---
# OS — TODO 12e: THE ONSET MEASURED AT α = 3 AND 4 (and C7 to 30 qubits)

TODO 12e was not a research question, it was a *compute* question: several
sharp predictions had been made and left unmeasured because the runs were too
slow. With the CUDA backend they are minutes. Files:
`experiments/experiment_c21_onset.py`, `experiment_c7_scale.py`,
`experiment_windowed_scale.py`. All three declare predictions before measuring
and carry must-fail controls.

### OS1 — C21's onset is now MEASURED at α = 3 and α = 4, not merely predicted

This was the single most valuable rerun: C21 ("the support locks at
n_exp = v₂(r) + 1") was exact at α = 1, 2 and, at α = 3 and 4, only
"still growing at the largest width we can reach, as predicted". Now:

```
  N=17 a=2  r=8  α=3      N=17 a=3  r=16 α=4      N=41 a=3  r=8  α=3
 n_exp  q   |support|    n_exp  q   |support|    n_exp  q   |support|
   1   20     255,104      1   20     255,356      1   23   2,070,878
   2   21   1,037,405      2   21   1,037,374      2   24   8,346,567
   3   22   2,093,137      3   22   2,093,202      3   25  16,766,478
   4   23   4,188,525 <    4   23   4,188,537      4   26  33,539,711 <
   5   24   4,188,525 =    5   24   8,379,626 <    5   27  33,539,711 =
                           6   25   8,379,626 =
```

**Locks exactly at n_exp = α + 1 in all three, with strict growth at every
step below it (7/7).** α = 4 is the deepest onset ever measured here, and the
two α = 3 rows are different moduli, so the rule is not a property of N = 17.

**Matched must-fail control, and it is the good kind: same modulus, same α,
same width, only β differs.** N = 41 has elements of order 8 (β = 1) *and* of
order 40 = 5·2³ (β = 5). The β = 5 row grows at every step —
2,072,174 / 8,346,761 / 16,766,484 / 33,539,777 / 67,086,624 — and never
locks. Note how close it is to its β = 1 twin at n_exp = 4 (33,539,777 vs
33,539,711, a gap of 66 in 3.4e7): **at any single width the two are
indistinguishable; only the growth separates them.** That is worth carrying
into the papers — the invariant is not visible in a cost measurement at one
size.

### OS2 — C24 strengthens from equal COUNTS to identical SETS (C43)

C23/C24 (iv) computes the surviving coefficient explicitly, and the expression
does not mention |I| at all — so the proof gives more than the constant-size
statement that was tested. Write each support element as (z_rest, tailflag)
with tailflag ∈ {0,1} saying whether z_I = 0 or 1_I. Then the *set* of such
pairs must be identical at consecutive widths, not merely the same size.

**Derived before measuring, then confirmed 3/3** (N=17 a=2 at n_exp 4 vs 5;
N=17 a=3 at 5 vs 6; N=41 a=3 at 4 vs 5) — bit-for-bit identical sorted key
arrays over supports of 4.2M, 8.4M and 33.5M elements. Logged as **C43**.

The two halves are individually constant as C24 says, and unequal to each
other (e.g. 2,094,285 vs 2,094,240) — which is expected: C24 predicts each
half is width-independent, not that they match.

### OS3 — METHOD: the must-fail control caught a vacuous test of my own

C2 **passed on the first run**, i.e. failed to fail. Diagnosis: the tail
confinement test "z_I ∈ {0, 1_I}" is vacuous at |I| = 1, because a single bit
*is* either all-zeros or all-ones. The β = 5 control was being asked a
question with only one possible answer — and so, silently, were three of P4's
six rows.

Nothing measured was wrong; the scope of what it could testify to was. Both
checks were re-scoped to |I| ≥ 2, which leaves P4 with 3 genuinely
non-vacuous confirmations (still 0 violations) and gives the control 33.5M
violations out of 67M, as required.

**Second time in the project a must-fail control has caught a vacuous
measurement** — the first was step 9 (§G), where synthetic blocks acted only
on b-qubits while the observable was `Z_x0`. Not the third: the project's
other two self-caught errors (a sweep drawing a new random table per t, and
degenerate random tables at small r) came from "vary exactly one parameter"
and from null-model hygiene, not from this rule. The rules are not
interchangeable and the tally should not be inflated.

### OS4 — C7's circuit series extended to 30 qubits

See `experiment_c7_scale.py`. Two series, each varying only N at fixed n_exp,
measured on the exact integer path.

**P0 first, because it gates everything else:** all six logged C7 rows
reproduce bit-for-bit on the exact (`!= 0`) support test — 15,493 / 15,539 /
127,936 / 1,037,322 / 8,347,241 / 8,346,759. So the thresholding artifact of
§GF does **not** touch the circuit-level numbers, which is what the arithmetic
says it should do: coefficients are integers/2^q, and at q ≤ 30 the smallest
nonzero magnitude is 2⁻³⁰ ≈ 9.3e-10, three orders above the 1e-12 threshold.
(It bites at function level, where n is much larger.)

```
  A: n_exp = 2, vary N              B: n_exp = 3, vary N (β>1 only)
   N   n   r    q   |support|  dens    N   n   r    q   |support|  dens
   5   3   4   15      15,493  .4728    7   3   6   16      30,712  .4686
   7   3   6   15      15,539  .4742   21   5   6   22   2,076,089  .4950
  15   4   4   18     127,936  .4880   33   6  10   25  16,766,640  .4997
  21   5   6   21   1,037,322  .4946   35   6  12   25  16,766,067  .4997
  33   6  10   24   8,347,241  .4975   77   7  30   28 134,188,859  .4999
  35   6  12   24   8,346,759  .4975
  77   7  30   27  66,952,166  .4988   slope  A 1.0054   B 1.0083
 143   8  20   30 536,271,623  .4994          pooled 1.0062  (was 1.008)
```

**q = 30, N = 143, a = 5: |support| = 536,271,623, density 0.49944.** Half a
billion Pauli terms, exactly counted, in 851 s on one A4500 — 64× the Hilbert
space C7 settled on. The slope stays at 1.006 bits/qubit over q = 15..30, and
the C30 linear structure w = b_msb ⊕ anc is still exactly a structure (0
support elements with ⟨z,w⟩ = 1) in 13/13 circuits.

Both controls behave: growing a β=1 instance in n_exp instead of N collapses
the density by exactly 4× per two qubits (0.4733 → 0.1183 → 0.0296, ratios
0.2500/0.2500), and an affine observable reads sparsity exactly 1 at every
width. So the ~½ reading is a property of these circuits, not of the
instrument or the register size.

**Open observation — RETRACTED 2026-08-08 by §DF, same day. Kept because the
way it was wrong is the point.** Write the *missing fraction* of the C30
hyperplane, 1 − 2·density. In series A it shrinks with almost perfect
regularity, ~2.1–2.3× per extra bit of modulus:
0.0544 / 0.0239 / 0.0107 / 0.00494 / 0.00234 / 0.00112 at n = 3..8. In series
B it does **not**: 0.0627 / 0.0100 / 0.00063 / 0.00022 at n = 3, 5, 6, 7.
Matched pairs (same N and a, only n_exp differing) appeared to split the same
way — N = 21 barely moving (0.01074 → 0.01004) while N = 33, 35, 77 each drop
8–11× — which read as a threshold in n_exp.

> **There is no threshold.** Sweeping n_exp by 1 at fixed (N, a) over 10–12
> widths (§DF) shows D decreasing **monotonically** in every instance, with the
> same step-ratio profile everywhere, and N = 21 the lone instance whose second
> drop lands one step late. The matched pairs above compared n_exp = 2 vs 3 —
> precisely the step at which N = 21 disagrees with everything else. **Two
> points straddling one instance's one anomaly.** The lesson is the one §I
> already taught and this failed to apply: a two-point comparison cannot
> distinguish a trend from a phase, and the fix is to sweep the parameter, not
> to collect more pairs. Also do not quote the "halves per qubit" reading; it is
> series A only.

### OS6 — C38's 2-periodicity holds to K = 5, i.e. three full periods

§WD stopped at K = 3, which is two periods — the bare minimum that can be
called a period at all, and exactly the shape that produced the retracted
"intermediate 2-adic law" (§I: three points of an oscillation read as a trend).
So this was worth two more points rather than none.

```
  K              0        1        2        3        4         5    q at K=5
  N=5 a=4   28,078   62,680   28,078   62,680   28,078    62,680       29
  N=7 a=6   13,563   62,570   13,563   62,570   13,563    62,570       29
  control   84,827 1,009,765 4,156,585 1.66e7 6.66e7  266,296,901      29
  (N=7 a=3, r=6, beta=3 — no tail window, so nothing is dead; ×4 per window)
```

Exact repeats, not near-misses. Two strengthenings beyond the counts, both
predicted first:

- **Sets, not sizes.** The K = 4 support restricted to the K = 0 coordinates
  equals the K = 0 support as a set, and likewise K = 5 against K = 1 — 4/4.
  If W² = id really makes the longer circuit the same permutation on a wider
  register, that is what has to happen, and it does.
- **C37's dead tail survives arbitrary tail length.** Live exponent bits are
  exactly {0…α−1} at K = 4 and 5 in both moduli; no tail bit ever enters the
  support, at any length.

### OS5 — infrastructure that came out of it (all gated in `test_accel.py`)

- **int32 permutation replay.** The replay, not the transform, is the binding
  memory constraint: `idx ^= ((idx >> c) & 1) << t` keeps the array plus two
  temporaries live, so int64 needs ~24 GiB at n = 30 and does not fit on a
  20 GiB card. Images are < 2ⁿ, so int32 is exact for n ≤ 30 and needs ~12
  GiB. **This is what makes q = 30 reachable at all**, and it is why the
  handoff's "n ≤ 30" figure — which was measured on the FWHT alone — was
  optimistic for circuit-level work until now.
- **`accel.pullback_stats`** returns count / density / per-mask GF(2) parity
  counts without moving the support to the host. At q = 30 the support is
  4 GiB as int64; a density sweep wants three scalars.
- **`lab.measure.support(..., exact=True)` and `lab.measure.stats`**, both
  cached under their own keys so exact and thresholded results can be
  compared rather than silently substituted.
