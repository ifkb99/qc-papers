---
code: I
title: "IS THERE AN INTERMEDIATE 2-ADIC LAW? (TODO step 6). No. Answer is negative."
outcome: negative
claims: [C15]
todo: [6]
---
# I — IS THERE AN INTERMEDIATE 2-ADIC LAW? (TODO step 6). No. Answer is negative.

N=323 (r = 144 = 16·9) had shown density 0.981 at t=16 rising to 1.000 at t=24,
which looked like it might be a quantitative law in α rather than C15's binary
split. It is not.

### Two methodology bugs of my own, both instructive

1. **Degenerate tables.** Testing `g(e) = h[e mod r]` with random h: for r=3
   there are only 8 possible tables and 2 are constant, so small-r rows were
   contaminated and appeared *sparse*. Fixed by rejecting constants.
2. **Varied two things at once.** I drew a *new* random table for each t, so
   t-dependence was confounded with table variance — a direct violation of the
   "vary exactly one parameter" rule in `METHOD.md`. Fixed by fixing h and
   sweeping t.

Both produced confident-looking numbers before being caught. Worth remembering
that the protocol's rules catch *my* errors, not just other people's.

### What is solid

**β = 1: sparsity exactly constant in t.** Consistent with the C15 proof.

```
  r= 4 alpha=2 beta=1: sparsity = 4 4 4 4 4 4 4 4 4 4 4 4   (t = 11..22)
  r= 8 alpha=3 beta=1: sparsity = 4 4 4 4 4 4 4 4 4 4 4 4
  r=16 alpha=4 beta=1: sparsity = 16 16 16 16 16 16 ...
```

**β > 1: sparsity is Θ(2^t)** — density bounded well away from 0 at every t.

### What is NOT a law

Density within β>1 does **not** organise by α. With a fixed table and
consecutive t it shows a strong period-`ord₂(β)` oscillation *plus* a slow
upward drift, and neither component is a function of α:

```
  r=12 (beta=3, ord2=2): 1.000 0.625 1.000 0.625 ...   clean period 2
  r= 5 (beta=5, ord2=4): t≡1 mod 4 gives 0.791, 0.815, 0.832 ...  DRIFTS
  r= 7 (beta=7, ord2=3): t≡2 mod 3 gives 0.834, 0.831, 0.850, 0.875 ... DRIFTS
```

So the periodicity hypothesis (H1) is **refuted**: residue classes are not
constant, they drift. Density is periodic-plus-drifting, and function-dependent.

**Conclusion: there is no clean intermediate law. The robust structure is the
binary β=1 / β>1 dichotomy that C15 already asserts.** The N=323 observation is
best explained as sampling aliasing — t = 16, 20, 24 hit residues 4, 2, 0 mod
ord₂(9)=6, so three different points of the oscillation were read as a trend.

This is a negative result that *protects* the paper: it rules out a complication
rather than adding one, and it means C15 should be stated as the binary split
without hedging about intermediate regimes.

### Scope caveat CLOSED — rechecked with the real modexp bit function

`experiment_c15_realfn.py` reruns the same fixed-function, consecutive-t design
using the genuine table `h[c] = bit_j(a^c mod N)` instead of a random one.

**P1 holds.** β=1 gives constant sparsity across t=11..22 (4 for r=4; 8 for
r=16, so ≤ 2^α without being tight).

**Step 6's conclusion is confirmed: α is not the controlling parameter.**

```
  r= 3 alpha=0 beta= 3:  1.000 0.500 1.000 0.500 ...   exact period 2
  r= 6 alpha=1 beta= 3:  1.000 0.500 1.000 0.500 ...   identical to alpha=0
  r=12 alpha=2 beta= 3:  all 1.000
  r=10 alpha=1 beta= 5:  0.770 1.000 0.791 1.000 ...   drift 0.054
  r=18 alpha=1 beta= 9:  ~0.93-1.000                    drift 0.053
  r=22 alpha=1 beta=11:  ~0.92-1.000                    drift 0.015
```

α=1 produces wildly different behaviour across r = 6, 10, 18, 22, while α=0
(r=3) matches α=1 (r=6) exactly. **α does not organise the data.** The period is
set by ord₂(β), the values by β. Same qualitative structure as the random-table
study, so that study's conclusion transfers.

### NEW: the real modexp table is measurably NON-generic

Same r, real table vs random table:

```
  r= 6 (beta=3):  real  1.000 0.500 1.000 0.500 ...
                  random 1.000 1.000 1.000 1.000 ...
  r=10 (beta=5):  real  0.770 1.000 0.791 1.000 ...
                  random 0.500 0.500 0.394 0.197 ...
```

The real function is *sparser* than random at r=6 (exactly 0.500 on even t — a
clean factor of two, suggesting a dead variable or parity constraint) and
*denser* at r=10. So modexp bit functions carry structure beyond "depends on
e mod r". Not pursued further; noted as an open observation, and it means
random-table results should not be used as a proxy for real densities — only
for the qualitative α question they were built to answer.

**→ Now explained: see §AF.** The exact-0.500 rows are all-ones linear/affine
structures with a checkable three-ingredient criterion.
