---
id: 2
state: done
title: Weight-truncation as Fourier tail mass
outcome: DONE, negative result
claims: [C14]
---
# Weight-truncation as Fourier tail mass

**Outcome: the exact model does NOT extend to weight truncation, and weight
truncation should not be used on reversible arithmetic.** See `NOTES.md` §W.

- **W1** The naive identity is false. PPS truncates *incrementally*, so a
  discarded term never branches and the result is not the truncated final
  operator. The Fourier tail describes **terminal** truncation only. Both are
  non-monotonic in k; incremental can be 400× better (k=3) or far worse (k=8).
  This also retroactively sharpens the δ story (C14): non-monotonicity is not
  only cancellation, it is truncation changing what subsequently branches.
- **W2** Signed sums per weight level are large and alternating
  (−0.500, −0.516, +0.117, …, −0.516, +0.297), cancelling only over all levels.
  Any cutoff slices the cancellation. No usable k exists.
- **W3** By contrast δ is exactly right: the `|c|>0.1` set (4 of 3086 terms)
  sums to **exactly** ⟨O⟩ = −1.00000000 and the remaining 3082 sum to
  **exactly** 0. Clean split.
- **W4** The four dominant coefficients sit at weights **1, 2, 8, 9** — not all
  low-degree, which is precisely why weight truncation cannot substitute for δ.
  (Also corrects a misreading: the "292× low-weight enrichment" is just those two
  0.5-magnitude coefficients; `k(99% mass)=11`, so the bulk is high-weight.)

Files: `experiment_weight.py`, `experiment_weight2.py`; `perm_pps.py` gained a
`max_weight` argument.

**Subsumes step 8** (decode the dominant coefficients) — done as part of W4.
