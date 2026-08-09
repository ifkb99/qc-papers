---
id: 7
state: done
title: Does C15 survive truncation?
outcome: DONE, mostly yes
claims: [C24, C28]
---
# Does C15 survive truncation?

See `NOTES.md` §T; files `experiment_c15_trunc.py`, `experiment_c15_trunc2.py`.
Predictions were written before measuring.

- **P1 confirmed exactly.** Terminal thresholding preserves the constancy at
  *every* δ (counts identical across n_exp=3,4,5 at all seven δ values, both
  β=1 moduli); the β=3 control diverges. Follows from C24/M2 — identical
  magnitude multisets, so any magnitude threshold keeps identical counts.
- **P2, the practically relevant half: peak cost survives.** Under incremental
  truncation (what PPS actually does) `N_max` is SAME at every δ tested
  (24369 / 24369 / 1028 / 34). That is the quantity bounding memory.
- **C28, the caveat.** `N_final` drifts at aggressive δ (8 → 16 → 32) and ⟨O⟩
  collapses to 0 at δ=1e-1 for the wider circuits while n_exp=3 stays exact.
  The wider circuit is *more fragile at the same δ* despite an identical exact
  spectrum, because it has more gates and so more incremental truncation events.

Net: "period-finding precision is free" holds for memory at all δ, and for
accuracy up to moderate δ. Must not be stated unqualified.
