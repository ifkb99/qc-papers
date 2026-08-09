---
id: 12f
state: open
title: Two-GPU split FWHT
outcome: only if a single run ever needs n = 31+
claims: []
---
# Two-GPU split FWHT

**Deliberately not done.** Memory doubles per qubit, so the second A4500 buys
exactly **one** more qubit (n = 30 → 31). Measured: float64 FWHT needs 12 GiB
at n = 30 (fits one card) and 24 GiB at n = 31 (does not).

The split itself is easy if it is ever wanted. Halve the array by its **top
bit**; then for every level h < size/2 the butterfly pairs stay inside one half,
so both cards run independently with no communication, and only the **final**
level pairs element i of half 0 with element i of half 1. That is one
peer-to-peer exchange of half the array (~8 GiB each way at n = 31), a couple of
seconds over PCIe.

**Do the cheap thing first.** `accel.wht_exact` already halves memory by using
int32 — legitimate because the FWHT of ±1 data is exactly integer-valued — which
buys the same +1 qubit with no split, and additionally makes the support test
exact rather than thresholded. Only build the split if something genuinely needs
n ≥ 32, and note int32 overflows past n = 30 (bound 2ⁿ vs 2³¹), so a split at
n = 32 needs int64 and therefore both cards anyway.

Meanwhile the second card is not idle-by-design: **run independent sweeps
concurrently**, `CUDA_VISIBLE_DEVICES=0` and `=1`. That is what TODO 12e wants
and it needs no new code.
