---
id: hk8
state: done
title: GPU acceleration
outcome: DONE 2026-08-08
claims: []
---
# GPU acceleration

**GPU acceleration — DONE 2026-08-08.** `accel.py` + `test_accel.py`
(suite 9). Both hot paths are memory-bound array passes, so a card with ~10×
the bandwidth wins by about that: measured **7.1× / 6.9× / 11.5×** on
permutation replay at q = 17/19/21 and **11.8× / 12.8× / 14.7×** on FWHT at
2^24/26/28 (RTX A4500). The speedup **grows with n**, which is the useful
direction. Opt-in via `LAB_GPU=1`; `walsh.py` stays the reference and the GPU
path is gated against it (exact equality on permutations, identical support
sets, plus must-fail controls for non-permutation input and oversized
registers). Capacity: n ≤ 30 on one 20 GiB card.
