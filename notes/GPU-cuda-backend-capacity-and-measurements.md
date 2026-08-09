---
code: GPU
title: "CUDA backend: what it buys, where it caps, and why the replay binds"
date: 2026-08-08
outcome: record
claims: []
todo: [12e, 12f]
---

# GPU — read this before running anything slow

The two hot paths (permutation replay, FWHT) are memory-bound array passes, and
a CUDA backend gives **7–15×, growing with n**. This is already built, gated and
wired in. It will save hours.

```bash
LAB_GPU=1 uv run python -m experiments.<name>     # that is the whole interface
```

- **`accel.py`** is the backend; **`test_accel.py`** (suite 9) gates it against
  the CPU reference and skips cleanly with no card.
- **Opt-in by design.** `walsh.py` stays the reference and is never silently
  substituted. `LAB_GPU=1` routes `lab.measure.support` through the GPU; nothing
  else changes behaviour.
- **Measured** (RTX A4500): permutation replay 7.1× / 6.9× / **11.5×** at
  q = 17 / 19 / 21 (29.3 s → 2.6 s at q = 21); FWHT 11.8× / 12.8× / **14.7×** at
  2²⁴ / 2²⁶ / 2²⁸ (18.1 s → 1.2 s at 2²⁸).
- **Capacity: n ≤ 30** on one 20 GiB card (`accel.MAX_QUBITS`); it raises rather
  than thrashing past that. Measured: float64 FWHT at n = 30 needs 12 GiB and
  fits; n = 31 needs 24 GiB and does not.
- **The binding constraint is the REPLAY, not the transform** (found 2026-08-08
  while doing TODO 12e). `idx ^= ((idx >> c) & 1) << t` keeps the array plus two
  temporaries live, so an int64 index array needs ~24 GiB at n = 30 and fails.
  Images are < 2ⁿ, so `accel._replay` builds it in **int32** for n ≤ 30 (~12
  GiB), which is what makes a q = 30 circuit-level run possible at all. Gated
  against the CPU int64 reference in `test_accel.py` [A]. A real q = 30 modexp
  now takes ~14 min end to end.
- **Does the second card scale it further? Only by +1 qubit, and there is a
  cheaper way.** Memory doubles per qubit, so 40 GiB buys exactly one more than
  20 GiB. A split would be genuinely easy — with the array halved by its top
  bit, every FWHT level except the last is local to a half, and only the final
  butterfly crosses cards — but +1 qubit is a poor return, so it is **not
  implemented** (logged as TODO 12f).
  **Use `wht_exact` / `pullback_support_exact` instead:** the FWHT of ±1 data is
  *exactly integer-valued* (verified bit-for-bit), so int32 gives the same
  answer in half the memory — the same +1 qubit, no complexity — **and it makes
  the support test exact (≠ 0) rather than a magnitude threshold.** That second
  property matters: it removes the thresholding artifact that otherwise makes
  measured densities drift below their true value as n grows (see `NOTES.md`
  §GF and `experiment_gf2law_scale.py`). Safe to n = 30, since intermediate
  magnitudes are bounded by 2ⁿ and 2³⁰ < 2³¹.
- The best use of the second card is **throughput**: two independent sweeps at
  once, one per device, via `CUDA_VISIBLE_DEVICES=0` / `=1`. Zero new code —
  and it is how 12e was actually run, two sweeps in parallel throughout. The
  measurement cache is content-addressed, so concurrent writers cannot collide
  and a later single-process run of the experiment replays everything for free.
  That is the pattern to reuse: **warm the cache in parallel, then run the
  experiment file once for the record.**
- **Install is machine-specific.** `pyproject.toml` pins `cupy-cuda13x` to match
  this machine's CUDA 13.3. On a CUDA 12 host swap to `cupy-cuda12x`; both were
  tested and perform identically. Without cupy everything still works on CPU.
- `accel.device_info()` prints what was detected; `accel.enabled()` tells you
  whether `LAB_GPU` actually took effect.

---
