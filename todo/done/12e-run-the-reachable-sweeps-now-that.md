---
id: 12e
state: done
title: Re-run the reachable sweeps now that the GPU makes them cheap
claims: [C7, C21, C24, C30, C38, C43]
---
# Re-run the reachable sweeps now that the GPU makes them cheap

**DONE 2026-08-08. All three sub-items closed; every prediction confirmed.**
See `NOTES.md` §OS; files `experiments/experiment_c21_onset.py`,
`experiment_c7_scale.py`, `experiment_windowed_scale.py`. C21, C24, C7, C30
and C38 all gained evidence; C43 is new.

- **C21's onset at α = 3 and 4 — MEASURED (§OS1).** Locks at exactly
  n_exp = α+1 in all three instances: N=17 a=2 (r=8) and N=41 a=3 (r=8) at
  n_exp=4, N=17 a=3 (r=16) at n_exp=5, with strict growth at all 7 steps
  below. Two moduli for α=3, so it is not an N=17 artifact. Matched control
  (N=41 a=6, r=40=5·2³ — same modulus, same α, same width, only β differs)
  never locks. **Worth carrying into Paper B: at n_exp=4 the β=1 and β=5 rows
  differ by 66 in 3.4e7, so no single-width cost measurement can see the
  invariant — only the growth separates them.**
- **C24 strengthened to set level (C43, §OS2).** The proof's step (iv) gives
  the surviving coefficient explicitly with no |I| in it, so the supports must
  coincide as *sets* under the (z_rest, tailflag) encoding, not merely in
  size. Derived before measuring, confirmed 3/3 bit-for-bit over supports of
  4.2M / 8.4M / 33.5M.
- **C7 extended to 30 qubits (§OS4).** Two series, slope 1.0054 / 1.0083,
  pooled 1.0062 bits/qubit over q = 15..30; q=30 gives |support| =
  536,271,623 at density 0.49944 in 851 s. All six original rows reproduce
  bit-for-bit on the exact integer path. `PAPER_A.md` §11.1's concession that
  the circuit series is modest can now be softened.
- **C38's 2-periodicity extended to K = 5 (§OS6)** — four full periods.
- **Method (§OS3): the must-fail control caught a vacuous test of mine.** The
  tail-confinement check has only one possible answer at |I| = 1; the control
  passed when it had to fail. Re-scoped, numbers unchanged, conclusions
  narrowed. Third instance of this failure mode in the project.
- **Infrastructure (§OS5):** int32 permutation replay (the *replay*, not the
  FWHT, is what binds memory — int64 needs ~24 GiB at n=30 and fails),
  `accel.pullback_stats`, and exact/cached `lab.measure.support(exact=True)` /
  `lab.measure.stats`. All gated in `test_accel.py`.

**Left open by this item, logged not claimed (§OS4).** The missing fraction of
the C30 hyperplane, 1 − 2·density, falls ~2.2× per bit of modulus in the
n_exp = 2 series and does *not* in the n_exp = 3 series; matched pairs show
N = 33, 35, 77 dropping 8–11× between n_exp = 2 and 3 while N = 21 barely
moves. Cheap follow-up: sweep n_exp at fixed N ≥ 33 with β > 1 and see whether
the deficit really has a cliff there.
