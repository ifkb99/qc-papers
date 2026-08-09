---
id: 12g
state: done
title: Does the hyperplane deficit fall off a cliff between n_exp 2 and 3?
claims: [C7, C21, C30]
---
# Does the hyperplane deficit fall off a cliff between n_exp 2 and 3?

**NO — closed negatively 2026-08-08, and §OS4's observation is RETRACTED.**
See `NOTES.md` §DF; file `experiments/experiment_c7_deficit.py`, which **fails
3 of 4 predictions and exits nonzero by design.** No new claim; one is removed.

- **D is monotone decreasing in n_exp**, over ~40 consecutive steps across four
  β > 1 instances at fixed (N, a). No oscillation, so no period, so the
  §I period-ord₂(β) structure — the prior that made this look worth doing —
  **does not transfer from function level to circuit level.**
- **The "cliff" was two points straddling one instance's one anomaly.** Every
  instance has the same step-ratio profile (tiny 1→2, one more sizeable drop,
  then a slow climb towards 1); N = 21 is the only one whose second drop lands
  a step late (.021 .936 .174 vs .045 .275 .734 elsewhere), and n_exp = 2 vs 3
  is precisely where it disagrees. Not α (N = 33 shares α = 1), not β (N = 13
  shares β = 3), and not pursued further.
- **Control clean:** β = 1 at the same modulus locks at n_exp = 2 (C21) and D
  climbs monotonically to 0.985, no oscillation — so "no oscillation" above is
  a fact about the circuits, not about the instrument.
- **Logged, not claimed:** the step ratios climb towards 1, so D appears to
  converge to a small instance-dependent constant rather than to 0 — density
  tending just below ½ as n_exp grows at fixed N. That is extrapolation from a
  converging sequence, which is the exact move §I punished. It does not touch
  C7, which is about growing N and is measured. The honest route if anyone
  wants it is a function-level n_exp → ∞ calculation, not more widths.

<details><summary>Original note, kept for the record</summary>

Opened by §OS4, and cheap. C30 confines the support to a hyperplane, so the
interesting quantity is how much of that hyperplane is *missing*:
1 − 2·density. Across the n_exp = 2 series it shrinks with striking regularity,
~2.1–2.3× per extra bit of modulus (0.0544 → 0.00112 over n = 3..8). Across the
n_exp = 3 series it does not. The matched pairs are the odd part: same N, same
a, only n_exp differing, and N = 21 barely moves (0.01074 → 0.01004) while
N = 33, 35 and 77 each drop **8–11×**.

Either there is a threshold in n relative to n_exp, or it is instance noise
across four points. **Design:** fix N ≥ 33 with β > 1, sweep n_exp by 1 (this
varies exactly one parameter, unlike the two series, which vary N), and look at
the deficit. q = 3n + 4 + n_exp, so N = 33 reaches n_exp = 8 within budget.
Must-fail control: a β = 1 instance, where C21 locks the support and the
deficit must therefore blow up rather than shrink.

Low stakes — it refines the constant in C7, it does not touch Θ(2ⁿ) — but it is
one sweep and the regularity in series A is too clean to leave unexamined.

</details>
