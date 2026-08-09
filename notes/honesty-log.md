---
code: honesty-log
title: "Honesty log"
date: 2026-08-08
outcome: record
claims: [C17, C24, C29]
todo: [12, 12e, 9]
---
# Honesty log

Things believed and then killed, in order. Keep adding to this.

- H1 "π/4 angles break truncation" — killed by F1 (it all cancels).
- H2 "uniform angle breaks the power law" — killed by paper App. C.
- F2 "the QFT is the PPS bottleneck in Shor" — killed by the real modexp
  (bit-reversal bug + a toy where the observable never met the arithmetic).
- F4/F6/F8 "Fourier lacks Z-closure; Toffoli is 8.5x better; truncation
  reverses it" — **all killed by the θ=π PPS bug.** Both compilations have
  exact Z-closure; Fourier is cheaper, not dearer; the ⟨O⟩ values were wrong.
- F1's *explanation* — right conclusion (adder collapses to 1 term), wrong
  mechanism (permutation-ness). Real reason: the output bit is affine.
- F7 "Walsh sparsity is the driver, compilation sets the ceiling" — half right.
  Walsh sparsity is exactly the driver (F9), but there is no separate
  compilation ceiling; the apparent one was an ancilla-count confound.
- C17 "permutation-native PPS halves peak memory **exactly**" — **rounding
  artifact**, corrected 2026-08-08. The test table printed `%.1f`, so 1.9997
  displayed as "2.0x" and the phrase "exactly" was written from the display
  rather than the numbers. True: 2.000000 for adders, 1.9997 for modexp with
  `rot = 2·perm − 2` in 3/3. **The lesson is narrow and worth keeping: never
  write a precision claim from a rounded display.** Caught by external review
  noticing the abstract claimed more than the section delivered.
- C29 "V²=id is the entire condition" (§G, TODO 9) — **over-claimed, narrowed
  by §WD.** Not wrong about anything it tested; wrong about what it had tested.
  Every synthetic block in step 9 was controlled on exactly one fresh qubit, so
  the experiment could not separate "V is an involution" from "the control is a
  single qubit". `SelectModExp` has the first and not the second, and loses the
  invariance. **The generalisation was stated at the altitude of the proof
  sketch rather than of the evidence** — the proof genuinely uses only V²=id
  *given* the block form, and that conditional got dropped in the write-up.
- P5 (TODO 12) "restoring the j=0 branch makes the tail dead" — **half
  refuted, and the half that failed was informative.** Support flatness held
  exactly; the dead-bit half failed because the activation ancilla is itself a
  scratch qubit the pullback ranges over. Logged rather than re-scoped; the
  diagnosis is in §WD and is now a check in `experiment_windowed.py`.
- §OS4's "the hyperplane deficit falls off a cliff between n_exp 2 and 3"
  (TODO 12e → 12g) — **retracted the same day it was written, by the sweep it
  asked for.** D is monotone decreasing in n_exp in all four instances; the
  "cliff" was a two-point comparison landing on the one step where N = 21
  disagrees with every other instance. Believed because three of four matched
  pairs moved together and the fourth looked anomalous — when in fact the
  fourth was anomalous *at that step only*. **A two-point comparison cannot
  tell a trend from a phase**; §I had already taught this and it was not
  applied. Logged rather than quietly deleted because the observation was
  correctly flagged "not a law" when written, and the flag is what made the
  one-sweep check obvious.
- C24's tail-confinement test at |I| = 1 (TODO 12e) — **not wrong, but
  vacuous, and it had been counted as evidence.** "z_I ∈ {0, 1_I}" has one
  possible answer when the tail is a single bit, so the n_exp = α+1 rows say
  nothing and a β>1 control evaluated there *passes*. Caught only because the
  control passed. Re-scoped to |I| ≥ 2; no number changed, three of six rows
  stopped counting. See §OS3. **The general shape: a test can be sound,
  correctly implemented, and still carry zero information at some parameter
  values — and those are exactly the values a sweep starts at.**

**Lessons, in order of how much they cost:**
1. The toy sandwich (H → adder → QFT) is not a proxy for Shor. Use `modexp.py`.
2. **Randomised tests only cover the gates they sample.** The θ=π bug survived
   because no random circuit ever emitted an X. Enumerate the gate *set*, don't
   sample it.
3. **An independent exact reference is worth more than more tests.** Walsh
   caught what six suites missed, because it computes the same quantity by a
   completely different route.
4. **Precision-sweep to separate bugs from float error.** Same error in
   `longdouble` as in `float64` ⟹ it is a bug, full stop.
5. Matched-instance A/B needs matched *qubit counts*, not just matched logical
   function.
6. **State a generalisation at the altitude of the evidence, not of the proof
   sketch.** C29 said "V²=id is the entire condition" on the strength of blocks
   that all had exactly one control. When every instance you tested shares an
   incidental property, that property is a *hypothesis you did not vary* — list
   it before generalising. The cheap check is to ask what the must-fail control
   would have to look like to break the clause you are about to drop; here it
   was "a block with V²=id and a non-affine control", which took an afternoon
   to build and refuted the claim immediately.
7. **The pullback sees scratch qubits in every state, not just |0⟩.** Two
   constructions identical on the valid subspace can differ completely in Walsh
   support — that is the whole content of §WD, and it is also what refuted half
   of P5. Reason about the permutation on the *full* space, always.

---
