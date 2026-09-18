---
id: 60
state: open
title: "Is the omitted-carry wrong reference a discriminating control, or is its disagreement forced by parity?"
outcome: ""
claims: [C72, C92]
---
# Does the no-h control test the carry, or only that something changed?

## The question

`experiment_dirty_prefix_walsh.py` builds a `no_h` wrong reference by passing
h = 0 into the arithmetic while the outer character keeps the physical h, and
registers its disagreement as a must-fail control (P4); `experiment_merged_
prefix_formulas.py` uses the same construction. Note DM reads its zero as
directly exposing the reference defect. The question is whether that zero is
evidence about the incoming carry at all.

## Premise, unverified and from outside the board

Calibration trial 1 (`SWARM_REVIEW_2026-09-12.md`, beside section 5) produced,
as a by-product, two enumeration claims on the N = 7, a = 1, n_exp = 1, q = 3
instance: the no-h reference is identically zero on every input mask containing
the h wire, and every carry-sensitive mask contains h. Its own probe is in
`out/agent-board/reviews/calibration-2026-09-17-latent-verifier/work/`, which is
gitignored, unreviewed and produced with no board assignment. It is exploration
and a starting point, not evidence: derive or re-enumerate before using it.

If it holds, the control's disagreement is parity-forced. It then shows that the
two references differ somewhere on an h-containing mask -- which any such mask
would show -- and not that the incoming carry matters for the queried
coefficient. That is METHOD.md's "a test with only one possible answer", and the
discrimination the submissions relied on actually came from the true coefficient
being nonzero.

## First step

Cheap and exact, at the existing fixed instance, before anything else: confirm
or refute both enumeration claims from the construction, and decide whether any
carry-sensitive mask avoids h.

- **Forced.** Replace the control with one that is not parity-annihilated --
  misroute the carry rather than omit it, so the wrong reference is a different
  nonzero function -- register the new prediction, and check the wording of C72,
  C92 and note DM where they read the zero as carry evidence. No measured
  coefficient changes: the true values are unaffected.
- **Not forced.** Record the mask that separates them and keep the control,
  naming that mask as the discriminating one.

## Do not

Re-run the dirty-prefix campaign, extend it to new sizes, or treat the trial
referee's numbers as an accepted result.
