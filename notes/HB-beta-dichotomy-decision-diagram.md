---
code: HB
date: 2026-09-15
title: "The beta dichotomy through a decision diagram: prefix-product counting"
outcome: exponent-first ROBDD widths of C15's object bounded by 2^alpha D_mu(beta)(k-alpha) via C101's block inverse and a group counting lemma; attained on 11 fixtures including 5 out of sample; final Walsh support doubles per bit on the same fixtures; two referee rounds of write-up corrections, the second adding the alpha+mu crossover caveat
claims: [C103]
todo: [13]
---
# HB — How the μ(β) bound was found

C103 owns the statement. This note records the path.

## Why

C102 (TODO54) compared Walsh support with the ROBDD of one compiled
multiplier. The user then chose to keep researching before review. TODO13
asks whether a third simulation method keys on Paper B's r = β·2^α
invariant. The ROBDD of C15's exact object, on the full dirty space and in
exponent-first order, was the cheapest same-object test. It ran as board task
T58ddf9e611124d65, in proposal mode because Ta294b7d88a4c4015 held the write
scope.

## Path

1. **Before measuring**, the β = 1 (C23 parity) and β = 3 (C101 block
   inverse: every prefix product is M_c^s H) bounds were derived and
   registered (Mc3a0e636ff9741de). The first write task, T90e1b7284a674a6c,
   was cancelled over the write-scope conflict and replaced by the proposal
   task.
2. **run_v1** confirmed both bounds with equality. For β = 5 the widths were
   2(F_{k+2} − 1) at N = 11, a = 2. That was noticed post hoc, then
   predicted for a = 3 before its rows existed (M0b80d942f7e64a90). The same
   message contained an arithmetic slip (174/282/454 for 176/286/464), which
   was corrected before those rows existed (Md5afaf1c5e234e33). Both numbers
   stay on record.
3. **The counting lemma and μ(β)** were derived after run_v1 (derivation_v1).
   A garbled board restatement was replaced by M5572e29dc56e49d7. An early
   claim that β = 7 gets only 2^k was wrong: plain periodicity gives μ = 3.
4. **mu_run_v1** registered equality on five untested moduli
   (M108c659b84254901) and passed.

## Review round 1 (Vc7ce69ad6e4b4fe1, changes_requested)

The referee confirmed the lemmas, the theorem and every logged number, with
independent q ≤ 20 recomputation and an 800-case group check. It required
write-up fixes, all made in derivation_v2.md and draft_C103_v2.md:
* **R1:** the control-0 identity needs the cancellation argument, since the
  multiply half moves about 96% of states.
* **R2:** the step from Lemma 2 to D_μ needs the increment induction.
* **R3:** the β = 1 t-independence argument.
* **R4:** a false sentence on work-part scaling.
* **R5:** bound language, the finite-group ceiling, and precise
  "basis-dependent" wording.
* **R6:** relations to C19, C45, C48 and Paper B §11.1.
* **R7:** this note was missing.

It also observed that mu_run_v1's controls follow from M2's equality, so they
are not independent fault detectors.

## Review round 2 (V2e9e6c322bee4511, accept)

A fresh referee verified Lemma 0 against the actual gate list rather than the
prose, tracing that the constant reaches `cc_add_mod` only through
ctrl-controlled `_load` Toffolis and that `cswap` collapses to two equal CNOTs
at ctrl = 0, then confirming numerically at N = 7 and 11. It re-derived the
exponent-level widths for six fixtures without using `lab/bdd_count.py`, and
checked the prediction record precedes `run.start`.

Two corrections were required before integration, and both are in C103 rather
than in the frozen derivation:

* **The scratch is unloaded, not empty.** The v2 wording said "the adders of
  the empty scratch". On the dirty space the scratch holds arbitrary values;
  were it empty, C_0 would be the identity, which the derivation's own
  measurement (7872 of 8192 states moved) refutes.
* **The crossover was missing**, and it is the one caveat that changes how the
  headline should be read. Below it the bound degenerates to the trivial 2^k,
  and since μ(β) can reach (β−1)/2, at large β most of the measured range sits
  there. C103's Limits own the exact statement and the counts; they are not
  repeated here.

The second is a good illustration of the standing warning that a result
agreeing with the hypothesis is the dangerous kind: every fixture satisfied
the bound, the bound was attained, and the agreement was partly because the
bound had degenerated to 2^k.

It has a second lesson attached. Review `V9b418bf68edf4bfd`, refereeing the
integration, found the crossover stated one level early — D_μ(μ) = 2^μ
exactly, so level k = α + μ is trivial too and the first strict improvement is
k = α + μ + 1, making the trivial count exactly min(t, α + μ + 1). That
off-by-one came from the accepted review's own wording of the correction, and
the coordinator carried it into the claim without deriving it. A referee's
formula is no more citable-without-checking than any other number in this
repo; the corrected version was checked against all five out-of-sample
fixtures before it was written down.

## Integration

Board task Tf102392bd0124c1f, replacing Tbe30548d3e774494. The replacement was
needed because the board's recursive dependency check went stale on
prose-only drift in claims/C102.md (the accepted C102 v3 corrections); the
accepted result is bound here by hash-sealed inputs instead, with the override
and its cost recorded in Mf2ab5cb32c584940.

## Candidate directions, not opened
* A proof of equality: that the compiled blocks act freely on these words.
* The same bound in other variable orders, or with the clean-input
  restriction (C45's setting).
* Whether any μ-type invariant appears in the post-QFT output distribution's
  cost.
* Where the α + μ crossover actually sits for the β that matter to Paper B,
  and whether any fixture reaches t far enough past it to show λ_μ growth
  directly rather than by extrapolation.
