---
code: EF
title: "ROUND AFTER TL. C107-C109 integrated after five integration reviews; slate 4 maps the escapes: they need structure (a in ±<2>, cofactors of 2^M ± 1), not small order, and every ordered-linear bound leaves DNNF, LIMDD and GF(2)-frame diagrams open."
outcome: record
claims: [C107, C108, C109]
todo: [68]
---
# EF — Round after TL (2026-09-18): integration of C107–C109, and slate 4 on the escapes

The user chose all three items proposed at the end of note TL: the author's revision of
TX8 Theorem 1, fresh referees for the rank identity and the composite-modulus bounds,
and a slate on the remaining escapes.

## Integration: C107, C108, C109

* Revision v2 of the TX8 derivation (`S6252c91848bd4c90`) applied R1–R9; a fresh referee
  accepted it (`V6a01bd7d375d4b39`, Q1–Q7 at integration). The rank identity was accepted
  with D1–D6 (`Sd4b9709071a34341`); the composite-modulus propositions were accepted,
  their restatements corrected (`S57b82692609045c0`, D1–D8).
* Integration task `T49262ab09c0845b0` took six submissions and five reviews before
  acceptance. Every blocking finding was in the coordinator's integration prose, never
  in a proof:
  * v1: a referee's finding rewritten into a different finding (twice), an undefined
    constant that parsed as a claim ID ("C1"), dropped hypotheses, an unconditioned
    title, a one_line repeating a defect its own body fixed.
  * v2–v4: the same class recurred next to each fix (B6 → N1 → R1): sentences stating
    the status ("open", "small") of an unproved class, each false just past the previous
    detector's fixtures; the last came from transcribing an earlier review.
  * v4–v5: provenance: a wording that misdescribed a sample, an unarchived check cited
    as evidence.
  * The fix that converged: state only proved examples and an explicit no-claim, and
    detectors that sweep all files (status words; every "checked/measured/verified"
    resolves to a run, a log or a proof; every evidence cell quotes its source line).
* The board refused closure after an "execution-time" correction changed accepted
  bytes; the correction went through as v6 on a reopened task. Integrated bytes equal
  accepted bytes.

## Slate 4 (surveyor `S757c3612601d4a6a`, deriver `Se2325ec074c0487a`)

* **Premise corrections.** The bounds cover ordered *linear* representations; Pauli-LIMDDs
  and diagrams over GF(2)-transformed variables are outside the rank argument. Ponzio and
  Bollig–Woelfel do not apply to f as theorems (two free operands); only their methods
  transfer. The formula equals the circuit's clean bit only where x ⊕ 1 < N (C107's
  object states this).
* **Candidates (unrefereed).** Surveyor S1: a DNNF / rectangle-cover bound
  Ω~((r²/N)^(1/3)) under a hypothesis H-eq (would cover d-DNNF, FBDD, SDD); S2: a
  sign-rank bound at every balanced cut (vtrees, binary tree tensor networks); S3/S4:
  LTOBDDs and LIMDDs open, no lower bound known. Deriver S1: a digit identity
  (parity(2^j y mod N) = bit_{−j mod m}(D_N y)) and a cofactor-escape family; S2: every
  exponent-first diagram has ≥ r nodes and every x-first one ≥ the largest prime power
  dividing r, so escapes in those order classes would give order finding; S3: small order
  is not sufficient (order-3 multipliers, a = 2 at unstructured N, measured n ≤ 20).
* **Baseline caution.** Every derived exponent (1/3, 1/6, η ≤ 1/6) sits below the √r cost
  of baby-step giant-step: these bounds rule out polynomial size, not beating the generic
  group baseline.

## The sealed candidate

Hash posted before dispatch (`Mba07396fce81467b`), text revealed in `M97896ec01279459`.
K0 conjectured "polynomial representation iff small r". The "only if" direction matches
the deriver's S2 in exponent-first and x-first orders; the "if" direction is refuted by
its S3. Half survived, the first coordinator candidate to do so in three rounds.

## Process

* Integration prose is where errors entered this round, not the proofs; see above.
* The v4 referee disclosed that it proposed the v5 wording it then re-checked; a fresh
  referee was not used for that focused step.
* A surveyor's downloads landed in the session's tool-results directory, outside its
  write area, and the guard blocked the cleanup.
