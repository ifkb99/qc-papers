---
code: SL
title: "SLATE ROUND 1. Peak memory of exact PPS on Toffoli modexp: the output contract decides."
outcome: record
claims: [C24, C45, C101, C103]
todo: [55, 56]
---
# SL — Slate round 1 (2026-09-16): peak memory of exact PPS on Toffoli modexp

First use of METHOD.md "Where ideas come from" and SWARM.md phase 0. Two fresh
contexts that could not see each other generated candidates: a surveyor
(transfer, obstruction; board submission See53b1157d2844cd) and a deriver
(structure inventory, exact symmetries, sequences; Sc7408765323745c1). The
coordinator sealed its own candidate by hash before either returned.

## Premise: the output contract

"Same output as `propagate_perm`" was unstated, and both slates found it
decisive. The final operator of the full-space modexp pullback has about
2^(q−1) Walsh terms (C30's half density).
* **O1** materialized dictionary: every method floors at the output, and the
  peak term count (~0.748·2^q measured on N = 7) can shrink by at most ~1.5×.
* **O2** streamed exact output: the Heisenberg peak above the output is removable.
* **O3** exact compressed output: the output itself can shrink.

## Candidates and choice

Chosen **O2 with S1** (TODO 55); **S2** queued under O3 (TODO 56).

| id | candidate | origin | contract |
|---|---|---|---|
| S1 | exponent-slice decomposition; each slice a C45 reduced observable | both slates, independently | O2 |
| S2 | exponent-class form with Krawtchouk mixing, β ∈ {1, 3} | deriver | O3 |
| S3 | decision diagram with work variables first, width ≤ 2^(q_w) at exponent levels | surveyor (not independent, below) | O3 |
| S4 | carry C103 Lemmas 1–3 to the mid-propagation Walsh-coefficient diagram | coordinator, sealed | O3 |
| S5 | orbit compression under exponent-bit permutations (symmetry proved at t ≤ 3) | deriver | O1/O3 |
| S6 | exact integer tensor train, ranks from C48 cut ranks | surveyor | O3 |
| K2 | dense work-state recursion (fair baseline at fixed N; labelled B1 in this note until 2026-09-17, and K2 in the TODO 55 derivation, where B1 means per-slice C45 runs) | both | baseline |
| X1 | the ~3/4 peak cap as one broken linear structure at the top-carry Toffoli | deriver | O1 |

## Killed, with reasons

* **Multiplication OBDD lower bounds as an obstruction here.** Woelfel's bounds
  (read in the body) need growing operand width; at fixed N they are constants,
  and they bound a Boolean OBDD of a product bit, not Walsh-coefficient
  functions at a propagation peak.
* **Reordering exponent blocks to lower the peak.** Changes it by at most 0.2%
  at t ≤ 3.
* **S4 for peak memory under O2.** S1 reaches a linear-in-t peak with a proved
  premise and existing code.

## Unexplained (deriver)

The first item was closed as a coincidence of totals in note SN (2026-09-18).

* N = 9, a = 4 and a = 7 at t = 3 have equal support counts (255,973) and
  different supports (symmetric difference 3,616).
* N = 7, a = 2, 4 at t = 2: the prior step's support violates the b_msb ⊕ anc
  hyperplane, yet the 3/4 cover still holds.

## Process lessons

* **Convergence.** The top candidate was found independently by both
  generators; it had sat beside C45 through the depth-first chains.
* **Blinding leaked through examples.** The coordinator withheld decision-
  diagram paths, but METHOD.md's new examples named that work; the surveyor
  disclosed it, so S3 does not count as independent. Blinding must cover
  examples in recently edited permitted files.
* **Re-running a worker's script can overwrite its evidence** when the script
  writes to a hard-coded attempt path. The coordinator's reproduction was
  stopped before its final write; the frozen submission was unaffected.
  Re-runs go through a copy whose outputs point elsewhere.
