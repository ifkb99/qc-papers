---
id: 68
state: open
title: "Cut rank and bond at cuts that are not exponent prefixes: interleaved and work-first orders of the modexp pullback"
outcome: ""
claims: [C19, C48, C103, C106]
---
# Beyond exponent-prefix cuts

From a 2026-09-18 discussion. C106 bounds the bond of every representation linear
across an **exponent-prefix** cut and says explicitly that it covers no other cut.
Tensor-network memory (the one simulator family where entanglement itself is the
cost) is set by the worst cut of the chosen order, so other orders are open:
work-first (note SL's S3 predicts ROBDD width ≤ 2^(q_w) at exponent levels),
interleaving each exponent bit with the work register, and orders that place clean
scratch at the ends.

## Questions

* Is there a derivable upper bound on the rank at a cut (work bits | exponent bits)
  or at an interleaved cut, from the prefix-product structure of C103 Lemma 0?
* Does the order that minimizes the maximum cut rank also minimize bytes of a TT/MPO
  with the same exact output? (C106 bounds bond, not bytes.)
* How does TODO 66's clean-scratch restriction change these ranks? C48 already gives
  rank ≤ r for the ideal scalar function, so the restricted object may be easy and
  the full-space one hard; say which object each number is about.

## Premise check (plan)

Note TR measured every natural-order cut for several rows; note HB lists "the same
bound in other variable orders" as open. Read both before computing, and do not
re-measure natural-order cuts TR already owns.

## Premise check (2026-09-18, arithmetic on TR's saved `out/tensor_memory.json`)

* **Natural order: the cost sits in the dirty scratch block.** Qubits 0-7 are b and
  t. At full N = 7, t = 4 the ranks are 2, 4, ..., 128: every cut through b|t is at
  its maximum 2^k. The rank peaks at the t|x boundary (217 for a = 3, 93 for
  a = 6), then falls to 2-8 across the exponent register (consistent with C106). Any order that cuts through b and t pays this;
  the clean restriction of TODO 66 (note VR) removes b and t from the object.
* **Bytes.** A natural-order TT stores sum 2·r_i·r_(i+1) numbers: 58,124-143,840
  against Walsh support 15,539-64,353 (a = 3, 6; t = 2..4). That is 2.2-3.7x more
  numbers, so a TT is not a memory win at these sizes. For a = 3 it grows about
  1.5x per exponent bit against Walsh's about 2x, so a crossover in t is
  possible. It has not been derived or measured.
* **Interleaved orders are partly answered for one block.** C102's u_a ROBDD
  (order exp, c0, anc, then t_i, b_i, x_i interleaved) is about 10·4^n nodes, and by
  an elementary count (rank ≤ distinct subfunctions across the cut ≤ total ROBDD
  nodes; C106 Theorem A states the exponent-cut case) interleaving bounds every
  per-block cut rank by about 10·4^n, against 2^(2n+2) through b|t in natural order.
  Nothing is measured for multi-block modexp in that order.

**Narrowed question.** Across t blocks, does the per-cut rank in C102's interleaved
order stay about 4^n (bounded by the ROBDD widths) or multiply with t? The fair
comparison is TT bytes against C104's K2 in bytes for the same exact output.
Deprioritized behind TODO 66's block split: on the clean contract the object is
small in every representation, and on the full-space contract this refines a
constant.

## Checkpoint (2026-09-18, note TL)

Answered at the level asked, pending review: deriver 4k (`S11b9708e34d14abe`) proved a
rank identity for arc-translate matrices and, modulo TX8 Theorem 1, a rank lower bound
N^η/(2 C1 log N) at an x-internal cut in **every** variable order, so every TT/MPS/MPO
bond is exponential in n there for typical prime N (register TX17). Remaining here: a
referee of the identity; composite N via TX25; bytes rather than bond.
