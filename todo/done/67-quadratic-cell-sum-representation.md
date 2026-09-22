---
id: 67
state: done
title: "Store the Heisenberg operator as a sum of affine-support quadratic-sign cells instead of Walsh keys: how many cells does modexp need?"
outcome: "Closed at the premise check (2026-09-18): the ledger already answers it. C102/HD measure affine pieces on u_a (a constant factor below Walsh, beaten by the ROBDD), C83/C84 and TODOs 45-47 identify quadratic cells with known stabilizer decomposition, and DN shows diagram-native propagation loses in bytes"
claims: [C8, C25, C83, C84, C95, C102]
---
# Quadratic-cell sums as the memory unit

## Closed at the premise check (2026-09-18)

No compute was run. The ledger already covers each part:
* **Affine pieces through u_a.** Note HD step 2 used `lab.affine_pieces` with
  sibling merging: 732 and 6,081 pieces at N = 7 and 11 against Walsh 3,206 and
  31,022. That is a constant factor, the counts depend on the branching rule
  (C102), and the ROBDD in a good order beats them (HD step 3).
* **The adder.** The "2m cells" derivation planned below is C102's affine-piece
  bullet, already proved.
* **Quadratic cells.** C83's closure test equals generic stabilizer tableau rank
  (C84, TODO 45); coherent recombination after an escape is TODO 46/47 (C85, C86),
  and its isolated-adder novelty is closed by known baselines.
* **Diagram-native propagation.** Note DN: about 2x fewer nodes than a matched
  null but 12-41x more bytes than the dictionary.
The idea's one untested residual, a quadratic (not affine) cell count through u_a,
is not worth opening: C84 makes it a stabilizer-rank count under known splitting,
and nothing in the ledger predicts it beating the ROBDD. The original brief follows.

From a 2026-09-18 discussion. Walsh-key memory is Walsh sparsity (C8), which C25
bounds below by nonlinearity: bent circuits reach full support 2^n. Yet a bent
quadratic function is a single affine-support quadratic-sign cell. If the operator
were stored as a sum of C83 cells, memory would count cells, not Walsh terms, and
a function of low "quadratic rank" would be cheap even when its spectrum is flat.

## Premise check first (no compute before it)

C83 and note QC already built exact one-cell closure for Toffoli and streamed
selected queries; QC's outcome was a loss to the recognized formula for selected
queries. C95 shows native final-shift lifts need exponentially many convex sign
cells. Establish:
* whether QC/C95's objects are this one (a whole-operator cell sum under reverse
  propagation, with cells split on closure failure), or a selected-coefficient
  certificate and a convex-cell lift, respectively;
* what split rule C83 implies when a cell escapes (rank-3 cubic term): the number
  of cells a split creates bounds growth per Toffoli, and should be derived;
* the strongest baseline: C104's K2 and C103's ROBDD for the same output contract
  (note SL's O1/O2/O3), with preprocessing and bytes charged.

## Cheapest discriminating step

Derive the cell count after one Cuccaro adder and after one cc_add_mod from the
construction, and compare with the Walsh count there. A must-fail control is a
random reversible circuit of the same Toffoli count, where cell count should track
Walsh sparsity. Only if the derivation predicts a separation on modexp is a run
worth it.
