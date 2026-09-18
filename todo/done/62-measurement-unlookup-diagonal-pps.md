---
id: 62
state: done
outcome: "C105/UL: measurement-based unlookup of the lookup register acts on diagonal PPS as the pullback that clears it (reset or deallocate; drop is measure-and-leave); same full-space function as the unitary circuit on every input with a clean lookup register; dead tail proved, 2-periodicity derived; off that set the functions differ (witness traced by hand) and the channel-model support sizes are not reported; Paper B 7.3 / 12.1 / 13 item 3 updated"
title: "Does Gidney's measurement-based unlookup keep diagonal PPS exact and branch-free, closing the unitary-unlookup gap in the windowed analysis?"
claims: [C8, C36, C37, C104, C105]
---
# Measurement-based unlookup under diagonal PPS

Chosen from the SN slate (2026-09-18), surveyor candidate R4. WD and PAPER_B §7.3/§13
item 3 record that the windowed analysis assumes a **unitary** unlookup, and that
Gidney's measurement-based uncomputation (X-measure the lookup register, fix the phase
negations with a smaller diagonal lookup, reset) is "the one real gap".

Candidate (surveyor's sketch, unreviewed): the channel's Heisenberg adjoint maps a
Z-string Z^{z_a} Z^{z_d} O to Z^{z_a} ⊗ I_d ⊗ O. The fixup is diagonal, the X outcomes
sum to the identity, and the reset gives <0|Z^{z_d}|0> = 1. So the rule is "clear the
lookup-register key bits, keep the coefficient, merge collisions". Diagonal observables
stay diagonal, and C8 survives with a non-injective population map in place of the
basis permutation.

Open before any test: reset against deallocate (the latter drops terms rather than
merging them); whether C17, C82 or C36–C39 lean on step-wise bijectivity; whether C104's
Lemma S carries over (its "Roles" clause fails for the windowed circuit, since `_activate`
applies X to window qubits, `windowed_arith.py:94-102`); and a must-fail control showing
that the channel and unitary models differ off the valid subspace.
