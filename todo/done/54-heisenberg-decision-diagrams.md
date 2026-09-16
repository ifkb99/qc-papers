---
id: 54
state: done
title: "Is the Walsh (Pauli) basis an asymptotically wasteful store for the pulled-back arithmetic observable, compared with a canonical decision diagram?"
outcome: "C102/HD: adder Walsh 3*2^(m-1)-2 vs ROBDD 3m+4 proved; u_a Walsh ~2^(3n+3) vs median ROBDD 9.7-14.3*4^n measured to n=8 (N mod 4 dependent) with a slowly growing random-circuit control; no u_a bound, no escape from the order barrier"
claims: [C7, C8, C17, C19, C48, C82, C102]
---
# Heisenberg decision diagrams against Walsh support

Completed at bounded scope. C102 owns the proved adder counts and the
measured u_a scaling with its limits; HD owns the path, checks, preserved
failures and weak points. No follow-up TODO was opened; HD lists candidate
directions for the user to decide on. The pilot expectations below are the
original brief, superseded by C102.

A Claude-originated question, raised on 2026-09-15 at the user's invitation to
pick a direction on the CNOT/Toffoli memory problem. The user has since
authorized spending the remaining usage on it, with the arb board kept current.

C8 fixes PPS storage for a Z-observable through a basis permutation as the
Walsh sparsity of the pulled-back output bit. That is one representation of a
Boolean function. A reduced ordered binary decision diagram (ROBDD) of the SAME
full-space function is canonical and supports gate-by-gate composition, so it
is a same-object, same-exactness store. The question is how the two scale on
the repository's compiled arithmetic, and where the difference comes from.

Scratchpad pilots (not evidence; predictions for the registered run are
frozen on the board before measurement) suggested:
* Cuccaro adder top sum bit: Walsh grows exponentially in the width m,
  while the interleaved-order ROBDD and an affine-subspace decomposition grow
  linearly. A natural (register-block) order is exponential.
* `ToffoliModExp.u_a` (full dirty space, observable x0): the Walsh/ROBDD
  ratio roughly doubled per bit of N for n = 3..6, while a random reversible
  circuit with the same gate multiset kept a near-constant ratio.

Scope to keep: the full-space permutation, final support versus intermediate
peak, node and term counts are not allocated bytes, and no representation here
escapes the order/table-size barrier of C19/C48 for modular exponentiation.
Decision-diagram simulation of arithmetic is established technology (TODO13
names MQT DDSIM); any novelty is limited to the exact comparison with the PPS
object and its proofs, and needs a primary-source audit before being asserted.
