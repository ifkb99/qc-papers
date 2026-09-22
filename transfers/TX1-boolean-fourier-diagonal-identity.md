---
id: TX1
field: "Analysis of Boolean functions"
status: imported
effect: exact-cost-model
one_line: "Diagonal Pauli expansion = Walsh transform of the pulled-back bit function"
source: "Montanaro and Osborne, 'Quantum boolean functions', Proposition 9 (credited in C8 and Paper A §9)"
claims: [C8, C17]
notes: []
todo: []
---
# TX1 — Analysis of Boolean functions

## Dictionary

Their Fourier expansion of a Boolean function ↔ our PPS term dictionary for a Z-type
observable pulled back through a basis permutation.

## Hypotheses

Holds for any permutation circuit and any computational-basis observable; C8 states
the object. Nothing assumed.

## Consequence for the goal

Gives an exact cost model, so every later Walsh-side theorem applies to PPS. It does
not remove any exponential; it relocates the question to Walsh sparsity.
