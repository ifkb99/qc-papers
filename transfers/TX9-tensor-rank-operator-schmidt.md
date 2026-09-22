---
id: TX9
field: "Tensor networks and operator Schmidt rank"
status: imported
effect: lower-bound
one_line: "Cut rank is basis-invariant and lower-bounds every representation linear across the cut"
source: "Standard local-basis invariance (C48); Dang, Hill and Hollenberg arXiv:1712.07311 §§4, 5.2 (C19)"
claims: [C19, C48, C106]
notes: [TR, CR]
todo: [68]
---
# TX9 — Tensor networks and operator Schmidt rank

## Dictionary

Matricization of the sign tensor across a cut ↔ TT/MPO bond, automaton states, K2
leaves.

## Hypotheses

Exponent-prefix cuts only for C106; other cuts are TODO 68.

## Consequence for the goal

Clean-code scalar functions have rank at most r; full-space ranks saturate on dirty
scratch (TODO 68 premise check). Removes nothing on the full space; says the clean
contract is the place to look.
