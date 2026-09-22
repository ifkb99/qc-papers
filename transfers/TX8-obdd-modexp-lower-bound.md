---
id: TX8
field: "Branching-program lower bounds; exponential sums over multiplicative subgroups"
status: obstruction
effect: lower-bound
one_line: "Every OBDD of the modexp x0 bit is exponential in n in every variable order when kappa(<a>) is small (C107); small-order structured families escape"
source: "C107 (refereed derivation v2 S6252c91848bd4c90); external inputs verified by Sf8eb4afd15a84c9e"
claims: [C107, C109]
notes: [SL, SN, ST, TL, EF]
todo: []
---
# TX8 — Branching-program lower bounds; exponential sums over multiplicative subgroups

## Dictionary

Cut count of an OBDD at a split of x's bits, after restricting e and scratch ↔
distinct rows of 1_I(c(alpha+beta) mod N) (C107 Lemma P). The formula equals the
circuit's clean bit where x XOR 1 < N (C107's object).

## Hypotheses

As stated in C107 (hypotheses there; not restated here). Composite N: C109. Escapes:
C107 Theorem 2 and its composite-N examples; slate 4 (note EF) finds escapes need
structure (a in ±<2>, cofactors of 2^M ± 1), not small order alone.

## Consequence for the goal

Closes ordered bit-level diagrams for typical moduli. Remaining diagram routes:
TX27-TX30.
