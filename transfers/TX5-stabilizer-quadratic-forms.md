---
id: TX5
field: "Stabilizer simulation: quadratic form expansions and Toffoli splitting"
status: imported
effect: constant-factor
one_line: "Quadratic cells = known stabilizer decompositions; exponential in escaped Toffolis, constant-factor gains on modexp"
source: "de Beaudrap and Herbert arXiv:2109.08629 §§2, 6; Bravyi et al. arXiv:1808.00128 §4.3; Khesin and Ren arXiv:2108.02686 Thm V.1 (all read in C83)"
claims: [C83, C84, C85]
notes: [QC, PR, CP, ST]
todo: []
---
# TX5 — Stabilizer simulation: quadratic form expansions and Toffoli splitting

## Dictionary

Affine-support quadratic-sign cell ↔ one stabilizer-state term; a Toffoli escape ↔ a
stabilizer-rank split.

## Hypotheses

Exact; the local closure test equals generic tableau rank (C84).

Montanaro's hitting-set algorithm for the path-sum cubic (arXiv:1607.08473, Props 13-14, read in slate 3b) is this row's Toffoli splitting in another form; not a separate transfer.

## Consequence for the goal

Cost is exponential in the number of escaping Toffolis unless recombination is
structured. TODO 67 closed on this row: no predicted advantage over TX7 on modexp.
