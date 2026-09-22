---
id: TX2
field: "Symmetric cryptography: linear cryptanalysis"
status: obstruction
effect: lower-bound
one_line: "Published nonlinearity bounds PPS term count from below; PPS hardness = linear-cryptanalysis resistance"
source: "Nonlinearity/Parseval bound as recorded in C25; dictionary in C12"
claims: [C12, C25, C26]
notes: [X]
todo: []
---
# TX2 — Symmetric cryptography: linear cryptanalysis

## Dictionary

Walsh spectrum / linear approximation table of an S-box or cipher ↔ PPS support of
the pulled-back bit.

## Hypotheses

Applies to any permutation. The bound is weak away from the extremes (C26).

## Consequence for the goal

Obstruction for every Walsh-basis method on high-nonlinearity functions. It says
nothing about other bases, which is why non-Walsh representations (TX5, TX7, TX9)
are the escape routes.
