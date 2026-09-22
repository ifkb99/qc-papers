---
id: TX30
field: "Diagrams over linearly transformed variables (LTOBDDs, GF(2) frames)"
status: open
effect: unknown
one_line: "OBDDs after a GF(2) change of variables (C82's frames) are outside C107/C108; Sieling's lemmas would need statements about every GF(2) matrix"
source: "Sieling, LTBDD lower-bound lemmas 1 and 5 (surveyor slate 4); Gunther-Drechsler LTBDDs (not read)"
claims: [C82, C107]
notes: [EF]
todo: []
---
# TX30 — Diagrams over linearly transformed variables (LTOBDDs, GF(2) frames)

## Dictionary

OBDD in a GF(2)-transformed frame ↔ C82's lazy Walsh frame on keys.

## Hypotheses

No estimate. The directional-affine-extractor route fails because f = 0 on the
subspace {x = 0}. A frame escape would need <a> to act GF(2)-affinely, known only
for order <= 2n on N = 2^n - 1.

## Consequence for the goal

An open escape; a lower bound here would close the frame-based variants of the
propagator.
