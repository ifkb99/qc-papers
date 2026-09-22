---
id: TX24
field: "Large-domain #CSP: Mal'tsev and coset tractability"
status: barrier
effect: barrier
one_line: "Shor's constraints are coset (affine-type) only in an order-aware presentation, which is TX15 in CSP language"
source: "Cai-Chen arXiv:1111.2384 Thm 1 p.2, Lemma 19 p.20, §3.2 p.21 (slate 3b K4)"
claims: [C52]
notes: [ST]
todo: []
---
# TX24 — Large-domain #CSP: Mal'tsev and coset tractability

## Dictionary

Register-valued variables over Z_N^* ↔ controlled multiplication as the relation R_c
= {(e, y, y'): y' = c^e y}; coset constraints over a finite abelian group ↔ the
affine-type tractable class.

## Hypotheses

R_c is not closed under the group Mal'tsev operation unless c² = 1 (derived; other
polymorphisms unchecked). The coset form needs r | M: Shor's hidden-subgroup
structure. Mal'tsev counting is polynomial in the domain size, which is N here (to
confirm).

## Consequence for the goal

Any affine or holographic success on Shor must supply the order or a group
presentation.
