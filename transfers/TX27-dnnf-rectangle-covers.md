---
id: TX27
field: "Knowledge compilation: DNNF, d-DNNF, SDD lower bounds via rectangle covers"
status: open
effect: lower-bound
one_line: "Candidate: a DNNF of the modexp bit needs about (r^2/N)^(1/3) size via balanced rectangle covers and a bilinear exponential-sum bound"
source: "Bova-Capelli-Mengel-Slivovsky IJCAI 2016 Thm 6 p.1011 (read in the body by surveyor slate 4, S757c3612601d4a6a)"
claims: [C107]
notes: [EF]
todo: []
---
# TX27 — Knowledge compilation: DNNF, d-DNNF, SDD lower bounds via rectangle covers

## Dictionary

DNNF of size s ↔ balanced rectangle cover of size s; within a slice f = bit0(uv mod
N) with u, v on opposite sides.

## Hypotheses

Unrefereed sketch; needs hypothesis H-eq (bounded multiplicity of exponent-position
sums mod r) and prime N. Composite N needs C109-type conditions.

## Consequence for the goal

Would extend the obstruction to non-ordered classes (d-DNNF, decision-DNNF, FBDD,
SDD). Exponent below the sqrt(r) baseline.
