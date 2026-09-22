---
id: TX41
field: "Word-level arithmetic diagrams and zero-set extraction"
status: open
effect: unknown
one_line: "Compact integer multiplication diagrams do not yet supply efficient construction or counting of the fixed-product Boolean zero-set"
source: "Bryant and Chen (1995), section 3.1/Figure 4 and section 5"
claims: [C118, C119]
notes: [FC]
todo: [72]
---
# TX41 — Keep the arithmetic expression compact, then find its zeros

## Dictionary

The numeric word function F(p,q)=pq-N can remain a compact arithmetic object.
Factoring needs assignments where F=0, with bounds, or their constructive
conditional counts. Those operations are distinct from representing F and
from testing equality of two entire arithmetic functions.

## Hypotheses

[Bryant and Chen](https://www.cecs.uci.edu/~papers/compendium94-03/papers/1995/dac95/pdffiles/32_1.pdf),
section 3.1 and Figure 4, construct linear-size word-level sum/product *BMDs.
Section 5 explains that Apply need not have a polynomial number of recursive
calls. Thus compact representation does not alone charge all operations.
TX26's modexp pullback bounds concern another numeric function and are not
used to rule out compact multiplication. C118 obstructs only its specified
residual processing rules, leaving stronger zero-set operations open.

## Consequence for the goal

The next useful lemma must construct or simplify the bounded zero-set from N
without factor-dependent advice. Expression size or a known-factor diagram
is insufficient. This round establishes neither a new factoring method nor a polynomial
constructor for the zero-set.
