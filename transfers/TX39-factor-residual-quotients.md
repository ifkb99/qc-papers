---
id: TX39
field: "Constraint contraction and residual polynomial quotients"
status: obstruction
effect: lower-bound
one_line: "The specified low-bit multiplication constructor retains exponentially many syntactic residuals despite primitive normalization, operand swap and endpoint-product pruning"
source: "Direct residual-polynomial derivation in C118; Bryant (1991), section 5, for the distinct unrestricted multiplication-bit OBDD object"
claims: [C118]
notes: [FC]
todo: [72]
---
# TX39 — Constructing a fixed-product relation

## Dictionary

Conditioning factor bits gives a bilinear residual equation with inherited
integer intervals. Its coefficient record is a candidate contraction state.
C118 owns the recurrence, specified merging/pruning rule and scoped obstruction.

## Hypotheses

The chosen constructor conditions simultaneous low operand bits, identifies
only primitive coefficient tuples and their operand swaps, and prunes only
when the endpoint products exclude N. C118 states its surviving-prefix range.
Semantic equality of zero-sets is a stronger, unimplemented operation.

[Bryant](https://www.cs.cmu.edu/~bryant/pubdir/ieeetc91.pdf), section 5,
Theorem 4, bounds OBDDs for the middle bit of unrestricted multiplication.
Fixing every output bit to N changes that object, so its lower bound is not
transferred to the final fixed-product relation. Its hypothesis mismatch is
an escape worth investigating, not a constructive factoring result.

## Consequence for the goal

Stop scaling this particular frontier constructor. Stronger zero-set
identities, other elimination orders and more informative pruning remain open.
The final factor relation itself can have a tiny decision diagram; C118 does
not contradict that and does not prove a general factoring lower bound.
