---
id: TX6
field: "Matchgates, Pfaffians and the FKT algorithm"
status: imported
effect: removes-exponential-in-subfamily
one_line: "Planar-coupled carry chains reduce to Pfaffian / dimer counting in polynomial time"
source: "Valiant's matchgates; FKT, as used in C86 to C88"
claims: [C86, C87, C88]
notes: [MG]
todo: [48, 49, 50]
---
# TX6 — Matchgates, Pfaffians and the FKT algorithm

## Dictionary

Carry-chain signed queries ↔ weighted perfect matchings of a planar graph.

## Hypotheses

Needs planarity of the coupling; actual modexp coupling is not known to be planar
(TODO 50).

## Consequence for the goal

The one place the ledger removes an exponential for arithmetic queries, in a
restricted family. Its natural generalization is TX13.
