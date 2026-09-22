---
id: TX31
field: "Repetends of y/N and cofactors of 2^M ± 1 (Cunningham numbers)"
status: open
effect: classification
one_line: "Escapes need structure: multiplication by 2 rotates the repetend of y/N, so a in ±<2> with N a cofactor of 2^M ± 1 gives small diagrams"
source: "Deriver slate 4 (Se2325ec074c0487a) S1: digit identity proved; cofactor escape sketched; quadratic sizes measured"
claims: [C107]
notes: [EF]
todo: []
---
# TX31 — Repetends of y/N and cofactors of 2^M ± 1 (Cunningham numbers)

## Dictionary

parity(2^j y mod N) = bit_{-j mod m}(D_N y), m = ord_N(2), D_N = (2^m - 1)/N ↔ a
rotation after the change of coordinates y -> D_N y.

## Hypotheses

Digit identity proved (unrefereed); 'cofactor height' h(N) as the deciding quantity
is a conjecture; run count of D_N refuted as the quantity.

## Consequence for the goal

Characterizes the escape side; every escape found has r = O(n), so it respects TX15
trivially.
