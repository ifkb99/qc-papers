---
id: TX23
field: "Tensor-network contraction and treewidth (Markov-Shi)"
status: open
effect: lower-bound
one_line: "One cmult_mod gate graph has an m × 5n grid minor, so every contraction order of the native network costs 2^Theta(n)"
source: "Markov-Shi (abstract only; line-graph statement to read); grid minor derived in slate 3b K3 by reading toffoli_arith.py"
claims: [C45, C104, C106]
notes: [ST]
todo: []
---
# TX23 — Tensor-network contraction and treewidth (Markov-Shi)

## Dictionary

Native gate network ↔ graph; contraction cost ↔ exp(treewidth).

## Hypotheses

Grid minor: branch sets {MAJ_i, UMA_i} per adder instance, vertical edges along t_i,
horizontal along b_i (derived, not checked mechanically or refereed). Upper bound
q_w + O(1) by time slicing.

## Consequence for the goal

If confirmed, graph-uniform contraction of modexp is 2^Theta(n), so C45/C104's
2^(q_w+1) is optimal up to the constant in the exponent among such methods; only
value-level structure goes below. Treewidth is not a scalar lower bound. Cheapest
step: mechanical minor and min-fill at N = 7, 11.
