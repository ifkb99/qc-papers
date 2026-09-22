---
id: TX13
field: "Counting complexity: Holant and #CSP dichotomies, holographic algorithms"
status: obstruction
effect: lower-bound
one_line: "The Toffoli-arithmetic gate language is #P-hard in every Boolean dichotomy checked, planar included; only gadget-level, instance-level or number-theoretic routes remain"
source: "Cai-Fu arXiv:1603.07046 Thm 1.1 p.2, Thm 2.31 p.17; Guo-Williams arXiv:1212.2284 Thm 1.1; Cai-Fu-Guo-Williams arXiv:1505.02993 Thm 8.1 p.60 (read in slate 3b)"
claims: [C86, C87]
notes: [MG, ST]
todo: [50]
---
# TX13 — Counting complexity: Holant and #CSP dichotomies, holographic algorithms

## Dictionary

Walsh coefficient or amplitude ↔ Holant value on the circuit network; with fan-out
realisable, #CSP({T4, XOR3, ≠, [1,±1], pins}); path-sum form #CSP({CCZ, CZ, Z, H}).

## Hypotheses

T4 and CCZ lie outside the affine, product and Hadamard-matchgate classes; pins lie
outside the latter; a CNOT gadget realises the crossover, so planar and general
coincide (slate 3b §1.2, run Rc5b7bf01443648f7). Hardness is worst-case over the
gate language under Turing reductions, conditional on FP ≠ #P.

## Consequence for the goal

Closes every gate-uniform holographic method (one local basis change on all wires,
then affine, product or FKT evaluation). Says nothing about the fixed ToffoliModExp
family, gadget-level signatures (C86 escapes this way) or value-level algebra. Pins
outside the matchgate class explain why C86's route needs the full-space contract.
The earlier wording ('closes every holographic route') was wrong.
