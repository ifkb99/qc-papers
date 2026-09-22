---
id: TX21
field: "Fault tolerance without fresh qubits: non-unital noise and refrigeration"
status: escaped
effect: none
one_line: "Non-unital noise can refrigerate designed circuits; build_shor is not one, so its output still loses its dependence on (N, a)"
source: "Ben-Or-Gottesman-Hassidim 2013 (abstract); SYGY App. G pp.26-27 (slate 3a S5)"
claims: [C104]
notes: [ST]
todo: []
---
# TX21 — Fault tolerance without fresh qubits: non-unital noise and refrigeration

## Dictionary

Fresh-qubit supply by amplitude damping ↔ fault tolerance without resets.

## Hypotheses

Circuits designed to refrigerate; build_shor is not. Exponent coherence is damped by
sqrt(1-λ) per location (S5, sampling step conjectured).

## Consequence for the goal

Escaped obstruction: under amplitude damping the output becomes independent of (N,
a), as under TX20.
