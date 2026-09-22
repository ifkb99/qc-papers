---
id: TX19
field: "Quantum fault tolerance: threshold theorems"
status: obstruction
effect: barrier
one_line: "Any fixed-input noisy simulation valid for all circuits below threshold would simulate fault-tolerant Shor"
source: "Aharonov-Ben-Or quant-ph/9906129; Aliferis-Gottesman-Preskill QIC 6, 97 (2006); SYGY p.1 and fn 69 (abstract-level reading in slate 3a)"
claims: [C52]
notes: [ST]
todo: []
---
# TX19 — Quantum fault tolerance: threshold theorems

## Dictionary

A general fixed-input noisy simulator ↔ a simulator of an encoded, fault-tolerant
Shor circuit ↔ TX15.

## Hypotheses

Needs fresh qubits and an encoded circuit. build_shor has neither, so the row does
not bind it; it binds every general fixed-input transfer.

## Consequence for the goal

Prunes any noise-based route that claims generality beyond unencoded circuits. For
build_shor itself, TX20 applies instead.
