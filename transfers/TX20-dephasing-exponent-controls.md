---
id: TX20
field: "Noisy computation without fresh qubits: dephasing of control registers"
status: open
effect: classification
one_line: "Dephasing on exponent qubits commutes to the start, so noisy build_shor's output is near uniform once each exponent qubit has one dephasing event"
source: "Aharonov-Ben-Or-Impagliazzo-Nisan 1996 Thm 3, Cor 2 p.4 (entropy idea); the lemma is derived in slate 3a S2"
claims: [C104]
notes: [ST]
todo: []
---
# TX20 — Noisy computation without fresh qubits: dephasing of control registers

## Dictionary

Loss of an exponent control's coherence ↔ ABIN96's entropy increase; exponent qubits
only ever act as controls (C104 Roles), verified at rotation level (run
R7ace5f08b1564841).

## Hypotheses

Pauli noise with a dephasing component; every gate touching an exponent qubit block-
diagonal on it; output = exponent register only. Fails for compilations where an
exponent qubit is a target, for noise without dephasing (TX21), and for contracts
that also measure the work register. Derived, coordinator re-read, not refereed;
likely folklore (decoherence in Shor, Miquel-Paz-Zurek, not yet searched).

## Consequence for the goal

TV(noisy, uniform) ≤ 1 - Π_j(1 - (1-λ)^(g_j)), 13n ≤ g_j ≤ 12n²+n: trivial for λ ≳
ln(n)/n worst case, ≈ ln(n)/n² typical. With TX15's low-noise amendment, only a
vanishing window remains. A noise threshold, not an exponential removed.
