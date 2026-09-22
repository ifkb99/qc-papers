---
id: TX37
field: "Invariant subspaces and ancilla-aware simulation"
status: imported
effect: classification
one_line: "Clean arithmetic blocks restrict to logical modular permutations; certification and strong logical baselines determine whether this saves memory"
source: "toffoli_arith.py u_a; direct clean-domain certification and instrument derivation in C115"
claims: [C115]
notes: [VR, MC]
todo: [66]
---
# TX37 — Restrict complete arithmetic blocks to the clean subspace

## Dictionary

A certified complete arithmetic block acts on valid work labels as a modular
permutation while restoring scratch. Linear extension then acts on coherent
work amplitudes without explicitly storing scratch. C115 records the
restricted inverse-QFT sampling instrument and its bounded resource outcome.

## Hypotheses

The implementation certifies both control values on every x<N clean basis
input using the actual u_a gate list. Every allowed multiplier preserves this
subspace. Scratch cleanliness at block boundaries does not authorize discarding
live scratch inside a block, or arbitrary interleaved work rotations: Rx can
move amplitude to x>=N. TODO 64 is therefore not solved by this restriction.
This is a direct construction check, not an invocation of C105 for general
coherent channels or a new ancilla-removal theorem.

## Consequence for the goal

It avoids explicit scratch-state storage under the stated contract, which
existing logical sparse and semiclassical baselines already exploit. The
charged pilot in C115 loses to those strong baselines in total callable
runtime and total traced peak on the tested fixture; a two-buffer payload is
not total allocated memory. The restricted implementation is retained as an
experiment, with no default-backend change or asymptotic memory claim.
