---
id: 71
state: done
title: "Structured Grover: extend coherent selected probabilities to exact output sampling"
outcome: "Completed (2026-09-21): exact coherent sampler, rare128 validation and bounded rolling-versus-cached comparison in C117"
claims: [C114, C117]
---
# Structured Grover output sampling

## Completed

C117 owns the derived exact sampler, fixed-iteration scaling and measured
comparison through128 qubits. Note GS records independent design/result
review and execution; TX38 records the known MPS-sampling transfer. The
rare128 branches were checked exactly before resources. No full oracle-gate
peak or amplification-length run was part of this contract.

## Original proposed contract

The user seeks lower simulator memory with moderate extra runtime. TX35 and
C114 give a restricted query implementation for the explicit predicate
x0 AND no adjacent 11, with uniform initial state, Grover iterations and a
final H layer. Selected probabilities alone do not implement a sampler.

Derive prefix marginals using paired weighted-automaton contractions, retaining
coherent cross terms, and compare with direct sampling from the corresponding
small-bond MPS. Count construction, normalization, integer/rational bit growth,
RNG and retained outputs. Do not treat this known structured baseline as a
new speedup, or supply marked solutions/counts as free input.

Before execution, derive the method and build its tiny full-distribution
reference, then design the comparison and obtain independent design review.
The design must include omitted oracle, dephasing and broken local-clause controls with checked
preconditions; zero-probability prefixes; and an exact sampling contract
including the random-bit procedure. Only then execute the bounded memory/time
comparison. A general oracle, nonuniform input, noise, or native reversible
oracle peak is a different task. No experiment is scheduled by this record.
