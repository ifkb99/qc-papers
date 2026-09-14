---
code: SE
date: 2026-09-10
title: "Same physical computation with different Walsh cost; same Walsh cost with different physical computation"
outcome: established
claims: [C50]
todo: [15d]
---
# SE — Scratch-space operational equivalence

Fourth investigation in the requested sequence. The controlled construction
and its exact logical-code argument are in C50. The experiment is
`experiments/experiment_scratch_equivalence.py` and passed **36/36 checks**.

For each of N=7,a=6 and N=7,a=3, hold exponent width 4, total width 17,
layout and observable fixed. Vary only the number k=0..4 of appended
CCX(b_j,t_j,x_0) guards. The logical subspace contains every exponent and
every x<N with clean scratch: 112 basis inputs per variant. All variants
agree exactly on their full output images and restricted observables.

| Base | k=0 | k=1 | k=2 | k=3 | k=4 |
|---|---:|---:|---:|---:|---:|
| a=6: full Walsh support | 15549 | 15994 | 16030 | 16036 | 31982 |
| a=3: full Walsh support | 64353 | 64471 | 64454 | 64422 | 129152 |

For legitimate order-finding initialization x=1, the exponent-only sign
function is unchanged: Walsh support 1 for base 6 and 8 for base 3. The
complete post-inverse-QFT output distributions agree within 1.50e-15.
Gate-level evolution of a nontrivial coherent code state agrees up to one
global phase within 4.14e-16. The code identity is proved by the zero guard
controls; that proof, not a single random-state check, establishes the
coherent logical-isometry equivalence.

Arithmetic outputs are prepared using the independently gated circuit's
permutation replay; the guards and inverse QFT are evaluated with the existing
state-vector simulator. This is not an independent dense re-simulation of
every arithmetic gate for every variant. The original Fourier/Toffoli layout
confound is absent: all variants use exactly the same physical register.

## The useful converse

The deliberately invalid CNOT(e_0,x_0) changes legitimate arithmetic, but
exactly translates the Walsh coefficient vector. Full support stays 15549
or 64353 respectively. At a=3, the maximum difference in actual order-finding
output probability is **0.1015625**, despite the identical support count.
At a=6 those output statistics happen to remain unchanged (within numerical
error), though the logical arithmetic is wrong. The control was explicitly
predicted for the odd-order example, not for every changed computation.

Thus neither equal full-support counts nor differences in those counts certify
logical equivalence. This is consistent with Paper A's exact representation
identity and helps identify the task to which that identity applies.

```bash
OPENBLAS_NUM_THREADS=1 LAB_GPU=1 uv run python -m experiments.experiment_scratch_equivalence
```

Complete rows and provenance are in ignored `out/scratch_equivalence.json`.
The pilot is complete; discovering a general cheap clean-subspace quotient
or optimizing circuit extensions is a separate research problem.
