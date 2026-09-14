---
code: TR
date: 2026-09-09
title: "Walsh support versus tensor rank: tiny ideal functions, substantially larger scratch-space tensors"
outcome: mixed
claims: [C48]
todo: [15b]
---
# TR — Tensor rank and finite-state memory

Second investigation in the requested sequence. C48 contains the exact local-
basis and residue-factorization statements. The diagnostic is
`experiments/experiment_tensor_memory.py`; reusable tools are in
`lab/tensors.py`. These are dense measurements, not a compressed propagator.

The saved run passed **72/72 checks** over 29 rows. All reported ranks are
numerical over R, not GF(2), and remain unchanged across relative tolerances
10^-8, 10^-10, 10^-12 (with corresponding absolute floors). Representative
normalized cut singular-value vectors agree before/after Walsh transformation
within 2e-13. Every natural-order cut is measured for each row; only selected
circuit suffixes are sampled, not every gate.

## The promising half

| Ideal scalar function | Exponent widths | Walsh support | Maximum cut rank |
|---|---|---|---|
| N=7,a=3, r=6, LSB | 4,8,12,16 | 8,128,2048,32768 | 3 throughout |
| N=13,a=4, r=6, LSB | 4,8,12,16 | 1 throughout | 1 throughout |
| N=15,a=7, r=4, LSB | 4,8,12,16 | 4 throughout | 2 throughout |

The residue automaton reproduces the scalar function on the tested inputs.
The first row cleanly separates dense coefficients from tensor rank. The second
recovers the earlier scalar-period counterexample rather than attributing all
cost to the multiplicative order. A fixed seeded random 12-bit sign function
has support 4096 and maximum cut rank 64: the small-rank control fails as it must.

## The cautionary half

For the full N=7 arithmetic function (including all invalid scratch inputs):

| Base | t=2 support/rank | t=3 support/rank | t=4 support/rank |
|---|---|---|---|
| 6, r=2 | 15549 / 128 | 15549 / 128 | 15549 / 128 |
| 3, r=6 | 15539 / 128 | 30712 / 153 | 64353 / 217 |

The rank-3 ideal function does not make its full scratch-space pullback rank 3.
In the t=3 suffix scan, 645 reversed logical gates give support/rank 258/17 for
base 6 and 730/15 for base 3; one completed final block gives 2606/65 and
3142/64 respectively. Early sampled suffixes leave the observable unchanged
(rank 1), so those trivial rows are not evidence for general compressibility.
The full-space ranks already rule out the most optimistic tiny-memory guess
in this ordering. They do not prove an asymptotic lower bound for every ordering.

## Reproduction and interpretation

```bash
OPENBLAS_NUM_THREADS=1 LAB_GPU=1 uv run python -m experiments.experiment_tensor_memory
```

Raw profiles, singular-value gaps and provenance are in ignored
`out/tensor_memory.json`. CPU fallback uses the existing Walsh replay. Each
standard arithmetic builder is checked on all legitimate basis inputs against
modular exponentiation and scratch cleanup; the existing science suites supply
the independent state-vector gate for these unchanged builders. An initial run
was interrupted because it redundantly reran the full state-vector verifier
for every exponent value at every width; no prediction was changed.

The finite-state/tensor connection is established prior art:
[Kiefer, weighted-automaton minimization](https://arxiv.org/abs/2009.01217) and
[Li, Precup and Rabusseau, automata and tensor networks](https://arxiv.org/abs/2010.10029).
The useful project-specific outcome is the controlled ideal/full-space contrast,
not a claim that tensorizing a Walsh transform is new. Compression algorithms,
ordering searches and minimal weighted-automaton discovery remain open beyond
this first-pass diagnostic.
