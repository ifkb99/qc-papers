---
code: RO
date: 2026-09-10
title: "Reachable compiled sampling removes scratch tables; charged classical spectral sampling is stronger"
outcome: mixed
claims: [C51, C52]
todo: [14, 16]
---
# RO — Reachable orbits and the latent-eigenphase baseline

The user requested consolidation into the two paper manuscripts, followed by
further research. PAPER_A.md now integrates the operational-equivalence and
tensor-memory findings; PAPER_B.md integrates contraction, balanced controls,
conditional sampling and this follow-up. Claim files remain authoritative;
abstract workshops were not rewritten as duplicate manuscripts.

## Derivation and experimental design

C51 gives the descending-power reachable-set bound. It was derived before
measurement. The implementation exposes the SAME existing classical replay
kernel on selected inputs (`walsh.classical_images`), then memoizes compiled
branch actions in `lab/reachable.py`. No full work-space permutation table is
built; neither an order nor a supplied orbit goes into the compiled sampler.

The correctness gate replays every clean work input x<N for each distinct
compiled multiplier, on both control branches, checking arithmetic and scratch
cleanup. It uses separate wrappers/caches so that measured sampler setup and
cold replay do not benefit from a pre-enumerated clean code. Every queried
cached input and output is checked clean at block boundaries. The ordinary
noncommuting-gate counterexample remains a required negative control.

Two one-parameter series are kept separate: N=21,a=2 while varying exponent
width, then a=2,t=32 while varying N across the fixed list below. Changes in
order in the second series are consequences of the changing modulus, not a
claim to have held the arithmetic instance fixed. Fixed RNG seeds couple the
compiled, high-level modular, and orbit-vector samplers for exact path checks.

## Results

The completed experiment passed **141/141 checks**. Complete small output
distributions agree with an independent FFT, and the smallest is checked
against the existing Fourier-compiled state-vector simulator. Wider sampled
paths agree with geometric-series marginal probabilities and two classical
work-state baselines. Forced rare outcomes reach probabilities around 10^-19;
the largest positive-reference relative error is **1.78e-14**. No amplitude
tolerance truncation is used; `occupied_1e10` is only a diagnostic count.

At fixed N=21,a=2, stored work support rises from 4 at t=2 to 6 at t=3 and
stays at 6 through the checked widths up to t=32. At fixed a=2,t=32:

| N | r | Compiled work bits | Peak stored amplitudes | Cached transitions after sampled paths |
|---|---:|---:|---:|---:|
| 15 | 4 | 16 | 4 | 8 |
| 21 | 6 | 19 | 6 | 18 |
| 35 | 12 | 22 | 12 | 30 |
| 77 | 30 | 25 | 30 | 150 |
| 143 | 60 | 28 | 60 | 210 |
| 221 | 24 | 28 | 24 | 54 |
| 323 | 72 | 31 | 72 | 234 |
| 437 | 198 | 31 | 198 | 5202 |
| 667 | 308 | 34 | 308 | 4258 |

For the last row a scratch-dense vector would require 2^34 complex amplitudes,
or 256 GiB of complex128 payload alone. **That vector was not allocated.**
The sparse result stores 308 work amplitudes plus labels, temporary vectors,
gate lists, caches and history. This is NOT a measured process-memory reduction
ratio and not a generic "66-qubit simulation" benchmark: exponent width is 32,
work width 34, but the arithmetic structure determines which labels are
evaluated. The compiled helper's int64 labels currently cap the reused
work-plus-control circuit at 63 qubits.

## A useful negative result, then a stronger comparator

The first run passed 135 checks and showed that ordinary orbit-vector sampling
was faster than compiled sparse sampling. This is useful: removing the scratch
overhead has not established an advantage over an obvious classical method.

That prompted the C52 spectral baseline, derived and declared before a second
run. All conditional effects are polynomials in one modular-multiplication
unitary. Initial |1> gives a uniform distribution over its orbit eigenphases,
so one can sample a latent k/r and generate output bits by scalar Bernoulli
updates. The underlying argument is explicitly prior art in
[Cleve et al., §6](https://arxiv.org/abs/quant-ph/9708016), equations 6.2–6.4
and the measurement-commutation argument that follows. We verified the small
mixtures by summing conditional-on-k probabilities, not by mistaking one such
probability for the marginal.

The spectral comparator is charged for discovering r by repeated multiplication
until return to 1. That takes O(r) arithmetic steps and no orbit table. It then
needs no work vector to sample. The added diagnostic also verifies coherent
exponent/work Schmidt rank min(r,2^t), including a scalar-rank-1 example whose
joint rank is 6. These ranks are not output-sampling memory lower bounds.

Illustrative **single-run** timings at N=667,a=2,t=32 (Python 3.12, CPU):

| Work | Seconds |
|---|---:|
| Compile distinct arithmetic circuits | 0.959 |
| First compiled sparse sample, initially empty caches | 0.863 |
| Eight compiled sparse samples, after the first sample | 0.0157 |
| Enumerate orbit and construct orbit transition arrays | 0.000497 |
| Eight orbit-vector samples | 0.00215 |
| Discover order without storing orbit | 0.00000963 |
| Eight latent-eigenphase samples | 0.000503 |

Warm sampling may still populate missing cache entries; it is not advertised
as fully cached. Compilation time includes the existing builder's emitted
rotation representation as well as its logical trace. Times are illustrative,
not repeated competitive benchmarks; microsecond setup timings especially do
not support precise ratios. High-level modular and spectral baselines do not
simulate arbitrary dirty-scratch arithmetic, but they do solve the same ideal
clean-input output-sampling task under comparison. No MPS superiority is claimed.

## Reproduction and provenance

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_reachable_order_finding
```

Report: ignored `out/reachable_order_finding.json`, with declarations, checks,
profiles, cache/setup costs, rare outputs, Python/NumPy versions and timings.
Completed run: Python 3.12.10, NumPy 2.4.6, CPU, about 7.81 seconds overall.
The initial 135-check run succeeded under Python 3.14 with CUDA for the small
dense reference. An expanded 3.14 run exited 139 during the modulus sweep,
without a Python traceback; its partial log is preserved as
`out/reachable_order_finding_py314_failed.log`. The known interpreter instability
motivated the 3.12 rerun; this particular crash's cause was not diagnosed.
All predictions were retained. Do not compare the cross-interpreter compilation
times as an algorithmic speedup.

The fresh nine-suite science gate and ten-check documentation gate are recorded
in the consolidation handoff. This bounded reachable-state follow-up is done.
TODO 14 now asks where a controlled noncommuting perturbation makes the
spectral reduction inadequate; bigger ideal-orbit sweeps alone are postponed.
