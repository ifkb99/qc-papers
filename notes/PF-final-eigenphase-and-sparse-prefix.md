---
code: PF
date: 2026-09-10
title: "Final-eigenphase conditioning survives a genuine work mixer; sparse rejection removes the orbit table under an explicit input contract"
outcome: mixed
claims: [C54]
todo: [14, 18, 19]
---
# PF — Sample the eigenphase after the defect

The user requested continued research, with lower-cost agents doing initial
tests and potential breakthroughs explicitly flagged. Two `gpt-5.6-luna`
agents independently implemented the physical mixing experiment and the
conditional-prefix probe. The main agent derived and implemented the reusable
sampler, inspected their code/results, corrected test hygiene, and reran the
experiments. This was concrete progress on the active research goal, not a
claim that the open-ended goal is finished.

The qsim-research workflow again changed the comparison: failure of initial
spectral dephasing was only a control. The stronger baseline conditions on a
FINAL eigenphase after the defect and retains a short coherent exponent
prefix. C54 owns that derivation and its sparse-rejection bound.

## The physical experiment

`experiment_mixing_defect.py` fixes the same tiny arithmetic instance and
insertion point as the preceding experiment, replacing its diagonal gate by
the clean-orbit-preserving mixer between work labels |1> and |5>. The gate
is four commuting Pauli rotations built by the existing Circuit engine;
its generator is specified in the experiment. No new circuit propagator was
written. Only its angle varies in the main series.

The full physical unitary is checked on all work states and spectators before
the sweep. References are the existing Fourier-compiled arithmetic with
`statevec.run`, the time-ordered finite spectral effects, and a clean-basis
joint amplitude table followed by FFT. The physical arithmetic order is never
reversed. The independent agent's completed run passes **27/27 checks**:

- Full-unitary error is below 2e-16 and orbit leakage is zero.
- Spectral/circuit output errors are below 2.4e-14; clean-FFT errors below 2e-16.
- The old ideal scalar model fails by up to about 0.0667 per output, and initial
  physical-work spectral dephasing by about 0.0483.
- Moving the defect to the end fails by about 0.0357. Reusing the previous
  diagonal phase offset fails by about 0.122.
- Extended spectral precision changes the checked row by below 2e-16; an
  extended-state replay with the existing gate constants agrees below 3e-14.

This is a genuine work-label mixer, but those failures are not hardness
evidence. The weighted final-eigenphase baseline succeeds.

## Independent prefix and rejection verification

`experiment_prefix_probe.py` passes **18/18 checks**. It includes seeded
complex orbit unitaries at two periods and two prefix sizes, including an
early prefix longer than the orbit. Its independent route constructs the
conditional CONTROL state and applies the existing coherent inverse-QFT
engine. This differs from both the spectral-effects calculation and the
streamed-column implementation in `lab/prefix.py`.

The weighted mixture agrees with the spectral outputs within 8e-16; the
late/early sampling factorization agrees with coherent inverse QFT within
1.6e-15. The physical mixer gives nonuniform final weights, about
0.1306–0.2028. Uniform weighting fails by about 0.022, and erasing the early
amplitude phases fails by about 0.212. These controls prevent mistaking the
prefix for an ordinary classical probability table.

The agent also checks the sparse rejection envelope on generic COMPLEX
normalized columns, not only the arithmetic mixer. For supplied bounds
D=1,2,3 the enumerated mean acceptance is 1, 1/2, 1/3 up to floating-point
roundoff, and no envelope violation occurs. This is a Parseval/Cauchy-Schwarz
identity checked by enumeration, not an inference from an acceptance histogram.

Main-agent review removed an unnecessary tolerance cutoff in the probe oracle
and a tautological “orientation” diagnostic that recomputed its own expression.
The genuinely independent circuit/effect comparisons remain. All predictions
and meaningful controls still pass after those changes.

## Bounded scaling, with its assumptions visible

`experiment_prefix_scaling.py` passes **39/39 checks**. It holds the sparse
mixer and four-amplitude early prefix fixed. One series varies exponent width
through 2–8, 16, 32 and 63 at fixed period. Another uses supplied abstract
periods through 1,000,000,007 at fixed width. These large-period cases are
explicitly INDEXED CYCLIC-ORBIT ORACLES, not compiled modular circuits or a
large physical-qubit benchmark. No order-finding speedup is inferred.

Each row samples sixteen fixed-seed paths, compares their reported joint
(eigenphase, output) probability with a forced path, and checks it independently
using an integer-reduced finite geometric tail plus a direct four-point DFT.
Full small marginals agree with spectral effects within 3.1e-15. The largest
sampled independent absolute error is below 2e-16 and positive-reference
relative error below 6.1e-15. Every strictly positive reference is included
in the relative-error report; there is no probability threshold hiding rare
cases. These are the sampled paths, not a uniform precision certificate for
all possible outputs. None of these samples has zero reference probability.

Each sampled phase row has four complex amplitudes, 64 bytes of payload.
This is NOT peak process memory: rows, FFT/feedback arrays, sparse columns,
Python objects and caller-owned data coexist. The largest dense array in the
physical unitary check is 262,144 bytes; a pre-allocation guard enforces its
16 MiB limit. Full reference enumeration in the scaling probe is separately
guarded to the small period and width range. Setup column calls and every
sample's proposal/column counts are saved. The oracle does not own an orbit
table. The table-free result relies on the exact known-index input contract
in C54; arbitrary computational labels cannot silently substitute for indices.

### A caught allocation error and an imperfect artifact trail

The initial scaling probe mistakenly attempted a dense identity matrix in its
largest abstract-period row. NumPy rejected the shape during size preflight
with “array is too big”; the agent reports no evidence of material RAM
allocation. The corrected code never constructs dense reference data for
abstract rows and now has explicit guards before reference allocations.

The first failure log was overwritten by the agent's rerun. It cannot be
claimed preserved: `out/prefix_scaling_initial_failure.md` records the reported
traceback and this limitation, rather than inventing a replacement raw log.
This failure was in the experiment's reference setup, not the prefix sampler,
and is not scientific evidence for the scaling result. Main-agent inspection
and a fresh safe rerun confirmed the final guarded code and successful report.

## Reproduction and review

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_mixing_defect
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_prefix_probe
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_prefix_scaling
```

Matching JSON/log files are under ignored `out/`, including all declarations,
checks, failures-as-controls, measured probabilities and resource accounting.
The completed runs use Python 3.12 and NumPy 2.4.6. The main agent reran all
three experiments after review. The core, lab and claims science suites pass;
their logs are `out/mixing_test_core.log`, `out/prefix_test_lab.log`, and
`out/prefix_test_claims.log`. The other six suites were not rerun for this
additive helper; previous consolidation results remain separate evidence.

A second lower-cost agent audited zero/full prefixes, repeated orbit indices,
invalid columns, impossible latent phases and retry-cap failure. The saved
audit is `out/prefix_code_audit.md`. It found no numerical discrepancy in its
bounded cases. Main fixed its minor `None`-label error-type finding and added
a regression; Python bool-as-int acceptance remains consistent with surrounding
integer APIs. The audit is a dated report, not an unchanged-source certification.
The ten-check documentation gate passes. Papers and abstract workshops were
not modified in this follow-up, and no commits were made.

## Impact and what remains open

This was flagged during derivation as a candidate compression result. Initial
agent tests and main-agent checks now support the mathematical mechanism and
implementation. It is **not established as a research breakthrough**: the
construction uses standard conditioning and rejection sampling, and does not
remove order discovery, arbitrary orbit-index access, or an exponentially
large early prefix. C54 gives the precise positive result and limitations.

Primary-source review included Van den Nest §4 on computationally tractable
states/sparse operators and Schwarz–Van den Nest Theorem 1 on approximately
sparse OUTPUT distributions. See C54 for references. That output-sparsity
promise is distinct from this sampler's sparse-defect-column promise; neither
paper is cited as proving our specific implementation or its novelty.

TODO 18 closes this bounded experiment. TODO 19 asks whether a localized
finite-rank defect admits a representation that also removes the large-prefix
FFT. Cheap evaluation of a probability is not yet a method for sampling the
whole distribution: that distinction remains a required control on the next
direction. TODO 14 remains the broader frontier and the thread goal stays active.
