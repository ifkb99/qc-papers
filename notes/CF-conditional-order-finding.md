---
code: CF
date: 2026-09-10
title: "Conditional measurement reproduces actual order-finding outputs; commutation and exponential work cost remain essential"
outcome: established
claims: [C49]
todo: [15c, 14]
---
# CF — Conditional order-finding, beyond the averaged observable

Third investigation in the requested sequence. C49 states the instrument,
its commutation condition and resource scope. The helpers are in
`lab/semiclassical.py`; `experiments/experiment_conditional_order_finding.py`
extracts the branch permutations from the existing Toffoli arithmetic builder.

The completed run passed **43/43 checks**, including four negative controls.
For (N,a)=(7,6),(7,3),(15,7), all output probabilities at exponent widths
3,4,6,8 agree with an independent FFT of the ideal order-finding state to
within 2.23e-16. At width 3, the full existing Fourier-arithmetic quantum
circuit (`ModExp.build_shor()` through `statevec.run`) independently agrees
within 5.42e-14. This compares the same valid-input physical task across
different arithmetic layouts, not their full-space Walsh spectra.

At widths 4,8,16,32, sixteen paths per instance/width were sampled without
enumerating the other outputs. Every sampled path probability agrees with
the independent geometric-series formula to within 1.39e-16 in the saved run.
The period is used only by this reference and the diagnostic metadata, never
by the sampler or branch-map builder. The reference reduces phase arguments
with integer arithmetic before evaluating trigonometric functions.

| Instance | Stored complex amplitudes | Maximum occupied amplitudes at t=32 |
|---|---:|---:|
| N=7,a=6 | 8192 | 2 |
| N=7,a=3 | 8192 | 6 |
| N=15,a=7 | 65536 | 4 |

The allocated work-vector size is unchanged as exponent width increases.
Occupied support is counted at tolerance 1e-10 for diagnostics only; the
algorithm does not truncate amplitudes. These counts are not peak allocated
bytes: temporary child vectors, permutation maps and classical histories also
consume memory. The dense transition-table construction is exponential, and
no large-modulus scalability claim follows from these small moduli.

## Controls and the important audit finding

Omitting feedback fails on the odd-order case. Reversing output-bit order
fails. Replacing the measured instrument with fair bits fails. All are checked
against full distributions, not one expectation.

The fourth control was added during implementation review: measuring controls
in inverse-QFT order and moving their arithmetic interactions into that order
requires a commutation argument. It is available on the valid modular-
multiplication subspace, not for arbitrary full-scratch maps. Two differently
controlled, noncommuting CNOT work gates give a maximum output-probability
error of **0.125** under the naive reordering. This is now documented in both
the helper and C49; the prototype is not advertised as a general noncommuting
circuit simulator.

## Reproduction and prior art

```bash
OPENBLAS_NUM_THREADS=1 LAB_GPU=1 uv run python -m experiments.experiment_conditional_order_finding
```

The complete distributions, sampled paths, profile rows and provenance are in
ignored `out/conditional_order_finding.json`.

The terminating semiclassical QFT is established by
[Griffiths and Niu](https://arxiv.org/abs/quant-ph/9511007); its simulation
implications, including appropriate structured input states, are discussed by
[Browne](https://arxiv.org/abs/quant-ph/0612021). The contribution here is a
matched, tested reference in this repository, not the general method.

TODO 15c's pilot is complete. TODO 14 remains open for a useful compressed
representation/cost characterization through the inverse QFT. The immediate
engineering opportunity is to avoid dense full-work transition tables while
retaining the exact conditional interference and auditing setup cost.
