---
code: SD
date: 2026-09-10
title: "A noncommuting work defect breaks physical eigenphase dephasing, but phase transfer restores scalar sampling"
outcome: mixed
claims: [C53]
todo: [14, 17]
---
# SD — Detectable coherences are not a sampling-memory barrier

The qsim-research workflow led this follow-up to preserve the actual time
ordering, derive a discriminating test first, and add the strongest simple
comparator. That last step changed the interpretation before measurements:
the selected defect has an input-specific exponent-phase replacement. A failure
of the old physical-work mixture would therefore not establish a difficult
sampling problem. C53 gives the derivation and its scope.

## Design and independent checks

The main series fixes the modulus, base, exponent width, insertion point and
work gate from C53's example; only theta varies. Arithmetic powers are applied
in their original ascending order, with the defect after the first two blocks.
All exponent outputs are compared across three distinct calculations:

- finite spectral effects from the two time-ordered Fourier filters;
- the existing `statevec.run` engine on the existing Fourier-compiled arithmetic,
  inserting the one physical gate before the remaining arithmetic and inverse QFT;
- an independent computational-basis joint amplitude table followed by NumPy FFT.

Before the sweep, existing selected Toffoli replay validates every clean work
input on both branches of every needed multiplier, including scratch cleanup.
The spectral oracle is explicitly order/orbit-informed. The main circuit
reference is coherent throughout and is never replaced by descending-power
arithmetic. A second clean initial work label checks complex effect orientation,
not just the all-ones spectral density of |1>.

Controls include zero perturbation; a nontrivial commuting XXX rotation
(equal to a power of U on this orbit); an invalid move of the work rotation
to the end; and a separate earlier-insertion case where the gate is
noncommuting on the orbit but scalar on the actual reachable states. The
phase-adjusted scalar sampler is checked against the full coherent output
distribution, not merely against the original ideal distribution.

## Measurements

The completed experiment passes **57/57 checks**, with no harness warnings.
The largest full-circuit probability discrepancy is **2.35e-14**; the
independent FFT discrepancy is below **3e-16**. Selected rows of the main sweep:

| theta/pi | Physical-work dephasing TV error | Detectable eigenbasis pairs | Real coherence-response rank |
|---:|---:|---:|---:|
| 0 | below 2e-16 | 0 | 0 |
| 1/8 | 0.032836 | 12 | 11 |
| 1/4 | 0.068602 | 12 | 11 |
| 1/2 | 0.156990 | 12 | 11 |
| 1 | 0.251481 | 6 | 5 |

Negative angles are also tested. These pair/rank counts are finite numerical
diagnostics at absolute tolerance 1e-10, not a general rank theorem. At the
half-pi row, pairs (0,3), (1,4), (2,5) are undetected; the other pairs enter
the measurement. The eleven nonzero singular values are well separated from
the floating-point floor (smallest about 0.048, next below 4e-16).

The largest single-output error from discarding initial work-eigenphase
coherences is **0.0520833** (numerically 5/96), while leaving the original ideal sampler unchanged
gives a largest error of **0.15625**. Moving the half-pi defect to the end
produces a **0.078125** error. All required failure controls therefore detect
real differences. Yet the phase-adjusted scalar mixture recovers the coherent
distribution to the circuit's numerical precision. This is the principal
result: the original model fails, but this test case remains cheaply sampleable.

The separate earlier-insertion control changes the output by less than
3e-14 despite a nonzero orbit commutator. Noncommutation is therefore an
insufficient test even for whether the observed distribution changes.

An extended-precision spectral calculation changes the half-pi result by
less than 2e-16. Running the SAME coherent engine with an extended-precision
state agrees within 2e-14. Its pre-existing compiled gate constants remain
float64: this is not an arbitrary-precision certification. The math identities
are exact; all reported distributions are floating-point evaluations.

## Implementation and resources

`lab/spectral.py` provides only bounded dense finite-sum/effect diagnostics,
not another circuit propagator or a compressed sampler. Width and entry caps
are enforced before allocation, with Python-integer budgeting to avoid int64
overflow. `lab.semiclassical.eigenphase_path` gains optional product exponent
phase offsets. The default ideal path is unchanged and allocates no phase
vector; supplied offsets require O(t) storage in this implementation.

Each complex128 coherent reference state has 4,096 amplitudes / 65,536 bytes
of payload. The spectral effect array has 576 complex entries / 9,216 bytes.
These are individual array payloads, NOT process peaks or a memory-advantage
benchmark. Gate lists, temporary vectors, filters, SVD work and all-output
enumeration coexist. Orbit discovery, circuit compilation, compiled clean-input
validation and each reference/comparator timing are recorded separately in JSON.
No lookup cache is passed from validation into a claimed timed sampler.

The experiment is small: roughly a few CPU seconds, using Python 3.12.10 and
NumPy 2.4.6. It enumerates complete distributions solely for validation. The
phase-adjusted scalar *sampling* path is separately covered in `test_lab.py`
against coherent inverse-QFT circuits with arbitrary product phase offsets;
full-mixture enumeration timings are not advertised as single-sample timings.

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_spectral_defect
```

Raw report/log: ignored `out/spectral_defect.json` and `out/spectral_defect.log`.
The report retains all declared predictions, controls, per-output probabilities,
pair contributions, singular values, errors, payloads and setup timings.
Regression coverage is in the existing lab and claims suites; validation logs
are `out/spectral_test_lab.log` and `out/spectral_test_claims.log`.
The core science gate, affected lab and claims suites, and all ten documentation
checks pass. The other six science suites were not rerun for this bounded
helper/API extension; their prior consolidation runs are separate evidence.

## Interpretation and next step

This is a useful negative for the proposed hardness direction, and a sharper
measurement criterion. It is NOT a breakthrough in general classical quantum
simulation. The familiar distinction between coherence being present and a
measurement detecting it is explicit in Theurer et al., Definition 1,
Proposition 2 and Definition 3 of
[Quantifying Operations with an Application to Coherence](https://arxiv.org/abs/1806.07332).
We read those definitions/proposition in the primary paper; they are framing
and prior art, not evidence that our circuit-specific phase transfer is novel.

The bounded experiment is closed in TODO 17. TODO 14 carries the next question:
test a clean-orbit-preserving work rotation that genuinely mixes basis labels,
then compare with a rewritten-circuit/small-prefix baseline before building
any general tensor representation. The target should leave the early reachable
subset while remaining in the full orbit: a unitary confined to the early
subset can transfer its transpose to the equally weighted, perfectly correlated
early exponent register. This is a design observation, not another measured
experiment or a claim of hardness for the revised target.
Larger sweeps of this diagonal defect would
mostly reverify the identity already derived in C53. Neither manuscript nor
the abstract workshops was changed during this experiment.
