---
code: VP
date: 2026-09-11
title: "Verified exact-input prefix arithmetic and finite-bit sampling, with precision costs exposed"
outcome: mixed
claims: [C60, C61]
todo: [14, 24]
---
# VP — The oracle requirement is now implemented for specified exact gates

**Follow-up:** C62/§NG distinguish the observed doubling endpoint from the
smallest scanned sufficient precision and test an opt-in norm enclosure.
The execution numbers below remain historical; they are not minimum-precision
or worst-case-runtime results.

C61 owns the mathematical construction, input/trusted-arithmetic contract and
cost limits. This follows §PE: the existing rational planner now has an
exact-input ball-arithmetic oracle and an integer categorical kernel to use.
TODO 24 stays open for the uniform working-precision bound and certified
same-accuracy rejection comparison, not because these two components remain
unimplemented. Neither paper nor either abstract workshop was edited.

Main read the backend's primary documentation, implemented
`lab/verified_prefix.py`, factored the original finite component contraction
into a shared arithmetic-callback helper, and wrote the exact tiny transition
experiment/regressions. Three existing lower-cost agents investigated the
backend, independent full-r reference and integer kernel; main audited and
reran their tests. The qsim-research skill's exact-input, independent-control
and fact-ownership rules materially constrained the claims recorded here.

The old float sampler's public API, gate order and mathematical algorithm are
preserved. The shared-loop refactor can change last-bit floating rounding
(multiplication by a reciprocal square root replaces division); it is not a
bitwise-identical execution promise. No rounding certificate is added to that
sampler. No dependency was added to project metadata; verified experiments
explicitly request the optional pinned backend through uv.

## Backend and exact kernel probes

The backend probe passes 5/5 checks in
`out/ball_backend_20260911T033015739206Z.json`. Exact rational phase inputs,
binary midpoint/radius extraction and precision-context restoration work.
At the fixed cancellation target 10^-30, working precision 64 misses the
target and 128 reaches it. A 256-bit midpoint converted through Python float
leaves the original tiny ball, a decisive must-fail shortcut control.

The corrected integer-kernel report is
`out/dyadic_kernel_20260911T034344168055Z.json`, 6/6 checks. It includes nine
signed/tied rounding cases; four rational categorical fixtures; and a bounded
sweep of 1,680 integer-weight/bit-count cases totaling 10,416 words. Every word
is compared to an independent Fraction strict-CDF rule, not just bin counts.
All-zero fallback, exact-zero filtering, negative rejection and scripted
initial-integer rejection have explicit checks. Randomness quality itself is
assumed, not certified by those tests.

## Independent tiny amplitude diagnostic

Fix r=10,b=2,t=4, W1=Rx(pi/7), W3=Rz(pi/5), and reflections K2(q=0,pi/5),
K3(q=1,pi/5). Use the EXISTING direct full-r matrix product, not a second
generic propagator. At each arithmetic/background/reflection boundary and
each unique partial-QFT label, compare the dyadic oracle to that float
reference. Only requested coordinate accuracy varies.

The final report is
`out/verified_prefix_reference_20260911T034430917442Z.json`, 5/5 checks,
1,520 distinct labels. Maximum real/imaginary coordinate discrepancies at
p=4,12,24 are approximately 0.0312338, 0.000122011, 2.97802e-8. These are
float-reference diagnostics; the oracle guarantee comes from its verified
finite expression and radius-plus-rounding argument, not these measurements.
The 128-bit ball midpoints agree with the reference within 2.23e-16; that
does NOT mean the rounded float reference lies inside the narrower balls.

An analytically empty-circuit coordinate is exactly zero in both its ball
and dyadic representations. Numerical near-zero coordinates of the mixed
fixture are separately labeled. Reverse query order, ignored high exponent
bits and context restoration (including a diagnostic exception) pass.
Wrong route sign, moved insertion and deleted imaginary components fail the
prefix comparison. Moving the insertion is NOT reversal of W/K at one
insertion; those two supplied operations commute in this model.

The reference caps its Python label count at 2,000 and each dense array at
32 MiB; its largest full-r matrix payload is 1,600 bytes. This is not a claim
that the Python label list occupies the payload of a hypothetical NumPy table,
nor a total-process RSS limit. Diagnostic calls have a 512-bit cap; they are
not a replacement for the uncapped all-label mathematical oracle argument.

## Exact tiny transition laws and a supplied-wide run

`experiment_verified_sampling` enumerates the complete finite-bit Markov law
with exact Fractions, retaining at most 160 unique states after consumed
exponent bits are removed. Conditional probabilities come from exact integer
CDF counts. No histogram uncertainty is involved. Its TV upper bound against
the ideal final joint law, INCLUDING within-sector work labels, is evaluated
with outward 160-bit final-amplitude balls. The existing independent full-r
reference separately agrees with the marginalized ball midpoints to 6.94e-17.

The final report `out/verified_sampling_20260911T034157376089Z.json` passes
8/8 checks. Vary only requested target TV, with p,L chosen by the planner:

| Target TV | Coordinate bits p | Categorical bits L | Outward observed TV upper bound |
|---|---:|---:|---:|
| 1/4 | 10 | 6 | 0.020521609 |
| 1/64 | 14 | 10 | 0.001392361 |
| 1/4096 | 20 | 16 | 1.819685e-5 |
| 1/1,000,000 | 28 | 24 | 6.948777e-8 |

All laws sum to one EXACTLY and meet the rational planned budgets. The last
budget is 1/1,048,576, below the requested target. The independent float
reference is diagnostic; it does not enter that outward TV calculation.

The supplied-wide case has r=2^61-2, t=63, three specified background rotations
and three reflections (eight histories), fixed seed 624. At target TV=10^-6,
the planner requests p=61 and L=28. The observed run uses 135 vector queries,
69 stochastic updates and 38 precision retries; maximum working precision is
154 bits. Its rational budget is 89432351553/144115188075855872 (about 6.21e-7).
Elapsed sample time in this report is about 0.243 seconds, not an asymptotic
benchmark or a same-accuracy comparison against rejection.

This is concrete evidence that requested absolute accuracy bits and working
mantissa bits are different. The wide run has no orbit/sector/output tables
and no reference cache; it still has ball scalars, exact integers, phase
evaluation costs and backend allocations. Peak traced Python memory during
the tiny sweep was 558,223 bytes; tracing excludes native backend allocations
and was stopped before the wide run. Do not report it as the wide sampler's
total memory or RSS.

The coarse p=0 control keeps mass one with normalized all-zero fallback;
deleting those blocks instead gives mass zero. A diagnostic precision cap
raises, and the float comparator returns no numerical certificate. These
controls expose missing assumptions; they do not close the remaining
same-accuracy certified rejection comparison.

## Preserved failures and audit corrections

- `out/ball_backend_failure_20260911T032959873464Z.json`: the initial fixture
  expected unit norm from a single unit phase divided by sqrt(2); its norm
  squared is one-half. The corrected two-phase fixture replaces it.
- The first dyadic-kernel failure passed Fraction weights to an intentionally
  integer-only API. The corrected reference clears denominators exactly.
  Its first passing report, `out/dyadic_kernel_20260911T033622738308Z.json`,
  mislabeled a uniformization of unequal weights as “deleting zero bins.”
  Filtering exact zeros while preserving weights/labels is valid. The final
  control deletes a whole zero approximate block carrying actual mass.
  Main also fixed a negative-input test whose second failure flag could never
  turn false. Previous reports are retained, not silently overwritten.
- `out/verified_prefix_reference_20260911T033847946488Z.json` overstated float
  tolerance agreement as enclosure containment, called thresholded numerical
  zeros exact, double-counted ignored high exponent labels, and claimed a
  context test it had not executed. Main required explicit corrections. The
  intermediate corrected report `out/verified_prefix_reference_20260911T034131226450Z.json`
  still mislabeled some resource/zero counters; main fixed those too.
- Main's resource-label cleanup introduced a leftover variable NameError,
  preserved in `out/verified_prefix_reference_20260911T034343403566Z.json`.
  It was a verifier bookkeeping failure after measurements, not a failed
  amplitude prediction. The final report above follows the correction.

## Validation and reproducibility

Use the isolated scientific environment from the research root:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_verified_sampling
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_verified_prefix_reference
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_dyadic_kernel
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_ball_backend
```

Core passed before scientific edits (`out/verified_prefix_core.log`). Full lab
and claims suites passed WITH the backend (`out/verified_prefix_test_lab.log`,
`out/verified_prefix_test_claims.log`), including analytic Bell/H amplitudes,
exact-zero coordinates, deliberately inflated valid enclosures forcing a
refinement retry, and raising on capped exhaustion. The lab suite also passed
without the optional backend and explicitly SKIPPED the Arb tests
(`out/verified_prefix_no_flint_test_lab.log`); that skip is not a verified test.

After the shared contraction refactor, the coherent sampling, formula and
compiled-indexed-circuit experiments pass in
`out/verified_prefix_refactor_sampling.log`,
`out/verified_prefix_refactor_formula.log`, and
`out/verified_prefix_refactor_circuit.log`. The other six science suites were
not rerun for these helper additions. Documentation regeneration/checking is
recorded in `out/verified_prefix_docs.log`. No commits were made.

Impact: this closes an implementation gap for a narrow exact-input family;
it does not establish a broadly faster classical simulator. The remaining
research question is specified only in TODO 24, not duplicated here.
