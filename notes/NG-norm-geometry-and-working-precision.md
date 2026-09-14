---
code: NG
date: 2026-09-11
title: "Norm geometry explains excess enclosure precision, but extra bookkeeping loses the timing comparison"
outcome: mixed
claims: [C61, C62]
todo: [14, 24]
---
# NG — Working precision includes wrapping and retry-policy overhead

C62 owns the norm-error theorem and conditional precision model. This follows
§VP's requested-accuracy/working-mantissa discrepancy. The qsim-research skill
required independent controls, matched accuracy and explicit resource/input
semantics. Three lower-cost agents tested precision growth, isolated wrapping
and norm prefixes. Main implemented/audited the opt-in mode, corrected the
verifiers, and ran the complete-law/timing comparison and regressions.

## The original wide case did not intrinsically require 154 bits

The previous maximum was an ADAPTIVE DOUBLING ENDPOINT, not the minimum
sufficient mantissa. The initial value was P=77 for p=61 coordinate bits.
Main's earlier conversational shorthand that queries “needed 154 bits”
overstated what had been measured. §VP's report correctly records execution;
it must not be read as a minimum-precision result.

The first retried query was selected diagnostically from fixed seed 624,
then frozen while P varies by one from 64 through 160: sector
1038651310644072125, exponent 2246884356896903187, stop 63, reflection
boundary, measured 45, output 16994140961528. This is not unbiased label
sampling or a worst-case search.

`out/precision_growth_20260911T040055026245Z.json` passes 3/3 checks. This
query's P=77 radius is about 2.32550e-19, just above target 2^-62, while
its midpoint discrepancy against P=512 is about 9.49546e-27. P=78 passes.
Scanning all 38 retried queries on this fixed trace, first passing precisions
from P=77 upwards lie between 78 and 89. The doubling policy overshoots those
thresholds. No bound on unqueried labels follows.

The experiment also checks 1,520 distinct tiny labels at P=64,128 against
P=512 midpoints, as precision diagnostics. It uses no amplitude denominator.
Main removed a float-square-root/Fraction detour, strengthened a count-only
convergence predicate, and added the retry-threshold scan. The earlier agent
report `out/precision_growth_20260911T035545724644Z.json` is retained.

## A known numerical-analysis connection

Main read FLINT's complete enclosure-quality/precision guide and Acb's
rectangular definition, linked in C62. A rotation enlarges an axis-aligned
error rectangle even though Euclidean distance is preserved. Alternating
a rotation and its inverse isolates representation overestimation from
physical amplification or difficult normalization.

With initial axis radii eps=2^-20, P=256 and alternating rational-pi
45-degree phases, sequential radii grow approximately as 2^(n/2)*eps through
64 steps; the final inflation is about 4.29497e9. Applying the combined phase
once gives eps on even identity steps and sqrt(2)*eps on odd steps.
`out/interval_wrapping_20260911T040054469855Z.json` passes 5/5 checks, using
exact squared-radius comparisons and a declared small allowance for the
backend's fixed-precision radius rounding. Floats only display ratios.

Main strengthened loose upper-bound/symmetry predicates to test the stated
growth and removed a fictitious 64-byte scalar-memory guard. Scalar-object
counts are not bytes or RSS. The initial report
`out/interval_wrapping_20260911T035738698435Z.json` and earlier failure
artifacts remain. An initial agent derivation also treated complex phases
as nonexpansive in coordinate-radius norm; main's pi/4 counterexample fixed
that mistake before the derivation entered C62.

## Opt-in norm transport

`VerifiedReflectionCircuit(..., enclosure_mode="norm")` uses the shared
finite contraction's new arithmetic-step callback. It neither reorders gates
nor introduces another propagator. Each exact contraction acts on an exact
dyadic midpoint, adding its local Euclidean error bound to a separately
retained prior radius. Reinsert that radius as a coordinate box only once per
history, before verified coefficient multiplication and coherent summation.

The first prototype used the square root of summed squared radii. Main
replaced it by twice the maximum coordinate radius, conservative for four
real coordinates and cheaper to compute. All rounding/error terms remain
charged. Rectangular mode stays default; lower precision does not guarantee
faster execution or tighter bounds on every input.

`out/norm_prefix_20260911T040237278315Z.json` passes 3/3 aggregate checks:
both modes, all 1,520 tiny labels, p=4,12,24; independent full-r FLOAT
diagnostics; exact Fraction outward comparisons to higher-precision boxes;
analytic zero/Bell-H cases; and a midpoint-only control excluded by the
higher-precision enclosure. The high-precision rerun is not an independent
algorithm; the full-r product is the independent algorithmic reference.
Main corrected the probe's float-to-Fraction detour, added an outward error
predicate, and replaced hypothetical NumPy label-payload accounting by a
conservative bounded-cache allowance. The prior report
`out/norm_prefix_20260911T040013602667Z.json` remains historical.

## Same accuracy, fewer bits, slower code

Main reused the EXACT Fraction transition enumerator and outward final-ball
TV calculation. At target TV 10^-3,10^-6,10^-12, both modes' complete tiny
laws normalize exactly and obey the SAME planned budgets. Their exact
finite-bit laws coincide in these fixtures, not necessarily for arbitrary
inputs/rounding ties. The finest outward TV bound is about 6.65523e-14.

`out/norm_sampling_20260911T040240521257Z.json` passes 9/9 checks, including
descriptive timing-accounting checks, not nine independent discoveries.
The fixed wide circuit uses seeds 624,625,626, two repeats and alternating
mode order. Rectangular mode reaches P=154; norm mode P=77 and avoids the
recorded retries. Inputs, requested TV, seeds and integer kernels match.

Median times are approximately 0.26066 s rectangular and 0.34938 s norm:
norm mode is about 1.34 times SLOWER here. Extra local radius bookkeeping
outweighs saved precision/retries. These are bounded timings, not universal
performance laws or memory measurements. Keep the strategy opt-in; do not
sell this precision reduction as a general simulation speedup.

## Validation and reproduction

Use the command below, substituting `experiment_interval_wrapping`,
`experiment_norm_prefix`, or `experiment_norm_sampling` for the other probes:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_precision_growth
```

Core passed first (`out/precision_growth_core.log`). Full lab and claims pass
with the backend (`out/precision_growth_test_lab.log`,
`out/precision_growth_test_claims.log`). Existing verified and float coherent
complete-law experiments pass after the callback refactor
(`out/precision_growth_verified_sampling.log`,
`out/precision_growth_coherent_sampling.log`). Other six science suites were
not rerun. The lab suite also passes without the optional backend, explicitly
skipping verified/norm oracle tests (`out/precision_growth_no_flint_lab.log`);
those skips are not additional verified tests. Documentation regeneration/checks are in
`out/precision_growth_docs.log`. No paper/abstract edits or commits. TODO 24
alone specifies next work; the open-ended research goal remains active.
