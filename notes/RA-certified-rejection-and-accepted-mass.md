---
code: RA
date: 2026-09-11
title: "A same-accuracy rejection comparator works: certify the proposal and charge lost accepted mass"
outcome: confirmed
claims: [C59, C60, C61, C62, C63]
todo: [14, 24]
---
# RA — Certification need not resolve a rare probability relatively

Follow-up: C64/C65 and §FW now provide opt-in verified finite-work/scalar
proposals and clarify the late-mixer limitation of this fixed binary fixture.
The evidence below is the earlier prefix-component comparison, not a claim
that it remains the strongest implementation.

C63 owns the proof and input/resource contract. This follows §NG's precision
audit. The qsim-research skill required separate proposal certification,
complete tiny laws, must-fail controls and matched requested accuracy. Three
lower-cost agents supplied an exact-rational initial probe, a complete-law
verifier and a read-only proof/code audit. Main implemented the sampler,
corrected the verifiers and reran the affected checks. No paper/abstract
edits or commits were made.

## Why not immediately build an exact-real sampler?

Main read [Brassard, Devroye and Gravel, Section 4, Algorithm 4 and Theorem 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC7514202/).
Adaptive fair-bit comparisons can implement exact rejection with computable
probabilities. Their algorithm first constructs a proposal from approximations
at every output; its refinement bound uses that proposal's positive floor.
That is not a free table-free proposal for our circuit. We do not import
their communication/random-bit bound as a quantum-contraction runtime bound.
For the current finite-TV task, C63 instead under-accepts slightly and bounds
the total lost mass. This application makes no novelty claim for rejection
sampling, adaptive comparisons or accepted-measure perturbation.

The useful link to C60 is the weighting: a rare block can have large relative
conditional error while contributing little to the overall law. Rejection
has the analogous accepted-submeasure bound. We pay explicitly for that lost
mass; we do not label a quantized local decision exact.

## Implementation and meaningful baselines

`lab/verified_rejection.py` streams the same C59 finite histories for final
acceptance. The original verified oracle now exposes that existing component
iteration; the ordinary coherent sum and both enclosure modes retain their
math. A chosen proposal history contracts just one deterministic routed
component. Its proposals use C61's VERIFIED prefix sampler, with rational
angle inputs and finite-bit integer kernels. The float finite-work sampler
is not used inside a purportedly certified proposal.

The stronger float finite-work instrument's linear-depth arithmetic cost is
therefore NOT achieved here: component proposal work still rebuilds prefixes.
This distinction matters both for positioning the result and for deciding
what remains worth optimizing. Both compared algorithms receive the same
supplied indexed circuit and promise the same joint (gamma,y) TV tolerance;
neither discovers the order or constructs physical gate decompositions.

## Exact-rational edge probe

`out/rejection_deficit_20260911T042546841544Z.json` passes 4/4 checks.
The endpoint probe covers eleven fixed cases at five word precisions,
including exact-zero/cancellation, negative outward lower endpoints and
denominators down to 2^-1024. It checks the exact deficit inequality and
agreement with the production integer floor. Only L varies in the fixed
two-cell accepted-law sweep. Wrong proposals and uncharged mass deletion
give the required counterexamples, each with TV=1/2.

The initial failure `out/rejection_deficit_20260911T042223122558Z.json`
is preserved. The agent initially set scaled numerator P equal to proposal q,
then compared the resulting law to a different declared target p. This was
an inconsistent fixture, not a failure of the bound. The corrected agent
report `out/rejection_deficit_20260911T042417648171Z.json` remains too.
Main then added negative-lower-endpoint cases, explicit true-value validity
checks and comparison to the production floor. This avoids testing only
intervals whose lower endpoints are already probabilities.

## Full tiny laws, not histograms

The fixed fixture remains r10,t4 with the exact gates in
`experiment_verified_sampling.fixtures()`. All four component Markov laws
are enumerated by the EXISTING exact Fraction transition evaluator, mixed
using the implemented finite-bit history probabilities, then multiplied by
the actual acceptance threshold at each of the eighty joint cells. Every
component, history mixture and normalized accepted law sums to one exactly.

`out/verified_rejection_20260911T042746299715Z.json` passes 4/4 aggregate
checks, including separate proposal error and success bounds. Outward exact
Fraction TV bounds against verified ideal intervals are:

| Requested target TV | Proposal TV upper | Accepted-output TV upper | One-attempt success |
|---|---|---|---|
| 10^-3 | 5.10714e-6 | 2.06971e-5 | 0.3966184 |
| 10^-6 | 6.19037e-9 | 1.77592e-8 | 0.3966582 |

Displayed error bounds round upward; success values are descriptive decimals.
The full exact fractions are in the report. Both P192 and P256 references
pass; neither is an independent propagation algorithm. Independence comes
from the existing full-r float matrix product: ideal target midpoint cells
agree within about 6.94e-17, and all four deterministic-route components
also pass individual independent-product/normalization checks. Floats are
diagnostics, not the certificate. Skipping acceptance has an OUTWARD LOWER
TV bound above 0.1775, a genuinely separating control.

The first agent PASS report, `out/verified_rejection_20260911T042515574744Z.json`,
had two invalid inference checks: it used an UPPER TV bound to infer a
nonzero error, and its proposal check merely required TV>=0. Main replaced
those with a lower bound and the actual proposal budget, checked the success
lower bound and requested target, and compared IDEAL intervals to independent
products rather than using loose accepted-law agreement as the only reference.
The agent's serialization failures in logs
`out/verified_rejection_test_20260911T042409Z.log` and
`out/verified_rejection_test_20260911T042502Z.log` remain.
Main's first stricter rerun also exposed nested Fraction report serialization;
`out/verified_rejection_20260911T042713826848Z.json` retains that failure.
The final report uses the standard harness schema and exact-string conversion.

Reference allocations are tiny and capped before enumeration. The existing
transition cache allowance applies; extra retained dictionaries have a
conservative entry allowance for these frozen precisions. These are resource
guards, not measured process/native memory. No reference table is used by
either wide production sampler.

## Zero/near-zero controls and matched wide timing

`out/rejection_comparison_20260911T042906091975Z.json` passes 9/9 checks
(including descriptive timing bookkeeping, not nine separate discoveries).
At r6,t2, final K(q1,theta), gamma0,y1, the identity component vanishes.
The angle sequence is zero, pi/2^4, pi/2^16, pi/2^64 and pi/2^256; only that
angle changes. All absolute-width decisions terminate at P64 with no retries.
At the two smallest positive angles, the finite acceptance grid rounds a
strictly positive ideal acceptance to zero. The exact outward deficit bound
passes: this is charged approximation, not an exact local sampler.

The supplied-wide input is the SAME rectangular-mode r=(2^61)-2,t63 fixture
as §NG, target TV 10^-6, seeds 624/625/626, two repeats and alternating
algorithm order. Times include construction, planning/history setup, all
rejected proposals, acceptance, refinements and bit requests. There are only
THREE distinct random traces, not six independent seeds. Rejection takes
3,2,4 attempts on those respective seeds, repeated identically.

Median elapsed seconds are about 0.25137 for prefix sampling and 0.12704
for certified rejection: about 1.98 times faster in this bounded comparison.
The first report `out/rejection_comparison_20260911T042322683508Z.json`
had the same qualitative result; it is retained. This is not a crossover,
uniform runtime bound, optimal implementation claim or inference of native
memory usage. Prefix proposals contracting one history can outweigh the
retry cost here; the stronger linear-depth proposal is still uncertified.
Both existing sampler defaults remain unchanged; the comparator is explicit.

An exact two-cell control also shows that a one-attempt cap plus uncharged
fallback can change the law by TV=1/4. This does NOT say that a fixed iid cap
with explicit abort biases the accepted law conditional on success; C63
distinguishes that case from a fallback or query-dependent numerical failure.

## Validation and reproduction

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_verified_rejection
```

Substitute `experiment_rejection_deficit` or `experiment_rejection_comparison`
for the other two probes. Optional python-flint remains outside project
metadata. Seeds are diagnostics, not evidence of ideal unbiased randomness.

Fresh core gate: `out/certified_rejection_core.log`. Affected suite passes:
`out/rejection_test_lab.log` and `out/rejection_test_claims.log` (Arb enabled).
Full-law regressions: `out/rejection_verified_sampling_regression.log` (8/8)
and `out/rejection_coherent_sampling_regression.log` (6/6).
`out/rejection_no_flint_lab.log` passes with explicit backend-test SKIPs;
it is not a verified-arithmetic test. The other six science suites were not
rerun for these helper changes. Documentation gate: `out/rejection_docs.log`.
Current next work belongs only in TODO 24. No general simulation breakthrough
or backend-wide worst-case precision theorem is inferred from these rows.
