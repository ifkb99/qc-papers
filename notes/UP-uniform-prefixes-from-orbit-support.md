---
code: UP
date: 2026-09-11
title: "Exactly uniform output bits survive dense sector mixing, but suffix conditioning remains"
outcome: confirmed
claims: [C53, C57, C76, C77]
todo: [14, 34, 36, 37]
---
# UP — Which work histories can actually interfere in this marginal?

C77 owns the support-separation theorem, strict endpoint condition and exact
floor-sum implementation. This note owns the discovery, bounded tests and
limits. The direction changes the task from compressing the entire state to
predicting a specified low-output marginal. It neither contradicts C76's
dense sector coupling nor establishes an efficient full-output sampler.

Main derived the partial-Fourier/separation criterion before this turn.
Three lower-cost agents supplied an independent proof audit, an initial
integer-count verifier and two independent output routes. Main implemented
the arithmetic helper, audited and strengthened the verifiers, and reran them.
Only the stateless support helper was added to production code; existing
samplers, defaults and propagators were not changed.

## Frozen circuit and the sharper support argument

The NEW sparse-mixer fixture, frozen in TODO 36 before execution, is
N61,a2,r60,b3,t6. The order is independently checked. It has initial W(pi/4),
post-control W(-pi/10) at i4 and W(pi/11) at i5, no other mixers, and G_1
after the i4 mixer. It is not a retuned row from GS/UT. All control powers
keep their original ascending order. Only the requested output-prefix length
d=1,2,3 changes across the primary comparison.

The original conservative radius counts all three mixers, R=6. It certifies
d1 but not d2: the latter's minimum center distance touches 2R. Main noticed
that the terminal W can instead be stripped under the final work trace,
without commuting it through arithmetic. The resulting R=4 certifies d2.
Main derived this before reading the agent's output; the agent had already
completed its original R6 run. Both stages are retained rather than rewriting
the earlier prediction. At d3 neither radius certifies uniformity.

The initial sequential reference found precisely this useful distinction:
d2 was uniform even though the conservative certificate was inconclusive.
A failed sufficient test is not evidence against the property. The sharper
argument was independently audited, then added as an explicit main predicate.

## Exact arithmetic certificate

`out/uniform_prefix_counts_20260911T101413613709Z.json` passes 3/3 checks.
All 9,000 cases over r=1..40,t=0..8,d<=t,R=0..4 agree with direct modular
counts. The family performs exactly 193,600 direct q-loop terms under its
2,000,000 cap; the touching control separately uses three q values and 64
Gaussian-integer quarter-root terms. Width-63 cases use independent quotient
formulas rather than enumerating their exponentially many histories.
The test retains only aggregate rows, bounded edge cases and at most sixteen
mismatch examples; no orbit/prefix numeric arrays are created by the probe.
This does not bound total imported-library/process RSS.

The must-fail test is stronger than a mere off-by-one counter comparison.
Main added a normalized state allowed by the inclusive support promise whose
shared endpoint produces a nonuniform Fourier prefix. It is a different
state from the actual phase circuit, as C77 explains. This establishes that
the weak equality-accepting criterion can make a genuinely false uniformity
prediction, even though the actual circuit happens to pass a sharper test.

Main read the AtCoder floor-sum documentation, the Euclidean implementation
and its CC0 license. The exact-Python implementation credits that standard
recursion and does not import the source's overflow semantics. The experiment
checks the arithmetic, not novelty of floor sums or discovery of the order.

## Two independent full-output routes

The final sequential report is
`out/uniform_prefix_output_20260911T101723682629Z.json` (5/5 PASS).
It assembles full-r branch matrices and calls the EXISTING sequential_path
for each complete output, then groups those same probabilities by y mod 2^d.
The API has no prefix-stop argument, so no fictitious partial-run savings
are reported. The identity-final-shift control changes only the last
controlled translation, leaving its common work gate unchanged; its first
output bit is deterministically even, so the uniformity control fails.

The independent state-vector report is
`out/uniform_prefix_statevec_20260911T101250167022Z.json` (4/4 PASS).
It evaluates the fixed sparse schedule's closed boundary formula, places
orbit amplitudes on physical residue labels and uses existing Circuit/statevec
for the twelve-qubit inverse QFT ONLY. This is not validation of a compiled
N61 arithmetic circuit. Main replaced the initial branch-action wrappers by
the literal fixed matrix formula; no second reusable propagator was added.
Post- and pre-terminal branch states are retained in the same pass.

The complete 64-outcome laws and the control are compared in
`out/uniform_prefix_cross_reference.log`; both independent routes agree
within 2e-15 per outcome. The high-history Gram formula separately matches
every tested prefix within 3.9e-16. Removing the common terminal W changes
the computed Gram matrices by at most 4.5e-16, while preserving both R4
support and d2 orthogonality. Numerical support displays use a 2e-12 cutoff;
the exact support promise follows from the construction, not this cutoff.

Both first-bit and two-bit laws are uniform at floating precision. The
three-bit law is nonuniform with TV 2.317360332e-5 from uniform; this
uncertified outcome is retained, not a theorem that every failed certificate
must be nonuniform. The full 64-point law has uniform-approximation TV
.2346269487, so this is not an accidentally uniform entire output.

Uniform prefix bits are NOT independent of the remaining bits. A lower-cost
read-only audit of the released law found a product-of-marginals discrepancy;
main then included this as a reproducible negative control. Its full-law
TV is .0235642153, and the maximum pairwise TV between conditional suffix
laws is .0703971779. The diagnostic derives the marginals from the enumerated
law; it does not provide a free suffix oracle or benchmark an implemented
compressed sampler. Keeping the correct Fourier feedback is a stronger
next comparison than discarding these correlations.

## Resource and verifier audit

The sequential reference executes 128 calls, charging 165,888,000 units of
t*r^3. It checks the call cap BEFORE each call, not just afterward. Setup
counts separately reconcile six repeated matrices, eleven shifts, 120
modular phase values, fourteen construction products, forty explicit
unitarity checks and two initial matrix-vector calls. Its aggregate numeric
preflight is 4,497,664 bytes under 16 MiB, including both simultaneously
retained fixtures and the instrument's internal working estimate.

The boundary/QFT route evaluates 128 conditional branch columns total:
64 original and 64 identity-shift control. Original pre-terminal columns are
retained during the same pass, not evaluated as another 64 histories. Counts
are 921,600 dense work matrix-vector terms, 15,360 shift entries, 7,680
phase entries, 114,688 Gram terms and 6,272 Gram-probability terms. Two
90-gate QFTs charge 737,280 gate-entry updates under 5,000,000. Its aggregate
numeric preflight is 3,151,872 bytes. These are explicitly named arithmetic/
entry counts, not native FLOPs, bit runtime, total RSS or physical gate cost
for preparing the boundary. No timing advantage is claimed.

Initial reports are
`out/uniform_prefix_counts_20260911T100445237192Z.json` and
`out/uniform_prefix_output_20260911T100426158162Z.json`. Main capped mismatch
examples, removed an unused integer-memory-cap assertion, and tested exact
metadata/loop bounds. In the output verifier, main replaced a tautological
classification-string predicate with law/unitarity/cost tests, corrected a
preflight counting only one of two retained fixtures, and added the pre-call
cap plus exact even-parity requirement. Scientific inputs and tolerances
were unchanged. Intermediate strengthened reports remain preserved.

The first boundary/QFT run passed its numerical checks but failed while
serializing an internal NumPy orbit array. Its traceback is preserved in
`out/uniform_prefix_statevec_initial.log`; NO partial JSON was written because
serialization failed before writing. The corrected initial report is
`out/uniform_prefix_statevec_20260911T100706255837Z.json`. Main subsequently
fixed cloned cumulative control counters, an understated payload estimate,
an extra branch-recomputation pass, and weak retained-array predicates.
The final Gram products and their costs are instrumented inside the loops.
Independent lower-cost audits found no substantive remaining proof/code issue.

## Validation and continuation

Core passed first in `out/uniform_prefix_core.log`. Lab and claims pass in
`out/uniform_prefix_lab.log` and `out/uniform_prefix_claims.log`, including
the new helper's exact/brute and invalid-input regressions plus C77's endpoint
and Fourier checks. Backend-absent lab checks pass in
`out/uniform_prefix_lab_without_backend.log`, explicitly skipping only the
optional verified-arithmetic sections. The other six science suites were not
rerun. System Python 3.12.3, NumPy 2.4.6 and one BLAS thread were used; the
backend-enabled gates use python-flint 0.9.0. Small passing checks do not
settle TODO 34. No host settings, firmware, manuscripts, abstracts or defaults
were changed, and no commits were made.

The documentation gate is recorded in `out/uniform_prefix_docs.log` after
index regeneration. Reproduce from research/:

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_uniform_prefix_counts
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_uniform_prefix_output
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_uniform_prefix_statevec
```

The skill's insistence on the actual marginal and strongest simple baseline
led to both the terminal-W improvement and the correlation check. The exact
prefix certificate is useful progress, not a general breakthrough. TODO 37
alone owns the next conditioned-suffix/feedback comparison; the broader
research goal remains active.
