---
code: RV
date: 2026-09-11
title: "Certified reverse-vector rejection saves working state but loses the matched timing comparison"
outcome: confirmed
claims: [C63, C64, C65, C67, C68, C69]
todo: [14, 24, 28, 29]
---
# RV — Remove normalization before trying to certify it

C69 owns the finite-error proof, API, termination assumptions and costs.
Main implemented the sampler; lower-cost agents independently audited its
proof/code and supplied complete-law and comparison experiments. Main
strengthened predicates and reran both. The qsim-research skill's existing-
helper and stronger-baseline rules prompted the simpler unnormalized vector
construction: C64 and C63 already control its proposal and boundary acceptance.
No new normalized-trajectory precision theorem was needed.

`VerifiedReverseWork` is opt-in and does not change the default sampler or
add a coherent-history outer option. It retains the selected reverse vector
and current children, constructing exact branches on demand. All forward
states are absent. Precision refinement replays fixed bits and j; rejection
redraws j but keeps the initially chosen coarse sector. C69 explains why
these two labels have different redraw rules.

## Complete finite laws

`experiment_verified_reverse_work.py` fixes §RI's r9,b3,t4 circuit including
W0. It enumerates all four histories, three initial sectors, three boundary
labels and sixteen outputs at targets 1e-3 and 1e-6; target intervals are
checked at P192/P256. There are 48 accepted-law rows, not 24. Each
j-conditional proposal normalizes exactly as Fractions and meets its own
TV bound against enclosed terminal D. The boundary-averaged accepted law
normalizes by its ACTUAL accepted mass, not the ideal one third, and meets
the separately proved accepted-law budget. Success and expected-attempt
bounds are checked from the full finite law, not a sample histogram.

At target 1e-3, maximum outward accepted-law TV is about 2.67e-5 against a
plan bound about 6.10e-4. At 1e-6 it is about 2.88e-8 against about 6.02e-7.
The shared prefix target also agrees with an independent full-r floating
diagnostic. The latter is not the outward certificate. Skipping acceptance
has outward TV at least about 0.235673 in a conditional law. Width-zero
identity rows include exact ZERO ACCEPTANCE NUMERATORS with denominator
one; these are not evidence of a zero denominator. The separate lab tests
explicitly require at least one zero-denominator trajectory.

Authoritative report: `out/verified_reverse_work_20260911T062706832480Z.json`,
4/4 checks; log `out/verified_reverse_work_audited.log`. Main discarded the
unnecessary retention of every full attempt dictionary, added per-j proposal
checks and corrected row counts and zero-case labels. A structural 16 MiB
guard covers retained summarized laws, serialization and current small
dictionaries; it is not native RSS. The diagnostic enumerates a tiny tree;
the production sampler never allocates that tree.

## Matched cost result: storage wins, timing does not

`experiment_reverse_work_comparison.py` uses §UG's supplied wide r,b,t and
fixed five seeds 624–628 at target 1e-6. Timings include worker snapshot
construction plus an accepted sample; the supplied input circuit is built
once outside the timer. All rejected attempts are included. Different random
decompositions are NOT expected to return matched outputs for a shared seed.

| Method | Median seconds | Range seconds | Highest working precision |
|---|---:|---:|---:|
| Full storage | 0.00519 | 0.00514–0.00538 | 202 |
| Checkpoint k8 | 0.00782 | 0.00779–0.00792 | 202 |
| Reverse rejection | 0.01033 | 0.00313–0.02799 | 212 |

Reverse attempts are 5,3,5,2,1. Five timings do not establish an expected
runtime, asymptotic speedup or universal slowdown. They DO fail to justify
replacing either existing method for speed. The binary comparison includes
C65's stronger scalar sampler; its median is also lower than reverse here.

Retained forward matrices remain 64 for full storage and 15 for k8, with
persistent branch counts 126 and zero respectively. Reverse retains neither
and has a bound of four work vectors. Like-scoped conservative working
scalar allowances are 1845,360,240 respectively, INCLUDING reverse matrix
temporaries. Comparing only its vector entries to a baseline's total working
allowance would exaggerate savings. None of these is RSS; supplied gate data,
output bits and scalar bit lengths still count. No width limit was raised.

Tiny endpoint fixtures at widths 0,1,3 include W0 and routes at both ends.
Their complete normalized reverse laws agree with the certified finite-work
baseline within the sum of their budgets, with matching final-sector labels.
A separately forced t32 no-background probe records 30 zero-block fallbacks
across all checked paths; it is not a complete t32 law. Omitting acceptance
differs from the MATCHED independent tiny target by float-diagnostic TV about
0.168811. Resetting to e_j changes normalized conditional weights, not merely
their scale. The stronger outward omission check is in the full-law probe.

Authoritative report: `out/reverse_work_comparison_20260911T063147820197Z.json`,
5/5 checks; log `out/reverse_work_comparison_audited.log`. Main added complete
edge-law comparisons, equal-accuracy predicates, a like-scoped scalar-storage
comparison and an aggregate 32 MiB retained-report preflight.

## Preserved failures and corrections

Main's initial polynomial norm used Arb general exponentiation `x**2`.
Intervals crossing zero then gave NaN and triggered futile precision doubling.
The initial smoke process was stopped after confirming it live; its log is
`out/reverse_sampler_initial.log` (buffered, empty). The bounded diagnostic
`out/reverse_sampler_accept_diagnostic.log` records the timeout traceback.
The repair uses multiplication as in C63; a regression checks finiteness on
coordinate intervals crossing zero. The repaired smoke log is
`out/reverse_sampler_smoke.log`. No timeout-limited run was called certified.

The full-law agent's first independent-reference check failed due to a
normalization mismatch, preserved in
`out/verified_reverse_work_20260911T062505232146Z.json` and its paired failure
report. The intermediate reports stamped `20260911T062522793897Z` and
`20260911T062545497189Z` precede main's stronger predicates and allocation fix.

Comparison logs `out/reverse_work_comparison_test_20260911T062*.log` and
their reports preserve dependency/API and guard-definition failures. Early
passing reports also had a mismatched omission fixture, raw rather than
normalized reset comparison, constructor-excluded timing, misleading default
forward counts, a vacuous rare-case predicate and overwritten harness schema.
Main flagged these; the agent repaired them before its pass stamped
`20260911T062912361783Z`. Main's additional audit/rerun above is authoritative.
These were verifier defects, not evidence against the reverse identity.

## Validation and next boundary

Core passed first (`out/reverse_sampler_core.log`). Backend-present lab and
claims pass (`out/reverse_sampler_lab_final.log`,
`out/reverse_sampler_claims.log`); backend-absent lab passes with explicit
verified-arithmetic skips (`out/reverse_sampler_lab_no_flint_final.log`). Lab
checks child and terminal replay with deliberately widened intervals, exact
zeros, API/plan validation and forbids invoking the full forward builder.
Claims include the exact frozen-j counterexample. The other six science
suites were not rerun. Documentation gate: `out/reverse_sampler_docs.log`.
No manuscript/abstract edits or commits; the broader research goal is active.

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_verified_reverse_work
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_reverse_work_comparison
```

This is a validated memory tradeoff, not a general simulation breakthrough.
The tempting next link is local quantum-instrument error plus exact dyadic
state compression, to avoid the full-tree precision coefficient. That follow-up
is now completed by TODO 29 / C70 / §QD, not certified by C69 itself.
TODO 30 owns the next structural investigation; TODO 24's broader backend
bit-cost issue remains separate.
