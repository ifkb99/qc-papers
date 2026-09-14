---
code: FB
date: 2026-09-11
title: "Keeping Fourier feedback improves the baseline, but coherent high histories remain detectable"
outcome: confirmed
claims: [C56, C77]
todo: [14, 37, 38]
---
# FB — Uniform prefix is not a license to forget the conditional state

C77 owns the normalized conditional-boundary formula. TODO 37 froze the
comparison before measurement. This note owns its discovery and evidence.
The fixture is unchanged from UP; only the additive phase k=0,1 changes in
this experiment. No angles, insertion, width, base or modulus were retuned.

The prior product-of-marginals control discarded ordinary QFT feedback as
well as harder coherences. Main derived a stronger comparator: condition on
the uniform two-bit prefix, retain the product low-input feedback phase, and
dephase only the high INPUT histories. In each diagonal history term the late
work unitary is common and cancels under the work trace. The comparator then
reduces to the four-control periodic-sector background, irrespective of the
late physical phase. This is an approximation, not a uniformity corollary.

A lower-cost independent audit checked the bit convention, normalization and
work-unitary cancellation before the tests. The displayed conditional state
is normalized only because the sharper C77 certificate proves uniform prefix;
an arbitrary circuit would need its actual prefix probability divided out.

## Two independent numerical routes

Main's final reports are:

- `out/conditional_feedback_boundary_20260911T104152447926Z.json`, 5/5 PASS;
- `out/conditional_feedback_sectors_20260911T104127516005Z.json`, 3/3 PASS.

The boundary route uses the literal frozen work formula and existing
Circuit/statevec inverse QFT, both for the complete six-bit state and the
exact four-bit conditioned boundary. This is QFT-only validation of supplied
states, not a compiled physical arithmetic circuit. The second route uses
existing full-r branch matrices and sequential_path. It independently
reproduces the comparator by averaging existing three-state sector shifts.

The cross-reference record `out/conditional_feedback_cross_reference.log`
compares all complete target/omission probabilities and every candidate
conditional probability: agreement is within 1.3e-16 per entry. Exact
conditioned-boundary laws agree with complete-QFT conditionals within 5.1e-16
TV. Prefix masses are uniform at floating precision. Full-law comparisons
also equal the average conditional TV, checking y=z+4w interleaving.

The feedback-aware candidate's full-law TV is **.01252648884**, versus
**.03459335156** for omission of G_1. Its four conditional TVs range from
.0121187 to .0127832. It is better than the previous UP independence
comparison, but misses BOTH frozen requested tolerances 1e-3 and 1e-2.
These were open classifications, not predicates rewritten after failure.
Feedback explains some of the discrepancy, not all surviving interference.

The normalized synthetic flat low-input state exercises the SAME Q=64
feedback denominator at z=1. Its measured finite Fourier law matches the
independent geometric formula; omitting feedback gives a delta at w=0 and
TV .1887792. The sector route also tests the opposite sign, whose TV from
the correct law is .0674581. These are genuine convention/implementation
controls, not claims about physical circuit complexity.

## Verifier corrections and resources

The initial sector run failed JSON serialization of a NumPy boolean after
its checks. The raw traceback remains in
`out/conditional_feedback_sectors_initial.log`; it records an agent run on
the default Python 3.14 environment. Main's authoritative runs use explicit
system Python 3.12.3 and NumPy 2.4.6, as below. An initial (z,w)/contiguous-y
reshape mistake was also corrected; intermediate reports remain retained.

Main rejected the first synthetic control because it evaluated an independent
integer-shift formula without exercising the candidate's feedback code.
The replacement uses phase_pair and sequential_path, with the same fractional
feedback as the physical comparator. The first corrected-control report
`out/conditional_feedback_sectors_20260911T103554342088Z.json` retains its
failed control predicate. Subsequent checks compare the actual geometric laws,
not an inappropriate peak-location assumption. Main further instrumented
geometric terms, checked setup-count equality and required two-sided norms.
Joint arrays now explicitly distinguish (z,w) tables from actual y ordering.

The boundary route preserves initial failures in
`out/conditional_feedback_boundary_20260911T103342147537Z.json` (orbit vector
versus padded work dimension) and
`out/conditional_feedback_boundary_20260911T103400470467Z.json` (missing
1/sqrt(Q) in the full reference). Its first passing report is
`out/conditional_feedback_boundary_20260911T103420067157Z.json`.
Main then strengthened finite-TV/key-presence predicates into law, interleaving
and actual-cost checks, increased the underestimated statevec buffer reserve,
required the omitted synthetic delta, and removed a wrong-shaped zero-prefix
fallback. Intermediate audit failures are also retained in the output folder.
Specifically, report `out/conditional_feedback_boundary_20260911T104021936177Z.json`
exposes a suffix-shape expectation of four instead of sixteen, a missing 1/H
in the interleaving verifier and an overcounted feedback-phase counter. Report
`out/conditional_feedback_boundary_20260911T104048504516Z.json` retains the
remaining wrong-shape audit after the other two repairs. The final predicates
check the corrected observable dimensions, not relaxed tolerances.
No scientific input or target accuracy changed during these verifier repairs.

The final sector route charges **1,520** forced instrument calls and
**221,322,432** named width*dimension^3 units. It explicitly enumerates all
twenty sectors for the tiny law; a sampled sector is not confused with a
free complete mixture. Setup separately records full-r constructors, initial
matrix-vector work, sector shifts, phase-pair entries, omission products and
768 independent geometric terms. Its conservative aggregate numeric preflight
is **6,917,024 bytes**, under 16 MiB. Calls and work are capped before execution.

The boundary route charges **128** work columns, **921,600** dense matvec
terms, **16** QFT calls and **1,425,408** gate-entry updates. Conditional
transforms, feedback phases, shifts and setup are separately instrumented
and reconciled. The numeric preflight is **1,771,520 bytes**, including
statevec's copied input, rotation/Pauli temporaries and retained boundaries.
These are named arithmetic/entry counts and numeric allocation reserves,
not native FLOPs, measured peak RSS or an end-to-end timing benchmark.

## Validation and next connection

Core passed before this science work at the C77 checkpoint. The affected
claims gate passes in `out/work_first_claims.log`, including C77 and the
neighboring C78 follow-up. Lab/backend checks from UP remain the scoped
production-helper validation; no shared production helper changed in FB.
The other six science suites were not rerun. Documentation is reindexed and
checked in `out/work_first_docs.log` at the combined handoff. No host setting,
firmware, manuscript, abstract, default or commit changed.

Reproduce from research/:

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_conditional_feedback_boundary
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_conditional_feedback_sectors
```

The skill's strongest-baseline and actual-output rules motivated the feedback
comparison and exposed its remaining error. While analyzing those coherences,
main noticed that WORK-first conditioning can leave short exponent rows or
few arithmetic progressions. C78/WF own that distinct positive result;
TODO 38 owns its still-unfinished scalable implementation/comparison.
