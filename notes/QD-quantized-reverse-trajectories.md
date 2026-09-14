---
code: QD
date: 2026-09-11
title: "Local trajectory error removes the tree precision factor, but exact integer compression loses on runtime"
outcome: confirmed
claims: [C60, C64, C65, C67, C68, C69, C70]
todo: [14, 24, 29, 30]
---
# QD — A precision improvement, not a speed improvement

C70 owns the proof, implementation contract, integer bounds and limitations.
Main implemented the bounded-grid sampler; lower-cost agents independently
audited the algebra/code and built the finite-law/comparison probes. Main
read their code, corrected predicates and reran both. The qsim-research
skill required complete finite laws, a strongest matched baseline and
preserved failures; these prevented treating fewer precision bits as a win
in runtime. Defaults and manuscripts remain unchanged.

## Frozen finite-law audit

Authoritative main report:
`out/quantized_reverse_work_20260911T065458560042Z.json` (5/5 checks).
The frozen RI family has r=9,b=3,t=4, initial W0, W1/W3, two route
insertions, all four deterministic histories, all three initial sectors,
all three boundary labels and all sixteen outputs. Requested TV is varied
between 1e-3 and 1e-6. P192/P256 vary the independent target enclosures;
the actual finite kernel is held fixed, not resampled with higher precision.

For EVERY boundary label, the quantized proposal is checked against outward
ideal norm-mass intervals from C69's raw vectors. Its compressed terminal
denominator is not the ideal proposal mass and is never used as such.
Accepted submeasures are normalized separately for each initial sector.
Both their normalized laws and success masses meet C70's plans.

| Requested TV | Largest per-j proposal TV upper bound | Proposal budget | Largest accepted-law TV upper bound | Accepted-law plan |
|---|---:|---:|---:|---:|
| 1e-3 | 9.45e-6 | 9.92e-5 | 2.82e-5 | 8.13e-4 |
| 1e-6 | 9.19e-9 | 9.69e-8 | 3.10e-8 | 7.93e-7 |

The independent full-r diagnostic agrees with the target interval midpoint
to at most 2.78e-17; it is a floating diagnostic, not the certificate.
Omitting terminal acceptance has TV lower bound about 0.235674. Fixing one
boundary proposal instead of mixing j changes the proposal by about 0.236.
The separate exact claims regression covers the stronger frozen-j-through-
rejection error. No empirical output histogram is used as a TV certificate.

## Matched costs

Authoritative SERIAL main report:
`out/quantized_reverse_comparison_20260911T065650763335Z.json` (5/5 checks).
It was run after the law and lab jobs finished; earlier overlapping timings
are retained but not used for this comparison.

The supplied UG fixture has r=3*(2^60-1), b=3,t=63, the same repeated
rational-pi mixer at every control and no routes. All methods request TV
1e-6; setup and every rejected attempt are timed for seeds 624–628.

| Method | Median setup + returned sample | Maximum observed backend precision |
|---|---:|---:|
| Quantized reverse | 0.0315 s | 64 bits |
| Raw reverse | 0.00947 s | 212 bits |
| Full forward/effect | 0.00495 s | 202 bits |
| Checkpoint spacing 8 | 0.00778 s | 202 bits |

Both reverse methods used attempts [5,3,5,2,1] on these seeds. The quantized
grid used p=34: maximum dynamic coordinate bits 35, initial-vector bits 36,
raw-child bits 70, weight bits 141 and terminal ratio operand bits 140.
No operator refinement or trajectory replay occurred in these finite rows.
These observations do not bound all backend calls or demonstrate a native
memory improvement. The complex-coordinate working allowance is 240, not
240 bytes or a count of fixed-width machine scalars. C70 explains its units.

The separate binary scalar fixture still favors C65 (median 0.000150 s).
The comparison also enumerates binary r=10 width-1/3 endpoint laws and
ternary r=9 width-0/1/3 laws. Each fixed sector is normalized before mixing;
TV against the finite-work kernel is bounded by the SUM of both certificates.
Controls show integer growth without compression and wrong conditional
weights if norms are computed after projective compression.

This settles TODO 29 at its bounded scope. Exact Python integer/rational
compression costs outweigh lower verified operator precision here. Do not
optimize this representation further or change defaults on these results.

## Failures and main-audit corrections

- Initial lab runs had a missing closing parenthesis in the newly added
  planner regression. `out/quantized_reverse_lab.log` and the corresponding
  `_no_flint.log` preserve that syntax failure; fixed/final runs pass.
- The first finite-law preflight counted discarded attempt payloads as
  retained storage and exceeded its guard. The corrected probe budgets its
  actual retained summarized-law/report payload, not every transient attempt.
  Failure: `out/quantized_reverse_work_failure_20260911T064652050406Z.json`.
- An incorrect expected integer helper value made the helper control fail.
  Failure: `out/quantized_reverse_work_failure_20260911T064717105213Z.json`.
- An aggregate-cost predicate multiplied an already aggregate block count by
  attempts again. Failure:
  `out/quantized_reverse_work_failure_20260911T065135012062Z.json`.
- Main found the endpoint verifier pooling accepted mass across initial
  sectors, which is not the implemented retry law. Per-sector normalization
  and the summed-certificate predicate replace that comparison. Earlier
  passing reports are historical, not evidence for the corrected endpoint law.
- Main separated the initial vector's finer grid/bit bound from the dynamic
  state's bound, added missing per-j and actual-child-size checks, and added
  binary endpoint cases. The agent's first proposal check covered only the
  uniform-j mixture. Final reports above supersede those weaker predicates.
- Main added positive rare-child, exact-zero child, local refinement without
  past-state changes, and no-hidden-forward-builder regressions. An initial
  audit worry about all-zero rounded weights was inapplicable: these child
  weights are exact integer norms. The finite CDF cannot select a zero weight.

## Validation and reproduction

Core ran first and passed: `out/quantized_reverse_core.log`. Scoped gates:
`out/quantized_reverse_lab_final.log`, `out/quantized_reverse_claims.log`,
and `out/quantized_reverse_lab_no_flint_final.log`. The other six science
suites were not rerun. Documentation gate: `out/quantized_reverse_docs.log`.
No commit or manuscript/abstract edit was made.

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_quantized_reverse_work
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_quantized_reverse_comparison
```

The useful next lead is structural rather than another precision sweep:
coherent history count may greatly exceed the number of backward-reachable
sector labels. TODO 30 owns its proposed derivation, first route-set audit
and unimplemented coherent-amplitude/sampling test. That lead is not a
result established by C70. TODO 24's backend-wide complexity boundary remains
open and deferred; the user's broader research goal remains active.
