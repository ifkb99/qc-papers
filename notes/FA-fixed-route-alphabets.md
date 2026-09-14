---
code: FA
date: 2026-09-11
title: "Fixed alphabets extend merged simulation beyond the history cap"
outcome: confirmed
claims: [C59, C60, C71, C72, C73]
todo: [14, 24, 32]
---
# FA — Alphabet size is a useful structural parameter

C73 owns the mathematical normal form, uniform support bound, input contract
and word-specific refinement. Main implemented the input-only structural class
and extended the shared support helper; lower-cost agents supplied initial
set, complex-prefix and full-law checks. Main audited the code/predicates and
reran them. TODO 32 remains open for the stronger-envelope integration and
wider task-matched comparison; those are NOT completed experiments here.

## Exact support audit

Main's final `out/fixed_alphabet_support_20260911T081239115867Z.json` passes
8/8 checks. It freezes M=1,000,003 and alphabet (0,1,101), cycles its labels,
and varies k by one from 0 through 32. At every reverse intermediate prefix
the direct reached set lies in the independently enumerated coefficient
cover, and production's tighter L1 count agrees with an independent integer
formula. Exhaustive selected histories are used only at k<=8.

Across the four tested target sectors, peak support at k=8,16,32 is 39,113,415;
the respective global production bounds are 256,1090,4226. Those peaks are
NOT maximized over all M sectors. Odd/even moduli, coincident labels, fixed
points and noninvertible differences include genuine three-label cases.
Growing alphabets give supports 3^m for m=1,...,8 in the frozen pair-word
control, refuting both an alphabet-independent linear and quadratic bound.
This is a counterexample about this support representation, not classical
simulation hardness.

The original support budget uses 1,922,770 of 2,000,000 declared candidate/
update operations. A separately frozen word probe uses 66,877 of 100,000.
Neither counter is a CPU-instruction or bit-operation count. Integer-set
and cumulative-prefix guards precede construction; no large-k path expansion
or modulus-sized scan is needed by these final tests.

## Genuine nine-insertion amplitudes and laws

Main's `out/fixed_alphabet_prefixes_20260911T081039922607Z.json` passes 6/6;
the crash-traced repeat `out/fixed_alphabet_prefixes_20260911T081154846518Z.json`
also passes. Fixtures have t=8, k=9, cyclic labels (0,1,2), pi/4 coherent
rotations, W0 and alternating noncommuting work blocks at all insertions.
They use r=6,b=2 and r=9,b=3, so the tiny full-r baseline is inexpensive.
Small M allows full reached support; these rows test correctness beyond the
old cap, not a storage or timing advantage over dense simulation.

Each fixture checks 10,743 distinct canonical complex prefixes against the
EXISTING independent full-r contraction. High unused exponent bits do not
affect a queried prefix; 87 additional noncanonical requests per fixture
check this aliasing explicitly, with zero observed discrepancy. This is not
enumeration of every redundant API spelling. Largest complex error is
2.17e-16 (binary) and 1.31e-16 (three-state).

The existing transition enumerator checks the entire joint sector/output law,
not a sampled histogram. TV against the independently calculated full-r law
is 6.77e-16 or smaller; all masses agree with one within 4.5e-16. Its test-only
canonical cache retains 8,699 vectors per fixture, not a production lookup
table. Separately returned UNMEMOIZED samples use 37 and 35 prefix queries
for binary and three-state fixtures, with zero rejection/history calls.
Independent reached-set bookkeeping verifies all prefix and uncached sample
operation counters. Legacy comparison uses the actual history sampler/law,
not two copies of the merged implementation.

The full three-state law changes by maximum cell difference .03475 when
intermediate/end work mixers are removed. Controls on the bounded t=4 fixture
give TV .06878 for wrong internal QFT phase and .09384 for wrong reflection
boundary. Default reference-width, excessive opt-in width and legacy k>8
calls are rejected. No claim of finite-TV certification follows from these
floating-point agreement diagnostics.

The reference helpers retain their default t<=4 guard and only accept t<=8
through an explicit keyword, with preallocation and work guards. Conservative
matrix-entry-touch estimates for the wider binary/three-state references are
94,379,616 and 212,362,884, under the declared 2^28 cap. Transition bookkeeping
has a separate bound; it does not include the oracle's contraction work.
The prefix cache allowance is 11,000,832 bytes, with ndarray payload 687,552
or 1,031,328 bytes. The law cache allowance is 8,907,776 bytes; simultaneous
frontiers and dense references are charged separately in the raw report.
These declared allowances are not native RSS measurements. Test-only caching
and dense references are never reported as production table-free storage.

## Stronger baseline found before wider timing

C73's E/O offset recurrence tracks the PARTICULAR chronological word. The
separated-label set probe uses M=2^40-1 and cyclic labels
(0,679535556937,314159265359). Its global envelope at k=4,8,12,16 is
14,48,104,182, versus the generic L1 bounds 16,256,626,1090. Main added actual
reached-set checks at the predicted disjointness witness, including an even
modulus, instead of accepting a witness-exists flag alone. All 19 product-
condition rows attain their envelopes at a checked sector.

The small control M=101, word (0,1,2) cycled through k=8, has ten even and
ten odd offsets. Sector zero reaches only 11 labels, while sector five reaches
20. Thus using the initial-sector count as a global rejection envelope is
unsafe. The additional even-M=1000 check also attains the larger count.

The connection is to additive difference sets: they determine where the two
affine orientation classes collide. The proof belongs in C73. The numerical
pattern is independently checked, but no production sampler uses this tighter
envelope yet. Following the skill's strongest-baseline rule, wider timing is
deferred until this inexpensive improvement can be included and charged.

## Audit trail and limits

- The original support report `out/fixed_alphabet_support_20260911T075625075106Z.json`
  passed while an edge check scanned every sector for fixed points outside its
  reported budget. Main caught the hidden O(M) work. It was replaced by solving
  the two-times-sector congruence directly, with at most two roots. Subsequent
  reports preserve each guard/bound refinement; no earlier report was deleted.
- The first prefix implementation had no memoizer in its expensive transition
  loop, a late cache guard, incomplete work accounting, and a legacy comparison
  between identical reference inputs. These were corrected before the agent's
  passing report `out/fixed_alphabet_prefixes_20260911T080631320535Z.json`.
  Main then added exact counter predicates, unused-bit alias checks, real legacy
  history comparison, uncached sample checks and stronger resource allowances.
- `out/fixed_alphabet_prefixes_main.log` records a main run that exited 139
  without a traceback or JSON report. Python 3.12.10 and NumPy 2.4.6 were
  independently confirmed. The next two runs passed, the latter with
  PYTHONFAULTHANDLER=1; the crash cause remains UNRESOLVED. This is not evidence
  that the documented Python 3.14 issue explains it. The agent's earlier empty
  aborted log is retained separately. The agent reports manually stopping that
  uncached run, after which its wrapper also reported 139 with a zero-byte log.
  That termination does not distinguish a crash from signal handling and is
  not an independently established second spontaneous failure.
- `out/fixed_alphabet_claims.log` preserves a missing-math-import NameError in
  main's new regression, corrected in subsequent runs. It was a test error,
  not a refuted scientific prediction.

The legacy eight-insertion constructor still rejects large history inputs.
Structural tests include snapshot ownership, both merged samplers, memory
rejection, endpoint/invalid requests and 64 identity rotations at the width
cap. Only the nine-insertion fixtures above have new full-law coverage; the
64-insertion identity check is NOT a generic deep-circuit validation.

Validation logs: `out/fixed_alphabet_core.log` (ran first),
`out/fixed_alphabet_lab.log`, `out/fixed_alphabet_lab_no_flint.log`,
`out/fixed_alphabet_claims_final.log`, `out/fixed_alphabet_claims_word.log`,
and `out/fixed_alphabet_legacy_law.log`. The backend-absent lab run explicitly
skips verified-arithmetic checks. The other six science suites were NOT rerun.
Documentation checks are recorded in `out/fixed_alphabet_docs.log`.

Reproduce the bounded experiments from research/:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -u -m experiments.experiment_fixed_alphabet_support
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -u -m experiments.experiment_fixed_alphabet_prefixes
```

Main read Evetts's integer-coordinate word construction in Section 2.1, not
the proofs of all growth-series theorems. C73 records that prior-art scope;
the mathematics here is not positioned as new group theory. No paper,
abstract workshop or sampler default was changed. The open-ended research
goal remains active; TODO 32 owns the remaining experiment.
