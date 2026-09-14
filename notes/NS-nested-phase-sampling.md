---
code: NS
date: 2026-09-12
title: "Nested-phase construction savings survive conditional rejection in a bounded fixture"
outcome: confirmed
claims: [C78, C79, C80, C81]
todo: [41, 42]
---
# NS — From progression rows to returned samples

The user resumed research from the paused TODO40 checkpoint using the original
qsim-research workflow. Main implemented `NestedPhaseProgressions`, preserving
the existing sample, conditional rejection and forced-joint methods. Lower-cost
initial testers audited bounds and supplied law/RNG/comparison drafts; main
reviewed, corrected and reproduced the evidence. C81 owns the implementation
contract and algebraic resource bounds. Current research artifacts and the
pre-edit backup are inside research/ (`out/` for temporary material); the
installed skill was only read externally.

## Corrections before accepting evidence

The independent bounds audit confirmed the selected-row reserve and stride,
but identified missing distinct-residue/candidate work in the total local cap.
Main included this, charged phase partitioning and precomputed cut geometry.
Right caching still does not make the selector a phase-cost optimizer; it can
increase calls in highly truncated rows.

The first law report had inadequate row-work allowances and calculated some
acceptance discrepancies without asserting them. The corrected audit verifies
public columns, independent conditional normalization, both normalized proposals
and accepted masses, and reconciles observed counters with pre-call bounds.
Main separated reservations from observed checksums and recorded the near-zero
row's work-weighted joint error. Earlier `out/nested_phase_sampler_*.json`
reports are preliminary unless identified below.

The RNG tester edited the law tester's file despite assignment to another
file, causing an intermediate missing-function failure. Reports remain. Its
recovered draft checked only selected paths, misaligned a retry RNG script,
and used a counter predicate as a negative control. Main replaced it with
complete independent weighted decision enumeration. The discarded draft is
`out/nested_phase_sampler_rng_agent_draft.py`.

Before executing the comparison, main corrected a reference that applied the
late mixer twice, missing counter subcategories and forced-query weighted
divisions, normalization/failure predicates and dense product accounting.
It strengthened the cached fixed-cut baselines and used common seed exponents
and work draws. The unexecuted initial draft is preserved as
`out/nested_phase_sampler_comparison_agent_draft.py`. These are verifier
corrections, not failed predictions about the physical circuit.

## Tiny laws and precision

Authoritative report: `out/nested_phase_sampler_20260912T154859442775Z.json`;
log: `out/nested_phase_sampler_laws_main_final.log` (five checks pass).
Fourteen fixture/construction cases include empty, initial, duplicate, endpoint
and distinct phases, r=b aliases, a non-divisor period with checked oracle
indices, revived cancellations, exact-zero rows and cached alternatives.
References use existing column engines. Public columns and streamed rows are
checked; auto rows exercise both proposals against complete tiny FFT laws.

Maximum amplitude error is below 2.23e-16, joint-coordinate error below
6.94e-18, and accepted conditional discrepancy below 1.12e-16 in positive
fixtures. NC's near-cancellation control retains positive computed norms and
conditional TV about 0.0896 despite amplitude error below 7.86e-17. Its
work-weighted joint half-L1 error is below 4.34e-33. This preserves the
conditional instability rather than declaring a tiny positive row zero.

Per-call reservations total 673,527 named units; reconciled observed counter
checksum is 183,254. These are distinct, overlapping diagnostic ledgers, not
native FLOPs. Additional reference/setup/FFT allowances are identified in the
report. Numeric payload reserve is 409,088 bytes, excluding Python headers,
RSS and arbitrary oracle internals. No finite-bit accuracy claim follows.

## Actual RNG paths

Authoritative report: `out/nested_phase_sampler_rng_20260912T154900507052Z.json`;
log: `out/nested_phase_sampler_rng_main_final.log` (four checks pass).
Main enumerates actual work/component/bit/lift/accept-or-reject transitions for
the exact-root literal columns in the source. Each of two cuts and two modes
traverses 176 weighted paths. Both gcd-lift regimes occur. The full sample
call performs a forced rejection followed by success with one work draw;
exhausted calls retain counters. Fixed-work retries recover the independent
joint FFT law within 8.68e-19. Restarting work instead changes its marginal
by TV 1/12.

The physical fixture has equal component masses, so a separately identified
unequal-mass conditional-row hook exercises different weights through the
actual RNG loop. Its accepted laws agree coordinatewise with an independent
FFT; it is not presented as another physical joint circuit. In total 1,364
sampler calls execute, with observed sampler checksum 102,058 against 725,416
reserved units. Reference/path/setup units are separately counted. Weighted
traversal does not establish the distribution of all finite machine RNG values.

## Matched returned-sample comparison

Authoritative report:
`out/nested_phase_sampler_comparison_20260912T154915832670Z.json`;
log: `out/nested_phase_sampler_comparison_main_final.log` (five checks pass).
The state is the previous NC long-row fixture. Six construction modes include
cached alternatives, the hybrid cut, automatic selection, right caching and
a supplied phase table. Each uses both proposals. Every automatic work row
selects the hybrid cut in this fixture.

These are work-marginal averages from constructed float rows and C80's
**uncapped mathematical envelope**, not averages fitted to four draws:

| Construction | Proposal | Expected attempts | Expected Fourier component terms | Expected progression marginal queries |
|---|---|---:|---:|---:|
| cut 5 or cut 12 | mass | 637.15 | 407,003 | 14,017 (cut 5); 2,549 (cut 12) |
| cut 11 / auto | mass | 149.34 | 22,360 | 2,389 |
| cut 5 or cut 12 | root_mass | 288.07 | 183,318 | 6,338 (cut 5); 1,152 (cut 12) |
| cut 11 / auto | root_mass | 67.35 | 10,045 | 1,078 |

The roughly eighteen-fold Fourier-component saving survives rejection.
Marginal-query savings against cut 12 are much smaller because its stride
has a smaller reduced Fourier register. The full report separates root-mode
square roots/divisions, column work, row construction, selector and phase costs.
Right caching lowers hybrid row phase calls from 578 to 302 and products
from 1,024 to 748. The table baseline charges its values, storage and lookups.

Forty-seven of forty-eight capped requests return. One cut-5 mass-mode request
exhausts its cap; all its work remains in totals. Per-request and per-return
counters including failures are recorded. Four seeds per configuration are
not a statistical runtime estimate; failed calls are not restarted or filtered.

Four complete-output marginal queries, summing every work label, agree with
the existing full-work sequential_path reference within 6.51e-19. Its omitted-
phase control changes checked outputs, and two sequential draws return.
This supplements complete tiny-law tests; it is not a full long-output TV
certificate. Dense matrix products, validation and QR shapes are charged
separately and cannot be equated with scalar geometric queries. No native
speedup over this baseline follows. Full-r cached/sequential approaches remain
relevant at fixed r.

Preflight reserves 11,948,224 numeric bytes, 59,266,012 scalar named units and
288,419,170 dense diagnostic units. Actual scalar budget usage is 26,563,649.
The dense checksum includes overlapping unitarity subcounts; neither checksum
is an additive native FLOP count. No Q-by-r archive or timing sweep is used.

## Validation and next boundary

The fresh core gate passed before implementation. `test_lab.py` passes in
`out/nested_phase_sampler_lab.log`, including distinct-oracle, revival,
varying-stride, pre-copy cap and invalid-input regressions. `test_claims.py`
passes in `out/nested_phase_sampler_claims.log`. An AST comparison with the
saved pre-task source confirms all pre-existing work-first functions/classes
are unchanged. The other six science suites were not rerun.
Indexes were regenerated and all ten documentation checks pass in
`out/nested_phase_sampler_docs.log`; `git diff --check` is clean.

TODO41 is complete at bounded float implementation/comparison scope. TODO42
owns cut selection for total expected sampling work without constructing every
competing row. Certification, arbitrary mixers, wider scaling and native timing
remain separate. No manuscripts, abstracts, existing defaults, host settings
or commits changed in this follow-up.
