---
code: CW
date: 2026-09-11
title: "Component masses refine the work-first rejection envelope"
outcome: confirmed
claims: [C59, C78, C79, C80]
todo: [39]
---
# CW — Count components, then account for their masses

After the C79 long-row checkpoint, main froze a square-root weighting pilot
in TODO39 before implementation. C80 owns the weighted-envelope proof and
its limits. It applies the existing C59 Cauchy-Schwarz argument to the disjoint
progressions; there is no new sampling principle or novelty claim.

An independent read-only audit confirmed both the raw-mass formulas and the
optimal-envelope/cost-model identities. The auditor initially reversed the
fixed-work retry rule, claiming unequal row envelopes would bias work even
when work was held fixed. Main rejected that warning using the conditional
geometric series; the auditor corrected it. A later response also confused
the older r3 work fixture with the new unequal-length Q8 row. Neither error
was incorporated into the proof or production loop.

The tiny pilot fixes one singleton at e0 with mass9/10 and one progression
e2,e4 with mass1/10. Their different Fourier laws make the mismatched-proposal
control meaningful: use the new component weights with the old acceptance.
Main reproduced the initial formula/progression-primitive pilot in
`out/component_weighting_main_pilot.log`. This stage alone was not a test of
the actual component-selection/acceptance loop. Its first failed report,
`out/component_weighting_20260911T124838531641Z.json`, compared an unweighted
union of component paths with a mass-weighted proposal; the corrected report
`...124909168989Z.json` retained component probabilities. Subsequent tests
exercise the shared production loop separately.

Main added the opt-in proposal mode and extracted the SAME conditional loop
into a private row method for isolated RNG tests. The full physical work draw
is validated separately. Existing defaults, manuscripts and host settings
remain unchanged. The evidence and final scoped gates are recorded below.

## Tiny conditional-loop audit

Main's final strengthened reproduction is
`out/component_weighting_20260911T130608099555Z.json`, with log
`out/component_weighting_main_audited.log`: all five predicates pass. It
compares all four fixture FFT laws coordinatewise to the normalized accepted
laws, and the independently weighted actual progression paths to the formula
proposal. The unequal fixture gives mean attempts 2 versus 1.6. Splitting
its two-point progression preserves the target but raises the new envelope
to 1.948528; the proposal TV changes by 0.0625. Exact-zero and equal-mass
controls behave as derived. The normalized deliberately mismatched law has
TV 0.0508644, so this is not a normalization-only control.

The actual shared conditional loop executes 58 calls / 60 attempts: both
proposal modes, component selection, interval bits, both gcd lifts,
accept/reject branches, two-attempt retries and finite-cap failures. Actual
counts are 240 real RNG draws, 60 integer lifts and 240 marginal queries.
Weighted proposal path error is at most 2.78e-17. The aggregate preflight
numeric reserve is 92,672 bytes against a 1 MiB cap; it includes retained
arrays, path records and conversion copies, not process RSS.

The preserved intermediate report
`out/component_weighting_20260911T125619406574Z.json` failed P4: its verifier
only forced acceptance on this strictly positive target, yet required a cap
failure. Main required explicit rejection branches whenever acceptance is
below one and a retry in BOTH modes. These are verifier corrections, not
evidence against the weighted-envelope identity. Main also replaced a
tolerance-based zero-path shortcut with the exact geometric-zero condition,
required observed rather than planned RNG counters, and strengthened the
all-fixture accepted-law and independent-mixture predicates.

## Physical same-law / actual-RNG comparison

Main's final reproduction is
`out/component_weighting_physical_20260911T1318763050Z.json`, with log
`out/component_weighting_physical_main_final.log`: all four predicates pass.
It retains r60,b3,t8,s7, the original two work blocks and phase oracles, varying
only insertion v0..7. Every work row for cycle/dual/auto agrees with the live
independent direct-column FFT reference (maximum conditional TV 2.96e-16).
Complete AUTO queries check both proposal modes at every output, including
q*A*E coordinatewise and EACH row's normalization, not just averages. Maximum
accepted-coordinate error is 4.17e-17. Explicit cycle/dual weighted formulas
check selected outputs; their full amplitude/FFT checks are separate.

The new envelope reduces work-marginal-averaged expected attempts in every
frozen insertion. For example v3 gives 18.6478 versus 7.9737. This is a
mathematical row-mass expectation, not a measured runtime. Both modes reset
seed 624+v and execute 32 draws per insertion, cap 512. All 512 calls finish;
default-mode work/output/attempt/component-count streams exactly match
`out/earlier_phase_sampler_20260911T123806945928Z.json`.

| Observed category, 256 draws per mode | mass | root_mass |
|---|---:|---:|
| attempts | 3,569 | 1,558 |
| Fourier component terms | 59,584 | 26,284 |
| progression marginal queries | 34,662 | 15,174 |
| component-weight preparation terms | 3,709 | 3,686 |
| additional square roots | 0 | 3,686 |
| additional weighted divisions | 0 | 26,284 |

Local/phase/setup work is retained separately in the report. These categories
are NOT interchangeable operation costs and must not be added into a fake
runtime total. Different RNG consumption changes later sampled work rows,
so even the same seeds do not give paired work-label trajectories. The
observed win is not a guaranteed seedwise comparison. The original mass-mode
sampler of the SAME conditional state is the direct strong baseline here;
this experiment does not rank every classical simulator.

The normalized wrong new-proposal/old-acceptance physical law has TV 0.0212528.
Main rejected the earlier passing verifier as insufficient: it compared the
unchanged conditional_probability field, rather than the new q*A*E law, and
averaged normalizations across work rows. Other corrections cached and charged
analytical roots/divisions, reconciled observed helper counters, retained
failed-call counters, fixed cover-specific phase bounds and reset both modes
to the frozen seed rule. Reports using seed 1624+v for the new mode are retained
but are not the final frozen comparison.

Preserved failures include `...1248303087Z.json` (accepted-submass compared
without its normalization), `...1258737994Z.json` and `...1259977937Z.json`
(full-query budgets omitted work/call factors), `...1309776595Z.json` and
`...1310302276Z.json` (sample phase reservation), and `...1311029046Z.json`
(incorrectly inflated phase preflight). All use the prefix
`out/component_weighting_physical_20260911T`. None is silently replaced by a
PASS artifact or interpreted as refuting C80.

Main's final aggregate numeric/verification reserve is 3,090,112 bytes against
16 MiB, including direct FFT/law temporaries, retained row/scalar records and
bounded prior-sample extraction. The large prior report is streamed; only its
bounded sample field is decoded, and prior reference-array views are freed
before the next insertion. This is a preflight reserve, not measured process
RSS. Production sampling receives no precomputed law or orbit table. Reserved
and observed counters are distinguished; total phase counts include the
early/late subcounts and must not be added to them again.

## Scoped regression gates

Core ran first and passed (`out/component_weighting_core.log`). The affected
lab suite passes with the optional backend
(`out/component_weighting_lab_backend_final.log`) and without it
(`out/component_weighting_lab_no_backend_final.log`). Both final lab runs
include the isolated-row actual RNG regression and eight independently
expanded full-joint cases across the two helper classes and proposal modes.
The claims suite passes (`out/component_weighting_claims_final.log`), including
exact-rational envelope/cost checks and the fixed-work retry distinction.
The other six science suites were not rerun. These are scoped float/proof
regressions, not a finite-bit sampling certificate or all-nine-suite claim.
After the C80/CW and TODO40/NC records were updated, reindexing and all ten
documentation checks passed (`out/component_weighting_docs_final.log`).
