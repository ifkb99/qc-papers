---
code: NC
date: 2026-09-11
title: "Several diagonal insertions share a nested-prefix partition"
outcome: confirmed
claims: [C79, C80, C81]
todo: [40, 41]
---
# NC — A phase-free binary interval suggests a useful cut

While completing C80's weighted-envelope audit, main noticed that the
phase-argument residues in C79 are nested binary prefixes, not independent
periodic variables. Main derived a hybrid remaining-history / residue-cycle
construction, then its largest-effective-gap identity. A lower-cost read-only
auditor independently confirmed both, including strict-left versus at-cut
phases, duplicate positions and the final endpoint. TODO40 retains the frozen
pilot record; C81 now owns the proof. No new production class or claim
of novelty has been made.

The lower-cost integer tester's first report,
`out/nested_phase_cuts_20260911T131221029090Z.json`, passes word/stride tests
but constructs its partition only by grouping the reference scan. Main
required an independent k/rho/N/z/start/count construction and comparison of
ALL keys and labels. The largest-gap prediction was frozen after that first
run, so it is a follow-up test, not retroactive evidence from the initial one.
The preserved `...131603348822Z.json` fails after the expanded construction's
preflight was undercounted; `...131626245476Z.json` corrects that reserve.
Main subsequently charged argument calls inside the independent construction,
all factor calls and the actual control work, instead of calling a fixed
reserved witness allowance an observed count.

Main's final report is
`out/nested_phase_cuts_20260911T131801918073Z.json`, with log
`out/nested_phase_cuts_main_final.log`: four predicates pass. All 7,200 frozen
cases have exact original/cut phase words and identical independently
constructed progression partitions. For positions (3,6), the cut-6 structural
factor is 4 versus 16 at cuts 3 and 7. The minimum factors for second insertion
3,4,5,6,7 are respectively 2,4,8,4,2, agreeing with the prederived largest-gap
formula. The explicit unrefined control has phase words (0,0) and (4,0), which
its wrong period-one construction collapses to the same word.

The integer-visit preflight is 3,037,884 against a 5,000,000 cap; actual named
relation, argument, construction and comparison counts are separate in the
report. Numeric payload is conservatively reserved at 65,536 bytes against
1 MiB, including one live case, summaries and conversion copies. There is
no amplitude array, quantum-law comparison, sampler or timing experiment.
This validates a structural mechanism, not actual construction savings:
the tiny fixture's finite-support cap can hide the untruncated bound's gain.
That checkpoint did not yet test phase-weighted coherent amplitudes; the
follow-up below does.

## Coherent-amplitude checkpoint

C81 now owns the independently audited coefficient construction, disjointness,
largest-gap identity and finite-truncation counterexample. Its mathematical
status comes from the local-path proof, not the finite measurements. The
permanent `test_claims.py` regression covers the integer identity, the exact
local-pair-cost counterexample and the Gaussian-integer revival numerator.

Main extended the existing ER `direct_column` reference with a sequence of
same-oracle insertion positions, preserving the single-insertion default.
The lower-cost tester implemented an experiment-local `cut_row`; no new
generic vector propagator or production sampler was added. All phase factors
are applied before summing candidate pairs, including candidates that would
have canceled without those factors.

The initial amplitude reports passed some useful predicates but were not
accepted as the final audit. Main found a global-period-60 reduction in a
constructor accepting an arbitrary period; divisor-period edge cases would
hide it. This was replaced by supplied-period reduction with finite unit-phase
checks, and the edge test was required to include a non-divisor period.
A failed preflight was initially explained by a false claim that complex work
blocks invalidated the geometric residue bound; main restored the proved
R/F bound. Setup/FFT work, pre-call row reservations, failed-call counters,
literal column/work normalization and the explicit C79 duplicate-phase
baseline were subsequently added. The baseline originally reserved the larger
of two covers while the workload preflight assumed its selected cover; main
made these agree and checks actual reference work against its reservations.
All historical failed and superseded JSON reports remain in `out/`.

Main's authoritative amplitude report is
`out/nested_phase_amplitudes_20260911T134901909125Z.json`, log
`out/nested_phase_amplitudes_main_final.log`: four predicates pass. This
compares all work labels for the frozen two-phase family and all selected
cuts against literal columns. Maximum amplitude error is 3.14e-16; selected
geometric Fourier error is 2.09e-15. Column norms and the complete work marginal
normalize, and the duplicate-phase C79 baseline agrees within 2.49e-16.
The wrong-left-coefficient control fails at work0, exponent60 with error
greater than 1.31. This is a named coordinate witness, not just unequal bounds.

The aggregate numeric-payload preflight is 1,744,384 bytes against 16 MiB.
The reserved named-counter checksum is 2,570,586 against 15,000,000, with each
category capped at 5,000,000. Reference, formula, setup, expansion, FFT and
baseline counters are reported separately. Overlapping counter totals are
diagnostic ledger checksums, not additive native FLOPs, timings or RSS.
Short rows confirm C81's warning: the largest-gap structural factor does not
by itself choose the cheapest actual construction.

Core ran first and passed (`out/nested_phase_amplitudes_core.log`). The new
claim regression passes (`out/nested_phase_amplitudes_claims.log`), and the
unchanged-default ER experiment passes after its reference extension
(`out/nested_phase_amplitudes_legacy_reference.log`, report
`out/earlier_phase_cycles_20260911T132833197145Z.json`). The other seven science
suites were not rerun; shared production helpers are unchanged at this stage.
No manuscript, abstract, default sampler, firmware or host setting changed.

## Edge coverage and a retained precision failure

The initial broad edge prototype exceeded its assigned scope, wrote a second
vector propagator, and did not exercise the actual cut constructor on the
canceled-candidate fixture. It was replaced by a small audit using only the
existing indexed literal-column reference. Main subsequently corrected the
per-row reservations, activated cumulative reference counters, reconciled
reservations with the whole preflight, and added column/work normalization.
Small positive rows are no longer classified as zero using a tolerance.

The positive fixtures cover empty, initial, duplicate and final insertions,
r=b aliases, and a non-divisor period with an oracle range assertion. The
complex three-cell matrices are Fourier blocks, with a column-phased middle
block; these are separate edge fixtures, not another point in the r60 sweep.
Actual cut-row amplitudes, disjoint norms and complete normalized conditional
FFT laws are checked at every cut and work label of those fixtures.

The retained failure `out/nested_phase_edges_20260911T134816675520Z.json`
used a ROW-phased Fourier middle block instead. Its r=b=3 work0 row is nearly
canceled. Changing to column phasing produced nondegenerate positive cases,
but main retained the original case as an explicit negative precision control
rather than deleting the failure. Small absolute amplitude error does not
ensure stable conditional normalization there; this is not appreciable error
in the overall joint law because that work event has tiny computed mass.

Main's final report is
`out/nested_phase_edges_20260911T135924003166Z.json`, log
`out/nested_phase_edges_main_roundoff.log`: five predicates pass. Positive
amplitude error is at most 1.25e-16 and conditional-law TV at most 1.60e-16.
Deleting the actual revived candidate changes its amplitude by about 0.7071.
The separate near-canceled control has amplitude error 7.86e-17, row norms
around 2–3e-31, yet conditional TV about 0.0896. These finite observations
are now in the reproducible report; they are not an exact-zero decision or
a numerical certificate. Preflight reserves 86,033 named units and 269,056
numeric bytes against the original 1,000,000-unit / 4 MiB caps. Actual named
units are 23,657; no native timing inference is made.

## Longer progressions: actual construction savings

After the tiny/edge audit, main froze the selected-row s12,t13 discriminator.
An independent integer audit predicted the exact grouping/coefficient counts
and geometric candidates before any amplitudes were measured. The independent
reference uses the initial/late-mixer displacement cone, not the cut formula,
to prove every unqueried column is zero for the selected work label.

The lower-cost tester built the experiment using the same cut constructor and
ER literal reference. Before execution, main required explicit-Q geometric
queries, own-normalized laws and TV, correction of doubled FFT reservations,
complete reservation reconciliation, and a work0-only control. Main further
added finite checks, charged scatter/comparison entries, checked observed
construction counts against the prederived integers, and reported actual
same-row comparisons rather than ratios to loose upper bounds.

The authoritative report is
`out/nested_phase_long_rows_20260911T140226518143Z.json`, log
`out/nested_phase_long_rows_main_final.log`: four predicates pass. All 8,192
amplitudes on each selected row agree within 3.52e-16; full normalized-law TV
is at most 1.74e-16. Selected geometric amplitudes agree within 1.61e-13.
The wrong-left control has a work0, exponent60 witness exceeding 1.25.

For work0, cut11 has 144 nonzero progressions of length 4–5, versus 614
singletons at cuts5 and12. Its grouping/coefficient visits are 36+288,
compared with 2304+1230 and 18+1230, respectively. Across work0,1,59 the
all-left-to-hybrid local-pair ratio is 3.83–3.86 and its phase-call ratio
4.23–4.26. The earliest-cut local-pair ratio is 10.88–10.91. These are the
named counts of the same naive constructor, not runtime ratios or superiority
over specialized caching, an orbit table, or a full-r sequential simulator.
The measured number of nonzero components is not an exact symbolic-zero
certificate; the independently proved/precounted geometric counts are separate.

Only 2,865 support-compatible literal columns were queried; three selected
Q-entry rows and their FFTs were retained, never a full Q-by-r table. The
aggregate numeric-payload reserve is 3,261,640 bytes against 4 MiB, including
rebindings, scalar records and conversion copies. The overlapping named-work
preflight is 2,546,732 against 15,000,000; category caps remain 5,000,000.
Actual reference block products are 515,700, distinct from row construction.

## An exact selector, not just a largest-gap heuristic

The independent geometry counts suggested a simplification to main: slice
occupancies differ by at most one, so an integer refinement period either
clips all residue classes or clips none. Main derived C81's exact T/S minima;
an independent lower-cost proof audit confirmed them, including aliases and
short slices. This removes the need to enumerate period cycles just to choose
a cut. C81 owns the proof, cost-model restrictions and complexity statement.

The initial selector pilot was reduced before measurement after its preflight
exceeded the frozen 100,000-visit cap. A subsequent verifier KeyError and
preliminary pass with incomplete accounting are retained. Main corrected a
monotonicity boundary, replaced a symbolic grouping charge by actual pair
enumeration, charged selector comparisons and retained counters on failure.

Main's final report is
`out/nested_phase_selector_20260911T140027539754Z.json`, log
`out/nested_phase_selector_main_final.log`: three predicates pass. All 20
fixture/phase families (110 cuts) match both integer closed forms, exercise
both clipping branches and retain a minimum among the phase-position/end
cuts. The known largest-gap misranking remains the must-fail control.
Preflight is 90,306 named visits and 196,608 numeric bytes; actual visits are
50,136. These are integer checks, not amplitude or phase-oracle evaluations.

The permanent claim regression adds another 270 small all-cut cases and
passes in `out/nested_phase_selector_claims_final.log`. Its first run failed
only the hard-coded expected case total (276 instead of 270); all preceding
identity assertions passed. That failed log remains
`out/nested_phase_selector_claims.log`. Core-first and legacy-reference gates
are recorded above. The other seven science suites were not rerun because
production helpers are unchanged. TODO40 is complete at this mathematical /
bounded amplitude scope; TODO41 owns actual sampler integration and its
stronger end-to-end comparison. No generic breakthrough or novelty is claimed.

Final documentation indexes were regenerated and all ten documentation checks
pass (`out/nested_phase_docs_check.log`); `git diff --check` is clean. The
qsim-research workflow kept the actual-count comparison and ill-conditioned
negative control in scope while deferring sampler and certification claims.
