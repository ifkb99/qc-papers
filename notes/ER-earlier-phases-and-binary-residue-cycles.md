---
code: ER
date: 2026-09-11
title: "Earlier-phase progression covers now support an opt-in sampler and long-row validation"
outcome: confirmed
claims: [C59, C77, C78, C79]
todo: [14, 39]
---
# ER — Test the gate-schedule boundary of work-first conditioning

After implementing C78, main chose a changed-schedule experiment instead of
another known-order size record or rejection-envelope micro-optimization.
An earlier diagonal phase is evaluated before the remaining low controls:
its argument is no longer fixed directly by the final work outcome. The
binary prefix of the exponent, however, walks through a finite modular cycle
along each orbit progression. C79 owns the elementary proof and its limits.
This is a useful connection to the project's 2-adic theme, not a claim that
the original arithmetic-support theorems automatically transfer to sampling.

Main froze the derivation and one-parameter test in TODO 39 before the
lower-cost pilot. A separate lower-cost algebra audit confirmed the cycle
and the two constant-coefficient cases. Main caught an important caveat in
the first informal count: the new phase can revive an old canceled
coefficient, so one must rebuild from geometric candidate residues rather
than the old nonzero components. The auditor independently confirmed a
Hadamard/parity-phase counterexample; C79 and exact claim regressions retain it.

## Bounded coefficient experiment

The frozen circuit is N61,a2,r60,b3,t8,s7,H2 with the former W0/W1 and late
G_1. Add one identical earlier G_1 after v=0,...,7 low controls. Only its
insertion position changes. The reference uses literal cyclic vector shifts,
pointwise diagonals and block-matrix multiplication in indexed work space;
the separately contracted formula constructs local columns and refined rows.
Neither formula reads a supplied orbit table. This is not a compiled physical
arithmetic reference or an additional generic propagator.

Main's final report is
`out/earlier_phase_cycles_20260911T112847299154Z.json`, **5/5 PASS**, with log
`out/earlier_phase_cycles_main_final.log`. It records system Python **3.12.3**,
NumPy **2.4.6**, with one BLAS thread. All **2,048** literal columns and
**480** refined rows agree across the eight insertion positions. Maximum
column/row errors are **2.29e-16 / 2.78e-16**, column normalization error
**6.67e-16**, and row norm-sum error **3.56e-15**.

| insertion v | generic cycle bound P | cycle used | maximum nonzero components |
|---|---|---|---|
| 0,1,2 | 1 | 1 | 10 |
| 3 | 2 | 2 | 20 |
| 4 | 4 | 4 | 25 |
| 5 | 8 | 8 | 25 |
| 6 | 16 | 16 | 25 |
| 7 | 32 | 1, by the endpoint identity | 10 |

The interior maxima are finite-row observations, not tight asymptotic lower
bounds. These rows contain only two or three terms per original progression,
so refinement can saturate at the available exponent labels before the period
bound is reached. No claim of large-P efficiency follows from this fixture.
The report retains explicit coefficient-pair witnesses for each interior
position, not just a Boolean variation flag. For example, at v=3,j=0,h=0,
the coefficients at exponents 0 and 60 differ in absolute value of their
complex difference by **1.32734**. They cannot be represented by one constant
coefficient on that original progression. Removing the earlier phase also
fails the amplitude comparison; that alone would NOT prove an output-law
difference, since some work-only phases are unobservable after tracing work.

## Audit, failures and resource limits

The initial failed reports are retained:

- `earlier_phase_cycles_20260911T112045001326Z.json`: an over-budget initial
  work estimate rejected before allocation; the bounded implementation and
  named categories were corrected, not silently exempted from the cap.
- `...112056701952Z.json`: missing local-term counter key.
- `...112105096375Z.json`: list-versus-scalar verifier comparison.
- `...112120213933Z.json`: both direct and column routes reversed the outputs
  of divmod(e,L), and thus AGREED on the wrong schedule. The independently
  reconstructed row caught the error. Its row error was not float noise.

The first corrected passing report `...112338829103Z.json` and stronger agent
report `...112554277299Z.json` remain available. Main subsequently removed a
redundant row recomputation, replaced a nested coefficient-label recovery by
direct quotient/remainder arithmetic, recorded explicit variation witnesses,
and increased the aggregate storage reserve for retained scalar lists and
old/new matrix views. The scientific fixture was not changed.

Final numeric preflight is **1,844,288 bytes**, below 16 MiB; Python container
headers and process RSS are not equated to numeric payload. Named actual
scalar/query terms total **505,824**, below the **804,096** preflight and
one-million cap. They include **18,432** column-local terms, **17,280**
row-local terms, **368,640** direct block-product terms, **19,104** reference
expansion entries and **55,584** pointwise modular powers. The **6,144** direct
vector shifts are counted separately; these are named arithmetic categories,
not a native instruction/bit-runtime total. Closed counts partition each
progression; the formula does not scan Q/r entries to determine class lengths.

`out/earlier_phase_cycles_claims.log` passes after adding exact residue-period,
truncated-partition and cancellation-revival regressions. The existing sampler
core/lab gates in WF cover the unchanged C78 helper. No further production
sampler change was made for C79. Documentation validation is
`out/earlier_phase_cycles_docs.log`. No manuscripts, abstracts, defaults,
host settings, firmware or commits changed.

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_earlier_phase_cycles
```

The skill's independent-reference and cancellation-control requirements
materially changed this audit. At this coefficient checkpoint, TODO 39's next
step was a complete-output comparison before sampler integration. That
comparison follows below; the coefficient result alone was not interpreted
as an observable advantage or hardness result.

## Complete-output checkpoint

Main retained the same eight-position fixture and froze complete-law
predictions before delegating three bounded tests. The cycle Fourier formula
uses closed geometric sums, the reference FFT uses literal columns, and the
selected Circuit/statevec check independently assembles those columns using
dense repeated blocks and pointwise phases. The latter embeds indexed labels
0..59 into a 64-state work register. It validates the QFT/output convention,
not a physical modular-arithmetic compiler.

Main's audited, reproduced reports all pass **4/4**:

- `out/earlier_phase_output_20260911T114909011094Z.json`, log
  `out/earlier_phase_output_main_final.log`: all eight complete laws from
  literal FFT and refined geometric components.
- `out/earlier_phase_feedback_20260911T115103487134Z.json`, log
  `out/earlier_phase_feedback_main_final.log`: dephase only the single high
  input history while retaining the early phase and exact Fourier feedback.
- `out/earlier_phase_qft_20260911T114642972207Z.json`, log
  `out/earlier_phase_qft_main_final.log`: existing Circuit/statevec references
  at v=0,2,3,7.

All use system Python **3.12.3**, NumPy **2.4.6**, one BLAS thread. The
geometric/literal FFT laws agree within **8.68e-18** per cell and
**1.28e-16** TV. The separately assembled Circuit-QFT laws agree with the
geometric route within **5.40e-16** TV. Literal columns normalize within
**6.67e-16**; target first-bit masses differ from one half by at most
**1.12e-16** in the all-position reference. C79/C77 own the exact support
argument; this is not an empirical claim of exactness from tiny errors.

The omission baseline removes ONLY the added early phase, yielding C78's
original schedule. It is computed once and reused. The feedback-aware
candidate retains the early phase and drops only coherence between the two
high input histories. Its zero-padded transform agrees with two feedback-aware
length-L conditional transforms within **3.47e-18** per output cell.

| earlier insertion v | omit early phase: full-law TV | dephase high history with feedback: TV |
|---|---:|---:|
| 0 | .02678914 | .01068526 |
| 1 | .04044824 | .01092837 |
| 2 | .08803871 | .01279596 |
| 3 | .25411263 | .03789582 |
| 4 | .28888720 | .02268610 |
| 5 | .32945628 | .02788206 |
| 6 | .31204759 | .01919770 |
| 7 | .00840446 | .00900860 |

Neither candidate reaches the frozen **1e-3** tolerance anywhere. Both reach
**1e-2** only at the final pre-mixer insertion. Retaining feedback and the
early phase improves seven of the eight comparisons, not all eight. This
stronger approximation explains much of the output variation, but leaves a
measurable coherent remainder. It does not establish that the remainder is
hard to compute, nor that constructing this candidate is generally cheap.

Main's closing baseline audit also identified the exact v=0 simplification
now stated in C79: absorb the early phase into the initial small unitary.
The existing C78 helper supports that transformed input without a new
schedule implementation. Thus the table compares TWO APPROXIMATIONS, not
the strongest exact baseline at every position. Their failure at v=0 is
not evidence that a new sampler is needed there. This endpoint reduction
is algebraic; no separate timing or new finite-RNG experiment is claimed.

Dropping coherent Fourier cross terms changes the full law by up to
**.239461** TV. Omitting the QFT or dephasing every exponent history gives
uniform-output controls that fail on the selected circuit laws; these are
closely related controls, not two independent algorithms. A synthetic flat
low-input state exercises the SAME feedback and missing-feedback functions
as the candidate. At nonzero prefix its conditional law is independently
checked by the geometric sum; missing feedback instead gives a delta, with
**.594695** TV. This replaced an initial control that merely duplicated the
standalone FFT expression.

## Full-output audit and interpretation

The initial output run's NumPy-Boolean serialization failure remains in
`out/earlier_phase_output_initial.log`. Intermediate output passing reports
are retained; their original threshold predicates and resource accounting
were weaker. Main/agent revisions added actual threshold classification,
cumulative counters, before-call bounds, two-sided marginal checks and
explicit storage for old/new views and converted scalar lists.

The feedback preflight failure is
`out/earlier_phase_feedback_20260911T114207329027Z.json`; subsequent intermediate
passes remain available. The final ledger removes fictitious FFT FLOP totals
and counts actual transform calls/input entries, local block products, shifts,
pointwise modular powers and chi entries. Main strengthened the synthetic
control and reconciled exact call/entry counts before the authoritative run.

Two aborted QFT reports, `...113905072355Z.json` and
`...113923793762Z.json`, retain failed check flags but no exception payload:
the first serializer saved summary rows instead of the full report. Their
causes cannot be reconstructed from those artifacts alone. The first passing
`...113933320015Z.json` likewise omitted the full-law evidence. The corrected
agent report `...114001016872Z.json` and main's final report retain all raw
probability vectors, counters and exceptions when present. A passing summary
was not accepted as a complete reference.

Final numeric payload reserves are **1,567,648**, **2,146,592**, and
**5,455,872 bytes** for output, feedback and QFT respectively, each below
16 MiB. These include conservative numeric temporaries, not measured RSS.
The output route performs **1,920,512** geometric-component terms against a
**4,147,200** preflight below the five-million cap; it also records **368,640**
block products, **9,720** row-local terms, **30,168** modular powers and
eight FFT calls. Feedback records **892,710** actual named units against
**895,788** preflight, including **44** FFT calls/**445,824** input entries
and **19,962** modular powers (the latter are also phase-query subcounts,
not double-added). QFT records **3,690,000** dense matvec terms and
**9,437,184** gate-entry updates, under separate five-/twelve-million caps.
These are diagnostic arithmetic categories, not end-to-end bit-runtime costs.

While examining the surviving two-history interference, main derived the
factor-two envelope now in C79; an independent audit checked ordering, Fourier
sign and the zero-mass case. This is the existing C59 Cauchy-Schwarz mechanism,
not a new sampling principle. It suggests another possible correction route,
but cheap proposal sampling and cheap coherence evaluation are still missing
premises. Do not turn a complete tiny probability table into a free sampler
oracle. Likewise, a uniform first output bit does not imply an independent
suffix or justify the dephasing approximation.

## A complementary cover, not a universally better representation

While the output comparison ran, main noticed that fixing the REMAINING low
history also makes the early phase constant. Main froze the dual construction
in TODO 39 before the pilot. Two independent proof audits checked aliases,
short intervals, cancellation handling and the balanced worst-position bound;
C79 owns the proof and its exponential-width limitation. This explains the
previous endpoint simplification as part of a broader two-sided partition.

Main's final dual report is
`out/earlier_phase_dual_20260911T120009418111Z.json`, **4/4 PASS**, with log
`out/earlier_phase_dual_main_final.log`. It reconstructs all literal amplitudes
and evaluates the dual FULL output law directly with existing geometric sums,
separately from FFTs of the expanded rows. Maximum dual amplitude error is
**2.78e-16**, row norm error **2.67e-15**, geometric/FFT TV **1.08e-16**, and
complete-law normalization error **1.12e-16**.

| insertion v | largest cycle-cover count | largest dual-cover count |
|---|---:|---:|
| 0,1,2 | 10 | 25 |
| 3 | 20 | 25 |
| 4,5 | 25 | 25 |
| 6 | 25 | 20 |
| 7 | 10 | 10 |

Neither cover dominates. The late interior insertion exhibits the predicted
benefit, while the early insertions favor the cycle cover. Most progressions
in this tiny fixture are short, so these counts do not demonstrate the full
potential gap or an implementation speedup. Dropping the remaining-history
shift from the phase fails at the fixed interior test v=3,j=0: the first
measured witness has h=0,k=7,a=2,e=58 and amplitude error **.102931**.
It was not the initially guessed k=1; the report retains the actual label.

Important retained failures:

- `earlier_phase_dual_20260911T115116525719Z.json`: a missing wrong-control
  expansion-counter key.
- `...115125183901Z.json`: a **misleading 4/4 PASS**, which main rejected.
  All three FFT routes used FFT/sqrt(Q) on unit-norm literal columns, so the
  reported laws summed to Q instead of one. Agreement alone missed the shared
  normalization defect. The first pilot also expanded every row and lacked
  the requested direct geometric-law check; it was not a sparse-formula test.
- `...115651670394Z.json`: after strengthening the verifier, it summed work
  rows rather than exponent columns for the column-norm audit. It also mixed
  overlapping geometric subcounts into a misleading aggregate cost check.

Corrected intermediate reports `...115824724024Z.json` and
`...115852379070Z.json` remain retained. Main's final reproduction adds explicit
structural component-bound checks, pointwise-query budgets and exact
reconciliation of local, storage, transform and geometric subcounts. All
laws are now explicitly finite, nonnegative and normalized. No probability
tolerance or circuit parameter was changed to make a check pass.

The final conservative numeric payload reserve is **7,951,360 bytes**, below
16 MiB; it deliberately includes oversized bounded scalar/container payload
reserves and is not measured RSS. The run uses **2,264,064** geometric terms
against a **3,686,400** preflight/five-million category cap, **24** FFT calls
with **368,640** input entries, **275,688** dual-local terms including the
wrong-phase control, and **18,068** dual phase/modular-power queries.
Geometric sum/phase calls are subcounts of the same component loop, not
three independent estimates of its runtime. Order/index discovery, precision
and efficient sampler oracles are not provided by this diagnostic.

## Final cross-check and scoped validation

`out/earlier_phase_output_cross_reference_final.log` compares all four final
reports directly: normalized output vectors, independent omission baselines,
parity interleaving, cycle/dual geometric laws and the proved factor-two
envelope. It checks complete finite laws, not Monte Carlo histograms or a few
selected probabilities. Sampled-grid ratio extrema are not promoted to a
general tighter rejection envelope.

Core passed first in `out/earlier_phase_output_core.log`. The affected
claims gate passes in `out/earlier_phase_output_claims_final.log`, including
exact support-separation, dual-partition and balanced-factor regressions.
Documentation is reindexed and checked in `out/earlier_phase_output_docs_final.log`.
No shared production helper changed; this turn does not claim a fresh lab
suite, full nine-suite run, numerical certificate or host-stability diagnosis.
No manuscripts, abstracts, defaults, host settings, firmware or commits changed.

Reproduce the new bounded experiments from research/:

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_earlier_phase_output
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_earlier_phase_feedback
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_earlier_phase_qft
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_earlier_phase_dual
```

The skill materially shaped the outcome: a stronger same-output baseline
reduced the apparent difficulty, and independent normalization/geometric
checks rejected a misleading consensus. The positive result is a proved
choice of representations with initial output validation, not a breakthrough
in factoring or a demonstrated broad simulation advantage. TODO 39 remains
OPEN for a costed sampler prototype and an untruncated comparison; the user's
open-ended research goal remains active.

## Sampler implementation and long-progression checkpoint

The preceding open implementation items are now complete at bounded FLOAT
scope. Main added EarlierPhaseProgressions to `lab/work_first.py`, sharing
the existing work draw, conditional rejection and progression kernel. C79
owns its schedule, choice-of-cover rule, counters and input/precision limits.
No default sampler or manuscript was changed. Three lower-cost initial testers
provided the independent experiments; main audited and reproduced them.

Final main reports, all **4/4 PASS**:

- Joint law: `out/earlier_phase_sampler_20260911T123806945928Z.json`,
  log `out/earlier_phase_sampler_main_final.log`.
- RNG/edges: `out/earlier_phase_sampler_edges_20260911T123729565062Z.json`,
  log `out/earlier_phase_sampler_edges_main_final.log`.
- Long rows: `out/earlier_phase_sampler_scaling_20260911T1234284124Z.json`,
  log `out/earlier_phase_sampler_scaling_main_final.log`. The historical
  scaling filename formatter omits seconds; use the exact path, not a guessed
  reconstruction of its timestamp.

The joint test retains the complete normalized work/output arrays for every
v=0..7 in the frozen r60,b3,t8,s7 fixture. Auto-law maximum entry error versus
literal-column FFT is **5.43e-19**. Explicit cycle/dual checks evaluate four
selected outputs for EVERY work row, with errors below **3.80e-19**; they
are not mislabeled complete explicit-cover laws. The exact v=0 folded-W0
C78 baseline agrees within **5.43e-19**. Archived omission/dephasing vectors
retain the earlier output-comparison provenance and are not supplied as free
proposal samplers. Dropping coherent acceptance gives per-insertion joint TV
between **.3737 and .4899**; the equal-mixture TV is **.287616**. Each law
is normalized separately before interpreting TV.

The **256** bounded seeded draws use **3,569** attempts, **59,584** coherent
Fourier-component terms and **34,662** interval-marginal queries. Every draw
keeps exactly one work selection. This is a path/counter check, not a histogram
accuracy estimate. This particular verifier aborts on sample exhaustion and
does not credit a failed call; the separate edge experiment captures failed
calls explicitly. Its complete-law and baseline work are separately charged:
the cover ledger records **321,144** row-local pair visits and **1,822,736**
Fourier-component terms, with the folded baseline charged separately.

### Actual RNG law and adversarial edges

For r3,b1,Q8,s2,v1 with an early parity phase, an independent literal reference
is phi_e=(-1)^e|e mod3>. Cycle rows have component counts (3,3,2), stride6,
two fair interval bits and two gcd lifts. The experiment enumerates all
**176** seed/component/reduced-output/lift combinations and executes the
positive accept/reject branches. Their path weights, not sample frequencies,
give accepted/rejected masses **.375/.625**. Restoring the original work
marginal after conditional rejection matches the independent full joint law
within **8.68e-19**. The actual scripted execution uses **332 sampler calls**,
**333 attempts**, **159 failed calls**, **1,664 real draws** and **665 integer
draws**. The extra attempt is a real partial-acceptance retry, not a second
work draw. An explicit finite-cap path raises without replacement.

The final edge verifier captures sampler counters even on exceptions,
including queried-invalid phases. Noncallable phases and structural caps are
tested with a matrix-copy trap. Zero rows, aliased cells and cancellation
revival are independently checked. A genuine count>1 stride control first
reconstructs the correct law (TV below **6.64e-16**), then substitutes stride3
for stride6. Spurious collisions change the wrong state's norm to **.5**;
after explicitly normalizing that DIFFERENT state, its TV is **.424228**.
The report also retains the unnormalized half-L1 difference, without calling
it a probability-law TV. The permanent lab regression separately checks the
weighted actual transition and detects the work-marginal bias from restarting
work after rejection.

### Long rows: a real cover gap, not a general simulator comparison

Only split s varies, in unit steps, under the frozen schedule t=s+1,v=s-1;
H=K=2 and both phases/matrices remain fixed. The reference calls the existing
literal column engine only on the union of nine support residue classes for
work0,1,59. The proven offset cone excludes every omitted entry. It uses a
Q-by-3 array, never Q-by-r. Both covers reconstruct amplitudes and work norms;
complete selected conditional FFT laws and four geometric output queries
are compared independently. Largest amplitude error is **2.49e-16**,
conditional-law TV **1.87e-15**, and selected geometric probability error
**6.60e-17**.

| split s | largest tested cycle component count | largest tested dual count | longest tested dual progression |
|---|---:|---:|---:|
| 7 | 23 | 20 | 2 |
| 8 | 40 | 20 | 3 |
| 9 | 88 | 20 | 5 |
| 10 | 170 | 20 | 9 |
| 11 | 343 | 20 | 18 |
| 12 | 680 | 20 | 35 |

These are maxima over three selected work labels, not all-work maxima.
Auto selects dual throughout. The independently proved all-row dual bound
is the one stated in C79, not inferred from this table. The reference performs
**2,427** column calls, **436,860** dense block-product terms and **18,597**
pointwise modular powers. All helper rows, including folded/omission controls
and sample rebuilds, charge **16,893** row-pair visits and **12,418** inclusive
phase queries. The **192** actual auto samples use **3,367** attempts and
**62,896** coherent-component terms. Fixed-r sequential simulation is ALSO
width-linear: this experiment establishes a choice-of-cover advantage, not
an end-to-end win over that strongest simple baseline. No timing comparison
or whole-large-joint-law test was performed.

### Verifier corrections retained, including misleading passes

Main did not accept the initial unanimous checks as sufficient evidence:

- Joint `...122053149462Z.json` failed from an output/work-axis reversal.
  The following `...122112525970Z.json` passed but used concatenated laws
  for a misleading TV; subsequent reports normalize the comparison per law
  or as an explicitly equal mixture. `...122824506835Z.json` then rejected
  an underestimated cumulative explicit-dual row budget. Both diagnostic
  covers, selected public calls, setup and marginal queries are now charged.
- Edge `...122006951089Z.json` expected arbitrary callable values to fail at
  construction, contradicting the documented query-time contract. A later
  `...122814529039Z.json` passed despite an extra row-norm divisor in its
  wrong-stride control and incomplete failed-call accounting. Even
  `...123339273177Z.json` computed totals before the final controls and called
  a nonnormalized mass difference TV. Final main reports correct these
  predicates, include zero-acceptance failures, and retain reference work.
- Scaling `...1224744346Z.json` rejected an overestimated component-times-Q
  expansion budget. `...1225185158Z.json` failed law/baseline checks; the
  source audit separately identified inconsistent split/width meanings and
  a reference phase inserted at the endpoint instead of one control earlier.
  `...1225801610Z.json` attempted a sampled-work lookup in a cache containing
  only three selected labels. Later passing agent reports still omitted
  baseline/control row work; main added cumulative reservations and retained
  actual per-call counters before the final reproduction.

These filenames use their respective full experiment prefixes above. The
failed/intermediate JSONs and logs remain in out/. Additional import/invocation
failures are retained as logs, not invented as numeric reports. Neither the
science fixture nor tolerance was retuned to turn a failure into a pass.

Final aggregate numeric reserves are **13,427,896 / 280,832 / 4,698,624 bytes**
for joint/edge/scaling, each below16MiB. They include bounded retained arrays,
numeric tuples/lists and live temporaries; the joint reserve also accounts
for parsed archived sources. These are not Python/RSS measurements, nor
end-to-end precision costs. Caps precede the relevant calls/loops. Inclusive
phase totals and early/late subcounts are not added as independent costs;
geometric component counts likewise are not native FLOP totals.

Core passed first in `out/earlier_sampler_core.log`. The affected lab suite
passes with optional python-flint in `out/earlier_sampler_lab_final.log` and
without it in `out/earlier_sampler_lab_no_backend_final.log` (optional checks
explicitly skipped there). Claims pass in `out/earlier_sampler_claims_final.log`.
The other six science suites were not rerun. Reindex/documentation validation
is recorded in `out/earlier_sampler_docs_final.log`. The small runs do not
resolve TODO34's host-instability issue.

Reproduce the new bounded experiments from research/:

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_earlier_phase_sampler
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_earlier_phase_sampler_edges
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_earlier_phase_sampler_scaling
```

The skill's actual-output, nonvacuous-control and strongest-baseline rules
materially shaped this checkpoint. A follow-up independent read-only algebra
audit checked C59's mass-weighted envelope for these components; TODO39 owns
its frozen next pilot. It has not changed the implemented proposal. The
current result is a restricted, costed float sampler plus a genuine long-row
representation gap, not a generic breakthrough. The broader user goal remains
active; no manuscripts, abstracts, host settings, firmware or commits changed.
