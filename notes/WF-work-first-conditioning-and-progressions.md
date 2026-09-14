---
code: WF
date: 2026-09-11
title: "Work-first conditioning trades dense sectors for sparse rows or a few coherent progressions"
outcome: confirmed
claims: [C54, C55, C57, C76, C77, C78]
todo: [14, 38, 39]
---
# WF — Sample a work outcome, then restore all compatible exponent histories

C78 owns the proof and precise oracle/schedule/cost contracts. This note owns
the discovery and bounded implementation evidence. The starting observation
was FB's residual high-history coherence, not a new dense-sector count.

Main first considered support-based bounds on the cross-history terms. A more
useful change was the contraction order: a final computational WORK measurement
commutes with the exponent QFT. Sampling its marginal is easy given a coherent
fixed-history column. Conditioning on the work index then preserves only the
compatible exponent histories, whose amplitudes must be reconstructed coherently.
The auxiliary exponent used to draw the work marginal must be discarded.

Initially the finite cone suggested a sparse row whose size includes Q/r.
C77's partial-Fourier identity then supplied a rejection-free grouped-prefix
sampler for that row. Main subsequently noticed that the literal late-phase
schedule has stronger structure: its long rows are unions of few arithmetic
progressions with constant coefficients. Existing C55 progression sampling
and coherent-component rejection apply. Their relative phases and the actual
work weights are essential; they are not a classical mixture approximation.

Two lower-cost audits independently checked the sparse-row identity and the
general late-only progression proof. The component rejection uses a mixture
of component FOURIER laws, not a uniform-output proposal; confusing those
proposals would invalidate the exact mean-attempt statement. C78 explicitly
separates the two methods. Known order/index, pointwise phase access, coherent
column evaluation and a small number of late HISTORIES are substantive inputs.

## Initial sparse-row pilot

The frozen row pilot is the same UP circuit, with only its common terminal W
stripped under the work trace. No physical gate is commuted across a control.
The final main report is
`out/work_first_rows_20260911T104128015399Z.json` (3/3 PASS), with raw output
in `out/work_first_rows_main.log`. It uses the existing dense closed boundary
as a CHARGED tiny oracle, not a scalable column implementation.

The complete row-mixture law agrees with the independent UP reference at
TV **2.41e-16**. Modular progression enumeration covers all nonzero row
amplitudes, and the resulting work marginal agrees with the direct column
Born marginal within **6.94e-18** TV. The work distribution is not uniform:
its TV from uniform is **.06944914**. Incorrect uniform-work weighting changes
the output by **.00463295** TV.

The must-fail control explicitly builds the conditional Born mixture by
discarding exponent-row interference. It yields the uniform output, differing
from the coherent law by **.23462695** TV. This tests the load-bearing
coherences, not just normalization of a positive distribution. Every tested
prefix agrees with the independently grouped reference and normalizes.

There are **576** candidate entries across sixty work rows, within the
conservative **1,080** preflight. Main found that the first implementation
scanned every empty residue group despite claiming a sparse grouping method,
and reported planned upper bounds as actual terms. The corrected implementation
groups only occupied residues, increments actual counters inside the loop and
guards cumulative budgets before calculation. Final terms are **8,064** for
the requested short prefixes, **36,864** for complete row laws and **36,864**
for the Born-mixture control. The dense oracle separately charges sixty-four
columns and **460,800** matvec terms, with column norm error below **4.5e-16**.

The final aggregate numeric preflight is **592,896 bytes**, including retained
final/pre-final rows, construction matrices/temporaries and laws, under 16 MiB.
It is not a process RSS bound. Reference support scans and dense oracle work
remain part of this tiny verifier, not costs removed by an implemented sampler.

An initial orbit-coordinate/physical-label indexing mistake was corrected
without changing the circuit. A serialization failure also occurred. The
file `out/work_first_rows_failure_serialization.log` is a RETROSPECTIVE agent
description, not its original raw traceback; do not cite it as a retained
raw run. The first retained passing report is
`out/work_first_rows_20260911T103444852474Z.json`, followed by the strengthened
agent report `out/work_first_rows_20260911T103811510460Z.json`. Main's final
run above uses explicit system Python 3.12.3/NumPy 2.4.6 and one BLAS thread.
The agent confirmed that the original indexing traceback was not retained;
the later passing runs reproduce the corrected law, but do not repair that
provenance gap. Retained agent reports also identify system Python 3.12.3.

## Progression extension and limits

For the literal sparse late-phase schedule, main derived the finite component
formula before requesting another lower-cost initial probe. C78 owns the
general b/r/H proof and its upper bound. It does not require Gauss-sum phases:
the diagonal phase is evaluated at a particular work index after conditioning.
This is why dense coupling in the old Fourier-sector basis need not govern
the full sampling task in this restricted schedule.

The final main progression report is
`out/work_first_progressions_20260911T105444407294Z.json`, **4/4 PASS**, in
`out/work_first_progressions_main.log`. Its primary literal sparse columns
and progression rows agree with the charged existing pre-final oracle within
**7.11e-16** amplitude error. The complete coherent law agrees with UP within
**2.16e-16** TV. The primary maximum is nine components; its envelope and
accepted-law checks match the proved row-dependent mean attempt count.

Main found that the initial prototype READ an orbit table even though it
allocated no such table inside its helpers. The corrected literal column/row
constructors evaluate physical phases with modular exponentiation at the
queried index. Only the tiny dense reference uses the supplied orbit list.
The fixed six-bit fixture also had singleton progressions throughout, so it
did not test the long-progression identity. Main froze a separate width-eight
row test, retaining the same N,a,r,b,angles,k and penultimate-control-relative
insertion rule. This is an explicit width family, not a retuned UP law.

All secondary rows now contain two- or three-term progressions. Their expanded
coefficients agree with literal sparse columns within **1.58e-16**; geometric
Fourier sums agree with direct finite-root sums within **2.78e-17**. Forty
work rows have nine components and twenty have ten, reaching the C78 upper
bound. Main also checks the genuinely NONUNIFORM component-mixture proposal
on these longer components: maximum acceptance is **.43649684**, error in
its mean acceptance 1/m is below **5.6e-17**, and its accepted law matches
the coherent law within **7e-18** per outcome. The secondary route is a
row/formula check, not an independent compiled arithmetic/QFT reference.

Counters are instrumented before loop operations. Both fixtures together
charge **2,880** column-local terms, **2,160** row-local terms, **1,492**
reference expansion entries, **162,432** geometric-component evaluations
and **305,664** direct-root terms. The combined preflight is **660,240**
of these named terms under a one-million cap. Pointwise phase evaluations
separately charge **2,800** modular-power queries, plus sixty reference label
queries. The dense primary oracle's work is separately checked. The final
numeric preflight is **702,848 bytes**, including retained secondary sparse
column values and temporary law arrays; Python object/report overhead is
bounded in object count, not asserted to equal this numeric payload or RSS.

The first progression report
`out/work_first_progressions_20260911T104249397048Z.json` preserves a float-
argument error at the integer geometric-sum API. Passing primary reports
`out/work_first_progressions_20260911T104309022089Z.json` and
`out/work_first_progressions_20260911T104345423424Z.json` retain the earlier
table-dependent, singleton-only implementation. A direct-script invocation
failed to import lab in `out/work_first_progressions_final_system312.log`;
reproduction uses the module command below. Report
`out/work_first_progressions_20260911T105127657345Z.json` preserves the
secondary control failure from leaking primary L=32 into the L=128 component
start. That parameter leak was fixed without changing either fixture.
The agent's subsequent reports pass but identify Python **3.12.10**, despite
their log filenames saying system312. Main's final report above records actual
system Python **3.12.3**, NumPy **2.4.6** and the explicitly single-thread run.
Neither set of small passes settles the separate host-reliability question.

Main removed a 1e-30 proposal cutoff, added explicit zero-proposal handling,
strengthened setup/work and law predicates, charged omitted secondary data,
and checked the long-row acceptance law without extra propagation. No cutoff
prunes a small positive component. Negative controls retain both the
incoherent-proposal and incorrect-uniform-work output discrepancies.
The lower-cost agent then read main's final source/report and found no
substantive remaining correctness or scope defect. These named operation
counters are not a native instruction count or an arbitrary-precision bound.

At the formula-only checkpoint, the existing interval/progression primitive
supplied the next implementation, but those diagnostics were not a working
end-to-end sampler. No random samples or timing comparison had been performed.
The subsequent implementation is recorded below. Small late-history count, gate setup,
order/index discovery, precision and failure handling must remain explicit.
General multi-defect schedules can destroy the constant-progression coefficients.
No general quantum-simulation/factoring advance or timing advantage follows.

Main read Van den Nest's CT definition and sparse-operator theorem body and
Schwarz–Van den Nest's sparse-output theorem assumptions; C78 links the exact
primary sources and distinguishes their contracts. These are standard
conditioning/Fourier/rejection ingredients. A complete priority study for the
specialization is still open; no novelty claim is made.

## Formula-checkpoint validation

At that checkpoint, the C77 core/lab/backend checks were the production gate.
New exact Gaussian-integer regressions in test_claims.py check work weights,
zero rows, sparse prefixes, full interference and component acceptance;
`out/work_first_claims.log` passes. The regenerated documentation gate is
`out/work_first_docs.log`. No production sampler, manuscript, abstract,
default, host setting or firmware was changed, and no commit was made.

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_work_first_rows
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_work_first_progressions
```

The qsim-research skill's actual-observable and coherent-reference rules led
to the work-first route and exposed the verifier's hidden dense loops. This
is the most promising restricted sampler connection in the present batch,
not a declared general breakthrough. TODO 38 subsequently completed the bounded
implementation and strongest-baseline comparison below; the broader goal stays active.

## Implemented sampler and actual decision paths

Main added `lab/work_first.py:LateWorkProgressions`, with supplied small
unitaries, known/indexed period and a deterministic pointwise phase callable.
Its scalar rejection loop reuses `progression_sample`; it does not allocate
orbit/output/exponent arrays or expand progression entries. Input/cost/precision
limits remain C78's contract. Queried phases and supplied unitaries are checked
at float tolerances, not globally certified. Exhaustion raises without a
fallback; failed calls must not simply be discarded.

Three lower-cost initial testers supplied the scientific/edge/comparison
experiments. Main's final reports are:

| test | report under `out/` | verdict |
|---|---|---|
| joint formulas and actual draws | `work_first_sampler_20260911T111839577336Z.json` | 5/5 |
| zero/cap/retry edges | `work_first_sampler_edges_20260911T111528371942Z.json` | 4/4 |
| strongest same-output comparison | `work_first_sampler_comparison_20260911T111839083292Z.json` | 5/5 |

All three main reports identify system Python **3.12.3**, NumPy **2.4.6**,
with one BLAS thread. Main logs are respectively
`work_first_sampler_main_final.log`, `work_first_sampler_edges_main.log`,
and `work_first_sampler_comparison_main_final.log`. Lower-cost sampler pilots
identify Python 3.12.10; edge/comparison pilots identify 3.12.3. A filename
or requested interpreter does not supersede the recorded runtime.

The sampler's primary complete JOINT law matches the dense pre-final oracle
within **1.09e-18** per entry. The independently checked work marginal differs
by **9.54e-17** TV. All conditional/proposal laws normalize; the normalized
accepted law matches the coherent law within **6.94e-18**. Both width fixtures
retain the former formula checks, including nontrivial progressions at width
eight. There are **512 actual draws per width**, with means **5.5215** and
**9.3359** proposals and maxima **48** and **72**. Each draw checks ranges,
positive target mass and actual local/phase/progression/prefix counters.
The **7,607** total proposals are below the explicit **262,144** worst-case
budget, with each individual call capped at 256. These draws are smoke tests,
not histogram evidence of a one-percent accuracy guarantee.

Main's independent `test_lab.py` FFT regression checks a small-period case
where residues alias. A second regression enumerates **all 128 actual
first-attempt RNG paths** for r=b=2,Q=4,W0=I,W1=H, including the primitive's
interval bit and gcd lift. The accepted joint submass is exactly 1/8 at
(j,y)=(0,0),(1,2), with total rejected mass 3/4. A deliberate no-acceptance
control gives the wrong joint law. This exercises actual sampling decisions,
not merely the forced probability formula. The independent edge experiment
adds one real reject-then-accept path with varying row component counts,
explicit exhaustion, exact-zero cancellation, zero work/proposal masses and
pre-copy allocation traps. It is not itself an exhaustive multi-retry tree.

The same-output comparison reproduces the existing 60-dimensional sequential
reference and obtains TV **2.30e-16** for the new helper's enumerated float law.
The cached, independently validated FB candidate has TV **.01252649** and
omitting G has TV **.03459335** against that reference; neither meets the
frozen 1e-3 or 1e-2 law thresholds. The comparison additionally returns eight
helper samples (42 proposals) and four reference samples. This establishes
the restricted exact-arithmetic route's relevance beyond those approximations,
NOT certified finite-RNG accuracy or speed superiority. FB is explicitly a
cached baseline law in this experiment, not an uncharged fresh sampler run.

The incorrect redraw-work-after-rejection law differs by **.00666690** TV;
fully dephasing exponent histories differs by **.23462695**. The former is
essential because rejection success varies with the work row. Initial edge
code had a vacuous control that repeated a positive counter assertion; it now
executes a deliberately wrong two-column redraw path and checks independent
RNG call counts as well. Main strengthened the exact-zero and allocation tests.

## Resource and failure audit for the implementation

The helper's frozen-fixture reserve is **6,272 numeric bytes**, including
conservative temporary storage, with **288 bytes** of retained matrices;
Python containers and the phase oracle's internals are separate. Its local
construction bound is 27 terms per draw, not a dense orbit scan. The full
sampler verifier reserves **1,053,072 numeric bytes** and charges **824,156**
actual named diagnostic scalar terms against **874,680** preflighted, plus
the separately budgeted dense oracle/setup and random proposals. The dense
Fourier reference charges all **245,760** root terms, the formula checks
**2,920** pointwise modular powers, and dense oracle matvec terms are retained
separately. The comparison's reserve is **3,846,784 bytes**, including the
existing sequential workspace. It charges **128** forced reference calls,
**four** reference draws, **6,048,000** dense setup product terms and **120**
pointwise phase queries. Counts of shared matrix references are not claimed
as fresh matrix allocations. The edge reserve is **56,064 bytes**, with
**64** actual local-column terms and **1,536** direct-root terms. None of
these numeric-payload bounds is a process RSS or bit-runtime measurement.

Preserved failures and superseded passing reports are not erased:

- `work_first_sampler_initial.log` retains a complex-JSON serialization
  traceback. `work_first_sampler_20260911T110724980686Z.json` failed the first
  zero-row control; the fixture was corrected to contain an unreachable row.
  `work_first_sampler_20260911T111413157986Z.json` retains a counter-predicate
  failure. Earlier passing pilots lacked the later complete accepted-law,
  resource and per-draw checks.
- `work_first_sampler_edges_20260911T110951977188Z.json` failed a tuple-versus-
  length zero-component predicate; `...111310306969Z.json` failed the scripted
  RNG call-count expectation. The corrected count accounts for the reduced
  Fourier dimension after taking the gcd. Main additionally replaced a
  planned direct-root count by a counter updated before actual evaluation.
- `work_first_sampler_comparison_20260911T111017335926Z.json` preserves the
  FB prefix-dimension mix-up (H=4 versus the sampler's late-history H=2).
  `...111037256587Z.json` preserves the wrong imported constructor's period
  (12 versus 60). These descriptions follow the actual retained tracebacks,
  not an agent's later reversed recollection. Main corrected shared-reference
  allocation accounting and strengthened returned-reference probability checks.

Main's final sampler verifier also fixes a leaked width-eight range bound
in the width-six sample check. No scientific fixture was retuned to repair
these verification defects, except replacing the invalid zero-row CONTROL
fixture. The full target schedules and fixed random seeds stayed unchanged.

## Implementation validation and next question

Core, lab (including the optional verified backend), and claims pass in
`out/work_first_sampler_core.log`, `out/work_first_sampler_lab_transition.log`,
and `out/work_first_sampler_claims.log`. The other six science suites were not
rerun for this opt-in helper. Documentation validation is recorded in
`out/work_first_sampler_docs.log`. No existing sampler default, manuscript,
abstract, host setting, firmware or commit changed. These small passes do not
resolve TODO 34's host issue.

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_work_first_sampler
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_work_first_sampler_edges
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_work_first_sampler_comparison
```

The qsim-research skill's independent-control, actual-output and charged-cost
requirements materially strengthened this implementation audit. TODO 38 is
complete at this bounded scope; TODO 39 owns the next, genuinely changed
earlier-phase schedule and its possible binary-residue-cycle connection.
The broader research goal stays active.
