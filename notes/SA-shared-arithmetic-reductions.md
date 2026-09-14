---
code: SA
date: 2026-09-12
title: "Shared arithmetic has exact scalar reductions; the raw cubic certificate redirects before execution"
outcome: two-macro full average implemented and checked; clean rotation and raw-lift limits proved; growing dirty prefixes remain open
claims: [C90, C91]
todo: [50]
---
# SA — Shared-variable reductions and a rejected raw encoding

C90 owns the causal-lift theorem and source baselines. C91 owns the exact
arithmetic reductions, mask scope and clean/dirty distinction. This is progress
within TODO50, which remains open; no general simulation breakthrough or
best-known memory advantage is established. The broad user goal remains active.

## How the two lines of investigation changed the next action

The side-board proposals in `topic:cnot-math-brainstorm` suggested global
quadratic restriction and parity-rank decompositions. Sol read the cited
primary theorem bodies. Root independently specialized the degree-drop kernel
to the factored causal lift, and a second Sol worker checked the native
control planes and affine-branch converse. The first carry stages already
settle rejection of the raw representation, while its successful one-parity
cases have an affine baseline. The worker also derived a nonzero easy endpoint
that the raw lift rejects. No native size experiment was warranted.

Root had created task `T2e51fb2255364667`, attempt `Ad603d1740e204919`, and
registered prediction `Me2e8a7915464497a` for a small recognizer comparison.
The independent proof settled its proposed discriminator before any source
was written or run started. The task was cancelled with that reason. The
prediction is unexecuted, neither confirmed nor refuted by measurement.

Meanwhile Astra derived the two-macro autocorrelation and its constant-state
full outside average. Root independently derived a nonzero full-space probe;
Astra checked it before execution. Its order independence exposed a possible
vacuous chronology test, so Astra supplied a different conditional query that
does detect order. Root implemented both comparisons against the actual gate
trace. The full outside average and the fixed-dirty-fiber chronology query
have deliberately different masks and normalization contracts.

Astra additionally proved the three-domain clean-scratch conjugacy and its
native a=1/a=2 scalar reduction. Root checked the branch equations, bijection
and bounded carry/borrow construction. These latter results are proof-only:
there is no implemented chart decoder, clean scalar API or benchmark.

## Actual-gate experiment

Source: `experiments/experiment_shared_control_average.py`. Reproduce with:

```bash
timeout 60s uv run python -u -X faulthandler -m experiments.experiment_shared_control_average
```

Prediction `M7822d1cb9ef14f7d` preceded the original run
`Rd56eed62abd84e80`. All eight scientific checks passed with terminal exit
zero, but the source called `finish()` without `report_path` and put rows and
metadata on unused attributes. Its two control results also used the generic
check method. The numerical comparisons were sound; the JSON report was
absent and the control labels needed correction.

Original source and raw log are preserved as
`out/shared_control_average_revision1.py` and
`out/shared_control_average_revision1.log`. Root corrected only report
arguments and control reporting. Prediction `Mc07d35ef2b194ad4` preceded
the identical-fixture rerun `R339f6c762a054d2d`, which again terminated with
exit zero. The corrected report is `out/shared_control_average_report.json`,
with raw log `out/shared_control_average_main.log`. All eight checks pass
without warnings, including both explicitly marked must-fail controls.

| exact comparison | scope | result |
|---|---|---|
| full shared-variable DP vs actual gate sums | N=3,a=2, two native macros, all 128 combined outside masks | zero errors; 18 nonzero coefficients |
| conditional low-bit autocorrelation vs actual gate sums | all 64 t,h,e1,e2 assignments in that fixture | zero errors |
| order-independent low-bit family vs reversed word map | same 64 assignments | zero errors |
| signed carry DP vs explicit t,h,w integer sums | N=9,c=(3,6),m=5; 24 mask/enable cases | zero errors; negative carry coefficients exercised |
| order-sensitive conditional query vs actual gates | N=3,a=1,t=2,h=0,u=1, both macro orders | C91's two distinct predicted values agree |

The first fixture has eleven qubits and 216 actual logical operations. Each
full-space reference coefficient sums all 2,048 labels. Three such small
gate maps are constructed across the two fixtures; the larger-word arithmetic
check constructs no seventeen-qubit gate map. Missing incoming carry destroys
C91's nonzero full-space probe; treating disabled dirty macros as identity
changes the conditional chronological probe as predicted. All constants,
masks, controls and caps were frozen before execution. These are correctness
fixtures, not a width-scaling or timing experiment.

The direct side uses existing `walsh.classical_permutation` replay of logical
X/CNOT/Toffoli operations. The carry DP has a different scalar construction.
C89's word map is used only for explicit wrong/reversed references; it is not
the main full-space reference. The test does not independently verify the
stored Clifford+T decompositions as unitaries. Exact integers/Fractions avoid
floating-point tolerances. Allocated bytes, RSS and comparative runtime were
not measured; the harness elapsed time is descriptive only.

## Reviewed evidence and retained corrections

| board task | accepted submission | accepting review |
|---|---|---|
| Td63294941a95401f, Astra arithmetic proof | Sa0255906e6c04810 | Veeff7565b3ec45f2 |
| T8f87ef33ee7f48c4, Sol primary-source audit | S73d41f041b4742c3 | V88cd6475777446c4 |
| T6260d3aeb7394589, Sol native-lift proof | S073d9402a64342bc | V9305195d5915428f |
| Ta65a555030304995, root implementation | S0721470decb1493a | V45550e55b4e846ab |

Proof artifact: `out/agent-board/workers/Afccf9c8c16cc42d9/signed_enable_reduction.md`.
Source audit: `out/agent-board/workers/A00093affbfac4b39/cubic_degree_drop_rankwidth_audit.md`.
Native-lift audit: `out/agent-board/workers/A63afbd07db3341ef/audit.md`.
Root read the full actual artifacts and checked archived hashes/current
identities. Astra independently read both implementation versions, the exact
reporting diff, both logs, corrected report and core evidence before accepting.
Workers made no canonical edits.

The native-lift worker's initial submission `Sad8dfbdf3bd84024` had corrupted
summary/limitation text despite a valid evidence artifact. Review
`V8c7fbd4520124948` requested clean metadata and qualification of an endpoint
Gauss-sum argument that had implied unverified C83 intermediate API closure.
The accepted revision corrects both; the original submission is retained.
No scientific execution failed in this investigation.

The primary audit includes Carlet's Boolean multiplication kernel,
Montanaro's supplied hitting-set/linear-change algorithms and de Colnet's
exact quadratic sign-edge rank-width baseline. Root read the relevant primary
definitions and algorithm/theorem bodies before using those comparisons.
The rank decomposition is supplied and its width belongs to the expanded
representation; ordinary Gauss elimination already handles Boolean quadratic
slices irrespective of that width. No favorable native decomposition has been
constructed here. C90 owns the corresponding source links and claim limits.

## Validation scope

Core ran before new science and after the reporting correction; both observed
terminal exits were zero and their logs end in ALL TESTS PASSED:
`out/shared_average_core_initial.log` and `out/shared_average_core_final.log`.
The affected new experiment passed as above. Production helpers were unchanged
and the other eight science suites were not rerun. There was no timeout or
native crash in these bounded runs; TODO34 retains its prior unresolved issue.
Documentation indexing exits zero (`out/shared_average_reindex.log`), and
all ten documentation checks pass (`out/shared_average_docs_check.log`).
Scoped whitespace and new-source syntax checks pass. All four reviewed tasks
are accepted and closed; the redundant raw-lift pilot is cancelled. C90/C91
canonical transfers received independent review with no further correction.
No manuscript, dependency, default or host setting changed; nothing was
committed or published.
